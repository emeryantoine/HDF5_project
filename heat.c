#include <mpi.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//include the hdf5 library
#include <hdf5.h>
#define RANK 2

//my custom API for the HDF5 library
hid_t create_file(char* filename)
{
	return H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

hid_t create_space(int x, int y)
{
	hsize_t dimspace[] = {x, y};
	return H5Screate_simple (RANK, dimspace, NULL); 
}

hid_t create_dataset(hid_t file, char* step, hid_t dataspace)
{
	return H5Dcreate (file, step, H5T_STD_I32BE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}


/** A function to initialize the temperature at t=0
 * @param	  dsize  size of the local data block (including ghost zones)
 * @param	  pcoord position of the local data block in the array of data blocks
 * @param[out] dat	the local data block to initialize
 */
void init(int dsize[2], int pcoord[2], double dat[dsize[0]][dsize[1]])
{
	// initialize everything to 0
	for (int yy=0; yy<dsize[0]; ++yy) {
		for (int xx=0; xx<dsize[1]; ++xx) {
			dat[yy][xx] = 0;
		}
	}
	// except the boundary condition at x=0 if our block is at the boundary itself
	if ( pcoord[1] == 0 ) {
		for (int yy=0; yy<dsize[0]; ++yy) {
			dat[yy][0] = 1000000;
		}
	}
}

/** A function to compute the temperature at t+delta_t based on the temperature at t
 * @param	  dsize  size of the local data block (including ghost zones)
 * @param	  pcoord position of the local data block in the array of data blocks
 * @param[in]  cur	the current value (t) of the local data block
 * @param[out] next   the next value (t+delta_t) of the local data block
 */
void iter(int dsize[2], double cur[dsize[0]][dsize[1]], double next[dsize[0]][dsize[1]])
{
	int xx, yy;
	// copy the boundary values at x=0 (Dirichlet boundary condition)
	for (xx=0; xx<dsize[1]; ++xx) {
		next[0][xx] = cur[0][xx];
	}
	for (yy=1; yy<dsize[0]-1; ++yy) {
		// copy the boundary values at y=0 (Dirichlet boundary condition)
		next[yy][0] = cur[yy][0];
		for (xx=1; xx<dsize[1]-1; ++xx) {
			next[yy][xx] =
				(cur[yy][xx]   *.5)
				+ (cur[yy][xx-1] *.125)
				+ (cur[yy][xx+1] *.125)
				+ (cur[yy-1][xx] *.125)
				+ (cur[yy+1][xx] *.125);
		}
		// copy the boundary values at y=YMAX (Dirichlet boundary condition)
		next[yy][dsize[1]-1] = cur[yy][dsize[1]-1];
	}
	// copy the boundary values at x=XMAX (Dirichlet boundary condition)
	for (xx=0; xx<dsize[1]; ++xx) {
		next[dsize[0]-1][xx] = cur[dsize[0]-1][xx];
	}
}

/** A function to update the values of the ghost zones
 * @param	  cart_comm a MPI Cartesian communicator including all processes arranged in grid
 * @param	  dsize	 size of the local data block (including ghost zones)
 * @param[out] next	  the next value (t+delta_t) of the local data block
 */
void exchange(MPI_Comm cart_comm, int dsize[2], double cur[dsize[0]][dsize[1]])
{
	MPI_Status status;
	int rank_source, rank_dest;
	static MPI_Datatype column, row;
	static int initialized = 0;
	
	// Build the MPI datatypes if this is the first time this function is called
	if ( !initialized ) {
		// A vector column when exchanging width neighbours on left/right
		MPI_Type_vector(dsize[0]-2, 1, dsize[1], MPI_DOUBLE, &column);
		MPI_Type_commit(&column);
		// A row column when exchanging width neighbours on top/down
		MPI_Type_contiguous(dsize[1]-2, MPI_DOUBLE, &row);
		MPI_Type_commit(&row);
		initialized = 1;
	}
	
	// send to the bottom, receive from the top
	MPI_Cart_shift(cart_comm, 0, 1, &rank_source, &rank_dest);
	MPI_Sendrecv(&cur[dsize[0]-2][1], 1, row, rank_dest,   100, /* send row before ghost */
		&cur[0][1], 1, row, rank_source, 100, /* receive 1st row (ghost) */
		cart_comm, &status);
		
	// send to the top, receive from the bottom
	MPI_Cart_shift(cart_comm, 0, -1, &rank_source, &rank_dest);
	MPI_Sendrecv(&cur[1][1], 1, row, rank_dest,   100, /* send column after ghost */
		&cur[dsize[0]-1][1], 1, row, rank_source, 100, /* receive last column (ghost) */
		cart_comm, &status);
		
	// send to the right, receive from the left
	MPI_Cart_shift(cart_comm, 1, 1, &rank_source, &rank_dest);
	MPI_Sendrecv(&cur[1][dsize[1]-2], 1, column, rank_dest,   100, /* send column before ghost */
		&cur[1][0], 1, column, rank_source, 100, /* receive 1st column (ghost) */
		cart_comm, &status);
	
	// send to the left, receive from the right
	MPI_Cart_shift(cart_comm, 1, -1, &rank_source, &rank_dest);
	MPI_Sendrecv(&cur[1][1], 1, column, rank_dest,   100, /* send column after ghost */
		&cur[1][dsize[1]-1], 1, column, rank_source, 100, /* receive last column (ghost) */
		cart_comm, &status);
}

/** A function to parse command line arguments
 * @param	  argc	  number of arguments received on the command line
 * @param[in]  argv	  values of arguments received on the command line
 * @param[out] nb_iter   number of iterations to execute
 * @param[out] dsize	 size of the local data block (including ghost zones)
 * @param[out] cart_comm a MPI Cartesian communicator including all processes arranged in grid
 */
void parse_args( int argc, char *argv[], int *nb_iter, int dsize[2], MPI_Comm *cart_comm )
{
	if ( argc != 4 ) {
		printf("Usage: %s <Nb_iter> <height> <width>\n", argv[0]);
		exit(1);
	}
	
	// total number of processes
	int comm_size; MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	int psize[2];
	
	// number of processes in the x dimension
	psize[0] = sqrt(comm_size);
	if ( psize[0]<1 ) psize[0]=1; 
	// number of processes in the y dimension
	psize[1] = comm_size/psize[0];
	// make sure the total number of processes is correct
	if (psize[0]*psize[1] != comm_size) {
		fprintf(stderr, "Error: invalid number of processes\n");
		abort();
	}
	
	// number of iterations
	*nb_iter = atoi(argv[1]);
	
	// global width of the problem
	dsize[1] = atoi(argv[3]);
	if ( dsize[1]%psize[1] != 0) {
		fprintf(stderr, "Error: invalid problem width\n");
		abort();
	}
	// width of the local data block (add boundary or ghost zone: 1 point on each side)
	dsize[1]  = dsize[1]/psize[1]  + 2;
	
	// global height of the problem
	dsize[0] = atoi(argv[2]);
	if ( dsize[0]%psize[0] != 0) {
		fprintf(stderr, "Error: invalid problem height\n");
		abort();
	}
	// height of the local data block (add boundary or ghost zone: 1 point on each side)
	dsize[0] = dsize[0]/psize[0] + 2;
	
	// creation of the communicator
	int cart_period[2] = { 0, 0 };
	MPI_Cart_create(MPI_COMM_WORLD, 2, psize, cart_period, 1, cart_comm);
}

int main( int argc, char* argv[] )
{
	// initialize the MPI library
	MPI_Init(&argc, &argv);
	
	// parse the command line arguments
	int nb_iter;
	int dsize[2];
	MPI_Comm cart_comm;
	parse_args(argc, argv, &nb_iter, dsize, &cart_comm);
	
	// find the coordinate of the local process
	int pcoord_1d; MPI_Comm_rank(MPI_COMM_WORLD, &pcoord_1d);
	int pcoord[2]; MPI_Cart_coords(cart_comm, pcoord_1d, 2, pcoord);
	
	// allocate data for the current iteration
	double(*cur)[dsize[1]]  = malloc(sizeof(double)*dsize[1]*dsize[0]);
	
	// initialize data at t=0
	init(dsize, pcoord, cur);
	
	// allocate data for the next iteration
	double(*next)[dsize[1]] = malloc(sizeof(double)*dsize[1]*dsize[0]);

	//create filename to output the data
	char filename[42];
	sprintf(filename, "heat%dx%d.h5", pcoord[0], pcoord[1]);
	//printf("%s\n", filename);
	
	//create the file and in each of them a data space
	hid_t file_id = create_file(filename);
	//hsize_t dimspace[] = {dsize[0], dsize[1]};
	hid_t mem_dataspace = create_space(dsize[0], dsize[1]);
	hid_t dataspace_f = create_space(dsize[0]-2, dsize[1]-2);


	//EXO 2 : using hyperslab to remove ghosts
	hsize_t start[] = {1, 1};
	hsize_t count[] = {dsize[0]-2, dsize[1]-2};
	hid_t status_hyperslab = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, start, NULL, count, NULL);

	char step[42];
	hid_t dataset;
	herr_t status;

	// the main (time) iteration
	for (int ii=0; ii<nb_iter; ++ii) {
		//create the step name and dataset to save a dataset at each iteration
		sprintf(step, "/step%d", ii);
		dataset = create_dataset(file_id, step, dataspace_f);
		

		//EXO 1 with ghosts
		//status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cur);
		
		//EXO 2 without ghosts
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, mem_dataspace, dataspace_f, H5P_DEFAULT, cur);
		

		//close the dataset
    	H5Dclose(dataset);
		
		// compute the temperature at the next iteration
		iter(dsize, cur, next);
		
		// update ghost zones
		exchange(cart_comm, dsize, next);
		
		// switch the current and next buffers
		double (*tmp)[dsize[1]] = cur; cur = next; next = tmp;
	}

	//save the last value, at step = max+1
	sprintf(step, "/step%d", nb_iter);
	dataset = create_dataset(file_id, step, dataspace_f);
	

	//EXO 1 with ghosts
	//status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cur);
	
	//EXO 2 without ghosts
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, mem_dataspace, dataspace_f, H5P_DEFAULT, cur);
	

	//close the dataset
	H5Dclose(dataset);

	//close dataspace and file
	H5Sclose(dataspace_f);
	H5Sclose(mem_dataspace);
	H5Fclose(file_id);
	
	// free memory
	free(cur);
	free(next);

	// finalize MPI
	MPI_Finalize();
	
	return 0;
}