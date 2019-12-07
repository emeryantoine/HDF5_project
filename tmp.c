MPI_Init();

	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
		hid_t file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

		
		hid_t filespace = H5Screate_simple(2, dimsf, NULL);

			hid_t dataspace_f = create_space(50, 50);

			//EXO 2 : using hyperslab to remove ghosts
			hsize_t start[] = {1, 1};
			hsize_t count[] = {50, 50};
			hid_t status_hyperslab = H5Sselect_hyperslab(mem_dataspace, H5S_SELECT_SET, start, NULL, count, NULL);

			hid_t dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
				hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
				H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
				H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, data);


			H5Dclose(dset_id);

		H5Sclose(filespace);

	H5Pclose(plist_id);
	H5Fclose(file_id);

MPI_Finalize();