EXEC=heat
CC=h5pcc
CFLAGS=-W -Wall -Werror -ansi -pedantic

$(EXEC):$(EXEC).c
	$(CC) -o $@ $^

test:
	mpirun -n 4 ./$(EXEC) 4 100 100
	echo "dump 0 0"
	h5dump heat0x0.h5
	echo "dump 0 1"
	h5dump heat0x1.h5
	echo "dump 1 0"
	h5dump heat1x0.h5
	echo "dump 1 1"
	h5dump heat1x1.h5

run:
	mpirun -n 4 ./$(EXEC) 100 100 100

show:
	echo "dump 0 0"
	h5dump heat0x0.h5
	echo "dump 0 1"
	h5dump heat0x1.h5
	echo "dump 1 0"
	h5dump heat1x0.h5
	echo "dump 1 1"
	h5dump heat1x1.h5



clean:
	rm -rf *.o

clean+:
	rm -rf *.o
	rm $(EXEC)
	rm -rf *.h5