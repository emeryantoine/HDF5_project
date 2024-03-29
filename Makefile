EXEC=heat
EXEC2=heat2
CC=h5pcc
CFLAGS=-W -Wall -Werror -ansi -pedantic
width=100
height=400
step=10


$(EXEC2):$(EXEC2).c
	$(CC) -o $@ $^


$(EXEC):$(EXEC).c
	$(CC) -o $@ $^

run:
	mpirun -n 4 ./$(EXEC) $(step) $(width) $(height)

run2:
	mpirun -n 4 ./$(EXEC2) $(step) $(width) $(height)

show:
	echo "dump 0 0"
	h5dump heat0x0.h5
	echo "dump 0 1"
	h5dump heat0x1.h5
	echo "dump 1 0"
	h5dump heat1x0.h5
	echo "dump 1 1"
	h5dump heat1x1.h5

show2:
	echo "dump parallel write"
	h5dump heat.h5



clean:
	rm -rf *.o

clean+:
	rm -rf *.o
	rm -rf $(EXEC)
	rm -rf $(EXEC2)
	rm -rf *.h5