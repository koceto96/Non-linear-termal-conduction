F95 = mpif90
OPTS = -O1 -Wall
LIBS=-L${BLASDIR} ${BLASLIB} -lpthread

OBJS = header.o q8.o laplace.o maximum.o newton.o func.o jacobian.o cg.o matmult.o sparsegather.o vecdot.o

q8: $(OBJS)
	$(F95) $(OPTS) $(OBJS) -o q8 $(LIBS)

%.o: %.f90
	$(F95) $(OPTS) -c $<

clean:
	rm -rf q8 *.o *.mod core solution.* slurm-*

