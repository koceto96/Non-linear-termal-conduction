F95 = mpif90
OPTS = -O3 -Wall
LIBS=-L${BLASDIR} ${BLASLIB} -lpthread

OBJS = header.o main.o laplace.o save_solution.o maximum.o newton.o func.o jacobian.o cg.o matmult.o sparsegather.o vecdot.o

main: $(OBJS)
	$(F95) $(OPTS) $(OBJS) -o main $(LIBS)

%.o: %.f90
	$(F95) $(OPTS) -c $<

clean:
	rm -rf main *.o *.mod core solution.* slurm-*

