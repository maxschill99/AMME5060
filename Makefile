MPI = mpif90

Prog: outputmod_MP.o varmod_MP.o nodemod_MP.o partitioning_MP.o residuals.o init.o jacobi.o main.o
	$(MPI) outputmod_MP.f90 varmod_MP.f90 nodemod_MP.f90 partitioning_MP.f90 residuals.f90 init.f90 jacobi.f90 main.f90 -o test1.exe

main.o: outputmod_MP.o varmod_MP.o nodemod_MP.o partitioning_MP.o residuals.o init.o jacobi.o main.f90
	$(MPI) -c main.f90
jacobi.o: outputmod_MP.o varmod_MP.o nodemod_MP.o partitioning_MP.o residuals.o init.o jacobi.f90
	$(MPI) -c jacobi.f90
init.o: outputmod_MP.o varmod_MP.o nodemod_MP.o partitioning_MP.o residuals.o init.f90
	$(MPI) -c init.f90
residuals.o: outputmod_MP.o varmod_MP.o nodemod_MP.o partitioning_MP.o residuals.f90
	$(MPI) -c residuals.f90
partitioning_MP.o: outputmod_MP.o varmod_MP.o nodemod_MP.o partitioning_MP.f90
	$(MPI) -c partitioning_MP.f90
nodemod_MP.o: outputmod_MP.o varmod_MP.o nodemod_MP.f90
	$(MPI) -c nodemod_MP.f90
varmod_MP.o: outputmod_MP.o varmod_MP.f90
	$(MPI) -c varmod_MP.f90
outputmod_MP.o: outputmod_MP.f90
	$(MPI) -c outputmod_MP.f90

