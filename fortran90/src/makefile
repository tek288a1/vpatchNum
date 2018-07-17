execute:
	@gfortran -O3 -o main.exe mymod.f90 main.f90 routines.f90 -framework Accelerate
	@./main.exe | tee output.log

clean:
	@rm -f *.o
