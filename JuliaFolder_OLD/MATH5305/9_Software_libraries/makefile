BLAS = blas

example1: example_my_fblas.o my_fblas.o
	gfortran -o $@ $^

example_my_fblas.o: example_my_fblas.f08 my_fblas.mod
	gfortran -c $<

my_fblas.o my_fblas.mod: my_fblas.f08
	gfortran -c -fpic $<

example2: example_my_cblas.o my_cblas.o
	gcc -o $@ $^

example_my_cblas.o: example_my_cblas.c my_cblas.h
	gcc -c $<

my_cblas.o: my_cblas.c
	gcc -c -fpic $<

example3: example_my_fblas.o libmyblas.so
	gfortran -o $@ $^ -Wl,-rpath=.

libmyblas.so: my_fblas.o my_cblas.o
	gfortran -o $@ -shared $^

example4: example_my_cblas.o libmyblas.so
	gcc -o $@ $< -L. -lmyblas \
	-Wl,-rpath=.

example5: example_system_fblas.o 
	gfortran -o $@ $< -l$(BLAS)

example_system_fblas.o: example_system_fblas.f08
	gfortran -c $<

example6: example_system_cblas.o
	gcc -o $@ $< -l$(BLAS)

example_system_cblas.o: example_system_cblas.c
	gcc -c $<

clean:
	rm -f *.o *.so *.mod
	rm -f example1 example2 example3 example4 example5 example6
