
example1: example1.f08
	gfortran -o $@ $<

example2: example2.o mystuff.o
	gfortran -o $@ $^

example2.o: example2.f08 mystuff.mod
	gfortran -c $<

mystuff.o mystuff.mod: mystuff.f08
	gfortran -c $< 

example4: example4.o quadform.o
	gcc -o $@ $^ -lm

example4.o: example4.c quadform.h
	gcc -c $<

quadform.o: quadform.c
	gcc -c $<

clean:
	rm -f example1 example2 example4
	rm -f *.o *.mod a.out *.s
