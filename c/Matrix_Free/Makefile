main:main.o sparseMatrix.o
	gcc -g $^ -o $@ -lm

run:
	./main

all:main run

sparseMatrix.o:sparseMatrix/sparseMatrix.c
	gcc -g $< -o $@ -c

main.o: main.c
	gcc -g $< -o $@ -c

.PHONY:clean realclean
clean:
	-rm *.o *.csv

realclean:
	-rm *.o *.csv main