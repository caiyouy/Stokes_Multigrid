sparseMatrix.o:sparseMatrix/sparseMatrix.c
	gcc -g $< -o $@ -c

main.o: main.c
	gcc -g $< -o $@ -c

main:main.o sparseMatrix.o
	gcc -g $^ -o $@ -lm