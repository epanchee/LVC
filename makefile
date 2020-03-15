all: LR.o AP.o main.o tipHans.o TP06-tab.o
	icc -std=c99 -o LVD LR.o AP.o main.o tipHans.o TP06-tab.o -fopenmp -lm -O2 -g
	rm *.o

LR.o:
	icc -std=c99 -c LR.c -O2 -g

AP.o:
	icc -std=c99 -c AP.c -O2 -g

main.o: 
	icc -std=c99 -c main.c -fopenmp -O2 -g

tipHans.o:
	icc -std=c99 -c tipHans.c -fopenmp -O2 -g

TP06-tab.o:
	icc -std=c99 -c TP06-tab.c -O2 -g

clean: 
	rm *.o LVD