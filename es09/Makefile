CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: main.exe
main.exe : GA.o main.o trajectory.o random.o
	$(CC) GA.o trajectory.o random.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp GA.h trajectory.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
GA.o : GA.cpp GA.h random.h trajectory.h
	$(CC) -c GA.cpp -o GA.o $(CFLAGS)
trajectory.o : trajectory.cpp trajectory.h random.h
	$(CC) -c trajectory.cpp -o trajectory.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe
