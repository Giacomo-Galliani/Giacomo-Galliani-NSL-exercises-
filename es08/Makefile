CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o SA.o
	$(CC) random.o main.o SA.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
SA.o : SA.cpp SA.h 
	$(CC) -c SA.cpp -o SA.o $(CFLAGS) 
clean :
	rm *.o main.exe seed.out
