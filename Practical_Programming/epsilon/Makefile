CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
CC = gcc

.PHONEY : all 
all : out.txt
	cat out.txt

out.txt : main
	./main > out.txt

main : main.o equal.o
	$(CC) $(CFLAGS) main.o equal.o -o main $(LDLIBS)

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

equal.o : equal.c
	$(CC) $(CFLAGS) -c equal.c -o equal.o

clean:
	$(RM) out.txt main main.o equal.o

test : echo $(CFLAGS)
