CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
CC = gcc

.PHONEY:all
all : out.txt
	cat out.txt

out.txt : main
	./main > out.txt


main : main.o nvector.o
	$(CC) $(CFLAGS) main.o nvector.o -o main $(LDLIBS)

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

nvector.o : nvector.c
	$(CC) $(CFLAGS) -c nvector.c -o nvector.o


clean :
	$(RM) out.txt main main.o nvector.o

test :
	echo $(CFLAGS)
