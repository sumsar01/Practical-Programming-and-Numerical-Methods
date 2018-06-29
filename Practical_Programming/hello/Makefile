CFLAGS = -std=gnu99 -Wall -O

.PHONEY : all
all : out.txt
	cat out.txt

out.txt : hello
	./$< > out.txt

hello: hello.o

hello.o:

clean:
	$(RM) hello out.txt
