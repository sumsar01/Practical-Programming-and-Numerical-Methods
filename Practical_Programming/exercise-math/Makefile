CFLAGS = -std=gnu99 -Wall -O

.PHONEY : default
default : out.txt
	cat out.txt

out.txt :a.out
	./a.out > out.txt

a.out : math.c
	cc $(CFLAGS) math.c -o a.out

.PHONEY : clean
clean:
	$(RM) a.out out.txt  
