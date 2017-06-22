# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall

main : main.o
	$(CC) main.o  -o main
main.o: main.cc
	$(CC) $(CFLAGS) main.cc
clean:
	rm -rf *o

