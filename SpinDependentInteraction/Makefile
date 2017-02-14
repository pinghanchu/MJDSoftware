# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall

main : spindependentinteraction.o
	$(CC) spindependentinteraction.o  -o spindependentinteraction
main.o: spindependentinteraction.cc
	$(CC) $(CFLAGS) spindependentinteraction.cc
clean:
	rm -rf *o

