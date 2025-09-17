CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors

symnmf: symnmf.o functions.o
	$(CC)  -o symnmf symnmf.o functions.o $(CFLAGS) -lm

symnmf.o: symnmf.c functions.c	
	$(CC) -c symnmf.c functions.c $(CFLAGS) -lm