CC = gcc
CFLAGS = -Wall -g -std=c99
LIBS = -lm 
INCLS = -I./ 


read_envi:	read_envi.o strparse.o
		${CC} ${CFLAGS} $@.o -o $@ strparse.o ${INCLS} ${LIBS}


.c.o: $<
		$(CC) ${INCLS} $(CFLAGS) -c $<

clean:
		\rm -f *.o *~ *%
