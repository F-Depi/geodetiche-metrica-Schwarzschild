CC = gcc
CFLAGS = -I${CURDIR}/include

OBJECTS_TEST = src/test.c src/TESI_fun.c
OBJECTS_MAIN = src/main.c src/TESI_fun.c


test: $(OBJECTS_TEST)
	$(CC) $(CFLAGS) -o test.x $(OBJECTS_TEST) -lm 

main: $(OBJECTS_MAIN)
	$(CC) $(CFLAGS) -o main.x $(OBJECTS_MAIN) -lm 
