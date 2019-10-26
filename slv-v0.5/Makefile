
# CC = /opt/local/bin/gcc-mp-4.6
CC = gcc -std=c99
CC = gcc -std=c99 -O3 -Wgnu -Wall -pedantic -funroll-loops -fomit-frame-pointer -fexpensive-optimizations  
CC = gcc -std=c99 -O3 -DNDEBUG  -Wall -Winline  -funroll-loops    

RM = rm


GMP_INC_DIR=.
GMP_LIB_DIR=.

LIB_PATHS = -L/opt/local/lib/ -L/usr/local/lib -L/Users/elias/local/lib -L$(GMP_LIB_DIR)
INC_PATHS = -I/opt/local/include -I/Users/elias/local/include -I/usr/local/include -I$(GMP_INC_DIR)

#CFLAGS = -O3 -Wall -ansi -pedantic -funroll-loops -fomit-frame-pointer -fexpensive-optimizations  
CFLAGS = -c -g -Wall -Wshadow -Wpointer-arith -Wwrite-strings
#-Wpedantic 

CFLAGS += $(INC_PATHS)

DEBUG = -g
LDFLAGS = -Wall $(DEBUG) $(LIB_PATHS) -lm -lgmp  -lmpfr -lmpfi

LIB_H=info.h interval.h poly_ops.h utils.h sys_queue.h vca_solver_1.h 	
LIB_SRC=info.c interval.c poly_ops.c vca_solver_1.c 	

LIB_OBJS=$(patsubst %.c,%.o,$(LIB_SRC))


all: slv



slv: slv.o $(LIB_OBJS)
	ls $(LIB_OBJS)
	$(CC)  $@.o $(LIB_OBJS) -o $@ $(LDFLAGS)

slv.o: slv.c $(LIB_H) $(LIB_SRC)
	$(CC) -c slv.c $(CFLAGS)


vca_solver_2.o: vca_solver_2.h  vca_solver_2.c vca_solver_1.h vca_solver_1.c  interval.h interval.c dbg.h
	$(CC) $(CFLAGS) vca_solver_2.c

vca_solver_1.o: vca_solver_1.h vca_solver_1.c  interval.h interval.c dbg.h
	$(CC) $(CFLAGS) vca_solver_1.c

poly_ops.o: poly_ops.h poly_ops.c interval.h interval.c dbg.h
	$(CC) $(CFLAGS) poly_ops.c


interval.o : interval.h interval.c dbg.h
	$(CC) $(CFLAGS) interval.c


.PHONY: clean

clean:
	$(RM) -rf  slv  *.o *~ 
