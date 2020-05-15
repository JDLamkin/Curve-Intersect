SLV_PATH = ../slv-v0.5

CPP = g++ -std=c++17 -O3 -DNDEBUG  -Wall -Winline  -funroll-loops

RM = rm -rf

GMP_INC_DIR=.
GMP_LIB_DIR=.

LIB_PATHS = -L/opt/local/lib/ -L/usr/local/lib -L$(GMP_LIB_DIR)
INC_PATHS = -I/opt/local/include -I/usr/local/include -I$(GMP_INC_DIR) -I$(SLV_PATH)

#CFLAGS = -O3 -Wall -ansi -pedantic -funroll-loops -fomit-frame-pointer -fexpensive-optimizations  
CFLAGS = -c -g -Wall -Wshadow -Wpointer-arith -Wwrite-strings
#-Wpedantic 

CFLAGS += $(INC_PATHS)

DEBUG = -g
LDFLAGS = -Wall $(DEBUG) $(LIB_PATHS) -lm -lgmp -lmpfr -lgmpxx #-lmpfi

SLV_LIB_H=info.h interval.h poly_ops.h utils.h sys_queue.h vca_solver_1.h 	
SLV_LIB_SRC=info.c interval.c poly_ops.c vca_solver_1.c 	
CCI_LIB_H=slvpp.h isolPnt.h polynom.h poly2.h utils.h rational.h
CCI_LIB_SRC=slvpp.cpp isolPnt.cpp polynom.cpp poly2.cpp

LIB_H=$(patsubst %,$(SLV_PATH)/%,$(SLV_LIB_H)) $(CCI_LIB_H)
LIB_OBJS=$(patsubst %.c,$(SLV_PATH)/%.o,$(SLV_LIB_SRC)) $(patsubst %.cpp,%.o,$(CCI_LIB_SRC))

.PHONY: clean all

all: cci

clean:
	$(RM) cci *.o *~

$(SLV_PATH)/%.o:
	cd $(SLV_PATH) && $(MAKE) %.o

%.o: %.cpp $(LIB_H)
	$(CPP) $(CFLAGS) $< -o $@

cci: cppTest.o $(LIB_OBJS)
	$(CPP) $^ -o $@ $(LDFLAGS)
