# CFLAGS for running
CFLAGS = -Wall -v -mtune=native -msse4.2 -O2

# CFLAGS for debugging
#CFLAGS = -Wall -v -O0 -g

CC = g++
RM = rm -rf
LIBS =-lboost_thread

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)


.cpp.o:
	 $(CC) -c $(CFLAGS) $*.cpp


all : mc-boost 


mc-boost: $(OBJS)
	 $(CC) -o  $@ $(OBJS) $(CFLAGS) $(LIBS)


clean::
	 $(RM) mc-boost
	 $(RM) *.o

