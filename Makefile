CC = gcc
LOCAL = local
GLOBAL = global
FITTING = fitting
GLOBALAFFINE = globalAffine
ORTHOLOGYFINDER = orthologyFinder
CFLAGS = -g -Wall -Wextra -std=c99 -lm

LOCAL_OBJS = local.c alignment.o
GLOBAL_OBJS = global.c alignment.o
FITTING_OBJS = fitting.c fittingFunctions.o alignment.o
GLOBALAFFINE_OBJS = globalAffine.c alignment.o
ORTHOLOGYFINDER_OBJS = fittingFunctions.o orthologyFinder.o alignment.o

all: local global fitting globalAffine orthologyFinder

local: $(LOCAL_OBJS)
	$(CC) $(CFLAGS) -o $(LOCAL) $(LOCAL_OBJS)

global: $(GLOBAL_OBJS)
	$(CC) $(CFLAGS) -o $(GLOBAL) $(GLOBAL_OBJS)

fitting: $(FITTING_OBJS)
	$(CC) $(CFLAGS) -o $(FITTING) $(FITTING_OBJS)

globalAffine: $(GLOBALAFFINE_OBJS)
	$(CC) $(CFLAGS) -o $(GLOBALAFFINE) $(GLOBALAFFINE_OBJS)

orthologyFinder: $(ORTHOLOGYFINDER_OBJS)
	$(CC) $(CFLAGS) -o $(ORTHOLOGYFINDER) $(ORTHOLOGYFINDER_OBJS)

alignment.o: alignment.c alignment.h
	$(CC) $(CFLAGS) -c alignment.c -o alignment.o

fittingFunctions.o: fittingFunctions.c fitting.h
	$(CC) $(CFLAGS) -c fittingFunctions.c -o fittingFunctions.o

clean:
	rm -f $(LOCAL) $(GLOBAL) $(FITTING) $(GLOBALAFFINE) $(ORTHOLOGYFINDER)
