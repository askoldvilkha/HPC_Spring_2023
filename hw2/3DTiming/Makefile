OBJS = driver.o looptest.o
PROG = driver
## You have to change the name of the PROGS file every time you try to create a different driver. 
## Make sure to keep the program names consistent with a driver type
## For example, for driver8zzrev you just change the PROG entry to driver8zzrev (PROG = driver8zzrev) (not driver8zzrev.c!)

## The following cryptic command means:
## If you have a file called BLAH.c and need BLAH.o then
## run the command gcc -Wall -O3 -ffast-math std=c99 -c BLAH.c
## make understands based on this syntax that BLAH.o depends 
## in BLAH.c . If BLAH.c changes, then it will be recompiled to
## BLAH.o, otherwise it won't
%.o: %.c
	gcc -Wall -O3 --fast-math -std=c99 -c $<


## This means that the file stored in the variable PROG (here, driver) depends
## on all the files stored in the variable OBJS (here, driver.o and looptest.o)
## And, to generate driver, make needs to run the command
##  gcc -Wall driver.o looptest.o -o driver

$(PROG): $(OBJS)
	gcc -Wall $(OBJS) -o $@


## .PHONY is a hint to make taht the commands don't produce files with that name
## i.e., "make clean" should do some cleanup, not produce a file called "clean"
.PHONY: clean
.PHONY: all

all: $(PROG)

## running "make clean" will cause the command "rm driver.o looptest.o driver"
## to be run. The extra "-" sign in front tells make to ignore errors
clean:
	-rm $(OBJS) $(PROG)
