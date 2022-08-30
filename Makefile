# Makefile for the extension.
#
# Andres Chamorro 
# This code is in the public domain.
CC = gcc
CFLAGS = -g -pthread -fwrapv -Wall -Wno-unused-result -Wstrict-prototypes -fPIC
LDFLAGS = -lz -lm -pthread -shared

# Point PYTHON_DIR to a place where Python headers
PYTHON_DIR =$(shell python -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())" ) 
INCLUDE = -I$(PYTHON_DIR)/include -I$(PYTHON_DIR)

all: ngsim.so

readmodule.o: readmodule.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $^ -o $@

ngsim.so: readmodule.o
	$(CC) $(LDFLAGS) $^ -o $@

clean:
	rm -rf *.o *.so build dist
