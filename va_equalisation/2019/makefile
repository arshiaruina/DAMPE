CC = g++
TARGET = va_equalisation
IDIR = ./inc
ODIR = ./obj
SDIR = ./src

all: $(TARGET)

DMPIDIR = $(DMPSWSYS)/include
DMPLIBS = $(DMPSWSYS)/lib
RCFLAGS = `root-config --cflags`
RLIBS = `root-config --libs`

CFLAGS = -Wall $(RCFLAGS) -I$(DMPIDIR) -I$(IDIR)
LIBS = $(RLIBS) -L$(DMPLIBS) -lDmpEvent -lDmpService

_SRCS = va_equalisation.cpp
SRCS = $(patsubst %,$(SDIR)/%,$(_SRCS))

_DEPS = va_equalisation.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJS = va_equalisation.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LIBS) $^ -o $@

$(OBJS): $(SRCS) $(DEPS)
	$(CC) -c $(CFLAGS) $(LFLAGS) $< -o $@

clean:
	rm -f eta ./*~ $(SDIR)/*~ $(IDIR)/*~ $(OBJS)
