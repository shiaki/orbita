
CC      = gcc
CFLAGS  = -c -Wall -O3
LDFLAGS = -lm -lgsl -lgslcblas

DEPS    = analyzer.h chaos.h disk.h fcsintp.h ferrers.h \
          integrator.h models.h orbit.h periodic_seeker.h poincare.h \
          potential.h runge_kutta.h sphintp.h utils.h
SOURCES = analyzer.c chaos.c disk.c fcsintp.c ferrers.c \
          integrator.c models.c orbit.c periodic_seeker.c poincare.c \
          potential.c runge_kutta.c sphintp.c utils.c

OBJECTS = $(SOURCES:.c=.o)
LIBRARY = orbita

all: $(SOURCES) $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	ar cr liborbita.a $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o
