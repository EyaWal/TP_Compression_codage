INCLUDE = ./

COMPILATEUR = g++

OPTIONS = -g -lm

SOURCES = ex2.cpp

OBJECTS = ex2.o matrix.o fichiers.o pred.o dct.o

EXECUTABLE = ex2

$(EXECUTABLE): $(OBJECTS)
	$(COMPILATEUR) $(OPTIONS) $(OBJECTS)  -o $(EXECUTABLE)

ex2.o: ex2.cpp
	$(COMPILATEUR) $(OPTIONS) -c ex2.cpp

