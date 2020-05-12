
all: o lib

o:
	g++ -c  -fpic liblombscargle.cc ComplexLS.c
	
lib:
	g++ -shared -o liblombscargle.so liblombscargle.o ComplexLS.o -lm
	
clean:
	rm -vf *.o *.so
	
