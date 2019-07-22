
all: o lib

o:
	g++ -c  -fpic liblombscargle.cc
	
lib:
	g++ -shared -o liblombscargle.so liblombscargle.o -lm
	
clean:
	rm -vf *.o *.so
	
