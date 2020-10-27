all: obj lib

obj:
	g++ -c -fPIC libfilter.cc -o libfilter.o


lib: 
	g++ -fPIC -shared -o libfilter.so libfilter.o -lm
	rm -v *.o

clean:
	rm -v  *.so

