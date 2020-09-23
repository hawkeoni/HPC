CC = /usr/local/Cellar/gcc/9.3.0/bin/g++-9
CCFLAGS = -Iinclude

all: parser folders src/main.cpp
	$(CC) src/main.cpp build/parser.o -o bin/main $(CCFLAGS)


parser: folders src/parser.cpp include/parser.h
	$(CC) -c src/parser.cpp -o build/parser.o $(CCFLAGS)


folders:
	mkdir -p bin build

clean:
	rm -rf bin build
