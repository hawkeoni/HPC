CC = mpic++
CCFLAGS = -Iinclude 

all: function src/main.cpp
	$(CC) src/main.cpp build/function.o -o bin/main $(CCFLAGS)

function: folders src/function.cpp include/function.h
	$(CC) -c src/function.cpp -o build/function.o $(CCFLAGS)

folders:
	mkdir -p bin build

clean:
	rm -rf bin build
