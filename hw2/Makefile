CC = mpic++
CCFLAGS = -Iinclude --std=c++11
BLUEGENE_CC = mpixlcxx_r
BLUEGENE_FLAGS = -Iinclude 

all: function src/main.cpp
	$(CC) src/main.cpp build/function.o -o bin/main $(CCFLAGS)

bluegene: src/main.cpp folders
	$(BLUEGENE_CC) src/main.cpp src/function.cpp $(BLUEGENE_FLAGS) -o bin/bluegene

function: folders src/function.cpp include/function.h
	$(CC) -c src/function.cpp -o build/function.o $(CCFLAGS)

folders:
	mkdir -p bin build

clean:
	rm -rf bin build
