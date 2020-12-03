CC = mpic++
CCFLAGS = -Iinclude -fopenmp
BLUEGENE_CC = mpixlcxx_r
POLUS_CC = mpixlC 

all: function src/main.cpp
	$(CC) src/main.cpp build/function.o -o bin/main $(CCFLAGS)

bluegene: src/main.cpp folders
	$(BLUEGENE_CC) src/main.cpp src/function.cpp $(CCFLAGS) -o bin/bluegene

polus: src/main.cpp folders
	$(POLUS_CC) src/main.cpp src/function.cpp $(CCFLAGS) -o bin/polus

function: folders src/function.cpp include/function.h
	$(CC) -c src/function.cpp -o build/function.o $(CCFLAGS)

folders:
	mkdir -p bin build

clean:
	rm -rf bin build
