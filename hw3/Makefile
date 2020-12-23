CC = mpic++
NVCC = nvcc
NVCC_FLAGS = -arch=sm_60
POLUS_CC = mpixlC 
CCFLAGS = -Iinclude -Isrc -O2

all: function src/main.cu
	$(NVCC) $(NVCC_FLAGS) -ccbin $(POLUS_CC) src/main.cu build/function.o -o bin/polus $(CCFLAGS)

function: folders src/function.cpp include/function.h
	$(CC) -c src/function.cpp -o build/function.o $(CCFLAGS)

folders:
	mkdir -p bin build

clean:
	rm -rf bin build
