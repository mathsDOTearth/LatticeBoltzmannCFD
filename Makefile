GCC=gcc
CFLAGS=-O3 -o 
LDFLAGS=-lm -lSDL2 -fopenmp
TARGET=lbm
SRC=main.c lbm.c 

build: $(SRC) lbm.h
	$(GCC) $(CFLAGS) $(TARGET) $(SRC) $(LDFLAGS)

run: build
	./$(TARGET)
clean:
	rm -f $(TARGET)
