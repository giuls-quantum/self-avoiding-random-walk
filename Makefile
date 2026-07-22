CC = gcc
CFLAGS = -Wall -Wextra -O3 -march=native
TARGET = bin/sarw

all: $(TARGET)

$(TARGET): src/selfa_random_walk.c
	mkdir -p bin
	$(CC) $(CFLAGS) -o $(TARGET) src/selfa_random_walk.c

clean:
	rm -rf bin/
