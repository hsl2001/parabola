CC = gcc
CFLAGS = -O3 -march=native -Wall -Wextra
LDFLAGS = -lm -lpthread -lz
TARGET = reverb
SRCS = reverb.c klib/kthread.c

.PHONY: clean static

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

all: $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

static: $(SRCS)
	$(CC) $(CFLAGS) -static -o $(TARGET) $(SRCS) $(LDFLAGS)

clean:
	rm -f $(TARGET)