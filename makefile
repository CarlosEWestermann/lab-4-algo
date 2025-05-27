CC = g++
CFLAGS = -O2
SRC =  hungarian.cpp
TARGET = hungarian

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)