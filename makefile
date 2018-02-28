CC = g++

SRCS = $(wildcard *.cpp)

CFLAGS = -Wall -Wextra -Werror -ggdb3 -O0 -std=c++11

PROGS = $(patsubst %.cpp,%,$(SRCS))

all: $(PROGS)

%: %.cpp
	$(CC) $(CFLAGS) -o $@ $< 
clean: 
	rm -f $(PROGS)
