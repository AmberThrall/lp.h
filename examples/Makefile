CXX      := g++
CXXFLAGS := -O2 -Wall -Wextra -pedantic -fsanitize=address -g

SRC      := $(wildcard *.cpp)
BINS     := $(SRC:.cpp=)

.PHONY: all clean

all: $(BINS)

%: %.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $< -o bin/$@

clean:
	rm -f $(BINS)

