CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic 

SRC      := $(wildcard *.cpp)
BINS     := $(SRC:.cpp=)

.PHONY: all clean

all: $(BINS)

%: %.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $< -o bin/$@

clean:
	rm -f $(BINS)

