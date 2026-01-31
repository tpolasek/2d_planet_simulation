CXX = /opt/homebrew/opt/llvm/bin/clang++
CXXFLAGS = -std=c++23 -Wall -Wextra -fopenmp -O3 -ffast-math
INCLUDES = -I/opt/homebrew/include -I/opt/homebrew/include/SDL2 -I/opt/homebrew/opt/libomp/include
LDFLAGS = -L/opt/homebrew/lib -L/opt/homebrew/opt/libomp/lib
LIBS = -lSDL2 -lSDL2_ttf

TARGET = pixel_gravity
SOURCES = main.cpp simulation.cpp utility.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)
