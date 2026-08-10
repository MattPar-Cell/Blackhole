# Plain make build, for when cmake is not around.
#   make          build
#   make test     build and run the physics validation suite
#   make images   render the sample images used in the README

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -Isrc
LDFLAGS  ?= -pthread

SRC  := src/main.cpp src/render.cpp src/spectrum.cpp src/image.cpp src/validate.cpp
OBJ  := $(SRC:.cpp=.o)
BIN  := blackhole

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

src/main.o:     src/render.h src/disc.h src/scene.h src/image.h src/kerr.h src/geodesic.h src/constants.h
src/render.o:   src/render.h src/disc.h src/scene.h src/image.h src/kerr.h src/geodesic.h src/spectrum.h
src/spectrum.o: src/spectrum.h src/constants.h
src/image.o:    src/image.h src/spectrum.h
src/validate.o: src/kerr.h src/geodesic.h src/disc.h src/spectrum.h src/constants.h

test: $(BIN)
	./$(BIN) --test

images: $(BIN)
	./$(BIN) --preset sgra --sun --width 1600 --height 900 --spp 3 --out gallery/sgra_sun.png
	./$(BIN) --preset stellar --width 1600 --height 900 --spp 3 --out gallery/stellar.png
	./$(BIN) --preset gargantua --width 1600 --height 900 --spp 3 --out gallery/gargantua.png

clean:
	rm -f $(OBJ) $(BIN)

.PHONY: all test images clean
