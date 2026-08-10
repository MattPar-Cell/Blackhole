# Plain make build, for when cmake is not around.
#   make          build
#   make test     build and run the physics validation suite
#   make images   render the sample stills in gallery/
#   make videos   render the sample films in video/ (slow: ~30 min)

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
	@mkdir -p gallery
	./$(BIN) --preset sgra      --width 1200 --height 675 --spp 3 --out gallery/sgra.png
	./$(BIN) --preset m87       --width 1000 --height 750 --spp 3 --out gallery/m87.png
	./$(BIN) --preset stellar   --width 1200 --height 675 --spp 3 --out gallery/stellar.png
	./$(BIN) --preset gargantua --width 1200 --height 675 --spp 3 --out gallery/gargantua.png
	./$(BIN) --preset sgra --no-disc --fov 42 --stars 400000 \
		--width 1200 --height 675 --spp 3 --out gallery/lensed-starfield.png
	./$(BIN) --preset sgra --sun --no-disc --stars 400000 \
		--tonemap log --log-decades 13 --glare 0 \
		--width 1600 --height 700 --spp 3 --out gallery/sun-comparison.png
	./$(BIN) --preset sgra --sun --width 1400 --height 620 --spp 2 \
		--out gallery/sun-with-disc.png

# Each film is 240 renders, so this takes roughly half an hour on four cores.
videos: $(BIN)
	@mkdir -p video
	./$(BIN) --preset sgra --inclination 60 --fov 34 \
		--hotspot 9 --hotspot-size 1.1 --hotspot-contrast 200 \
		--orbits 2 --duration 10 --fps 24 \
		--width 400 --height 300 --spp 1 --no-stars \
		--out video/hotspot-orbit.apng
	./$(BIN) --preset sgra --inclination 4 --inclination-to 89 --fov 34 \
		--duration 10 --fps 15 --width 400 --height 250 --spp 1 --no-stars \
		--out video/inclination-sweep.apng

clean:
	rm -f $(OBJ) $(BIN)

.PHONY: all test images videos clean
