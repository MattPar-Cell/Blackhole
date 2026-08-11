# Plain make build, for when cmake is not around.
#   make          build
#   make test     build and run the physics validation suite
#   make images   render the 20 stills in gallery/ (slow: ~3 h, 8 are 4K)
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
# The flagship frames are rendered at 4K.  Note that 3840x2160 at --spp 2
# samples *finer* than 1920x1080 at --spp 3 - the sample spacing is fov/7680
# against fov/5760 - so this is higher quality as well as bigger.
	./$(BIN) --preset ton618    --width 3840 --height 2160 --spp 2 --out gallery/ton618.png
	./$(BIN) --preset sgra      --width 3840 --height 2160 --spp 2 --out gallery/sgra.png
	./$(BIN) --preset m87       --width 2880 --height 2160 --spp 2 --out gallery/m87.png
	./$(BIN) --preset stellar   --width 3840 --height 2160 --spp 2 --out gallery/stellar.png
	./$(BIN) --preset gargantua --width 3840 --height 2160 --spp 2 --out gallery/gargantua.png
	./$(BIN) --preset ton618 --inclination 80 --fov 30 \
		--width 3840 --height 2160 --spp 2 --out gallery/ton618-edge.png
	./$(BIN) --preset sgra --inclination 5  --fov 26 \
		--width 1440 --height 1080 --spp 2 --out gallery/face-on.png
	./$(BIN) --preset sgra --inclination 60 --fov 30 \
		--width 3840 --height 2160 --spp 2 --out gallery/inclined.png
	./$(BIN) --preset sgra --inclination 90 --fov 30 \
		--width 1920 --height 1080 --spp 3 --out gallery/edge-on.png
	./$(BIN) --preset sgra --spin 0.0   --inclination 80 --fov 26 --disc-outer 22 \
		--width 1280 --height 720 --spp 2 --out gallery/spin-00.png
	./$(BIN) --preset sgra --spin 0.5   --inclination 80 --fov 26 --disc-outer 22 \
		--width 1280 --height 720 --spp 2 --out gallery/spin-05.png
	./$(BIN) --preset sgra --spin 0.9   --inclination 80 --fov 26 --disc-outer 22 \
		--width 1280 --height 720 --spp 2 --out gallery/spin-09.png
	./$(BIN) --preset sgra --spin 0.998 --inclination 80 --fov 26 --disc-outer 22 \
		--width 1280 --height 720 --spp 2 --out gallery/spin-0998.png
	./$(BIN) --preset sgra --retrograde --inclination 80 --fov 30 \
		--width 1920 --height 1080 --spp 2 --out gallery/retrograde.png
	./$(BIN) --preset sgra --inclination 85 --fov 13 \
		--width 3840 --height 2160 --spp 2 --out gallery/photon-ring.png
	./$(BIN) --preset sgra --no-disc --fov 42 --stars 500000 \
		--width 1920 --height 1080 --spp 3 --out gallery/lensed-starfield.png
	./$(BIN) --preset sgra --no-disc --distance 40 --fov 110 --stars 500000 \
		--width 1600 --height 900 --spp 3 --out gallery/wide-field.png
	./$(BIN) --preset sgra --sun --no-disc --stars 400000 \
		--tonemap log --log-decades 13 --glare 0 \
		--width 2000 --height 875 --spp 3 --out gallery/sun-comparison.png
	./$(BIN) --preset sgra --sun --width 1920 --height 850 --spp 2 \
		--out gallery/sun-with-disc.png
	./$(BIN) --preset stellar --sun --width 1600 --height 900 --spp 2 \
		--out gallery/sun-vs-stellar.png

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
