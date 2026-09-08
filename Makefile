# Plain make build, for when cmake is not around.
#   make          build
#   make test     build and run the physics validation suite
#   make images   render the stills in gallery/ (slow: ~4 h, several are 4K)
#   make stars    render the neutron star stills in gallery/
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

src/main.o:     src/render.h src/jet.h src/disc.h src/scene.h src/image.h src/kerr.h src/geodesic.h src/constants.h
src/render.o:   src/render.h src/jet.h src/disc.h src/scene.h src/image.h src/kerr.h src/geodesic.h src/spectrum.h
src/spectrum.o: src/spectrum.h src/constants.h
src/image.o:    src/image.h src/spectrum.h
src/validate.o: src/kerr.h src/jet.h src/geodesic.h src/disc.h src/spectrum.h src/constants.h

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
	./$(BIN) --preset quasar --width 3840 --height 2160 --spp 2 --out gallery/quasar.png
# The full length of the jets, with the Milky Way behind them.
	./$(BIN) --preset quasar --distance 150 --fov 75 --jet-length 220 \
		--width 3840 --height 2160 --spp 2 --out gallery/quasar-jets.png
# Seen from 45 degrees the approaching jet is beamed towards us and the
# receding one away: the classic one-sided jet, produced by the ray tracer
# rather than drawn in.
	./$(BIN) --preset quasar --inclination 45 --fov 55 \
		--width 1920 --height 1440 --spp 2 --out gallery/quasar-beamed.png
	./$(BIN) --preset quasar --inclination 60 --fov 55 \
		--width 1920 --height 1440 --spp 2 --out gallery/quasar-face.png
# Same hole, same accretion rate, no spin - and therefore no jet at all, since
# Blandford-Znajek has nothing to extract.  The honest control image.
	./$(BIN) --preset quasar --spin 0.0 --width 1920 --height 1080 --spp 2 \
		--out gallery/quasar-spin0.png

# Neutron stars.  Mass and radius come from a TOV solve of the chosen equation
# of state, so the only inputs are the physics.
stars: $(BIN)
	@mkdir -p gallery
	./$(BIN) --neutron-star --log-decades 18 --stars 500000 --inclination 70 \
		--width 3840 --height 2160 --spp 2 --out gallery/ns-quiet.png
	./$(BIN) --neutron-star --ns-spin 400 --caps --cap-tilt 60 \
		--log-decades 18 --stars 500000 --inclination 70 \
		--width 3840 --height 2160 --spp 2 --out gallery/ns-pulsar.png
	./$(BIN) --neutron-star --ns-spin 700 --caps --cap-tilt 75 --cap-temp 4e6 \
		--log-decades 18 --stars 500000 --inclination 85 \
		--width 3840 --height 2160 --spp 2 --out gallery/ns-millisecond.png
	./$(BIN) --neutron-star --eos ms1 --log-decades 18 --stars 500000 --inclination 70 \
		--width 1920 --height 1080 --spp 2 --out gallery/ns-stiff.png
	./$(BIN) --neutron-star --ns-mass 2.05 --log-decades 18 --stars 500000 --inclination 70 \
		--width 1920 --height 1080 --spp 2 --out gallery/ns-heavy.png
	./$(BIN) --neutron-star --caps --cap-tilt 10 --cap-radius 30 \
		--log-decades 18 --stars 500000 --inclination 20 \
		--width 1920 --height 1080 --spp 2 --out gallery/ns-poleon.png
	./$(BIN) --neutron-star --log-decades 18 --stars 500000 --distance 25 --fov 120 \
		--width 1920 --height 1080 --spp 3 --out gallery/ns-closeup.png

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
	./$(BIN) --preset quasar --hotspot 1.6 --hotspot-size 0.4 --hotspot-contrast 250 \
		--orbits 3 --duration 10 --fps 24 --width 420 --height 300 --spp 1 \
		--out video/quasar-isco.apng
# Swinging the camera from nearly down the jet round to edge-on.  Nothing in
# the scene changes: the jets are identical in their own rest frame throughout.
# Everything you see happen is Doppler beaming turning off.
	./$(BIN) --preset quasar --inclination 8 --inclination-to 90 --fov 62 \
		--duration 10 --fps 15 --width 420 --height 315 --spp 1 \
		--out video/quasar-jet-sweep.apng

clean:
	rm -f $(OBJ) $(BIN)

.PHONY: all test images stars videos clean
