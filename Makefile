CC = mpicxx
FLAGS = -D_GNU_SOURCE -Isrc/vectorclass -Isrc/hdf5 -Wall -Wextra -ansi -pedantic -std=c++11
SOURCE = $(wildcard src/*.cpp)
OBJECTS = $(SOURCE:src/%.cpp=obj/%.o)

STANDARD := -D_GNU_SOURCE -Wall -Wextra -ansi
INCLUDE  := -Isrc/vectorclass -Isrc/hdf5
LIB = -L lib_64/lib_silo lib_64/lib_silo/libsiloh5.a lib_64/lib_silo/libhdf5_hl.a lib_64/lib_silo/libhdf5.a lib_64/lib_silo/libsz.a lib_64/lib_silo/libz.a lib_64/lib_gsl/libgsl.a -ldl -lm

bin/aikef_mpi: $(OBJECTS)
	$(CC) $(INCLUDE) $(OBJECTS) -o bin/aikef_mpi $(STANDARD) $(LIB) -no-pie

obj/%.o: src/%.cpp 
	$(CC) $(FLAGS) -c $< -o $@

clean: 
	rm -f bin/aikef_mpi
	rm -f bin/process*
	rm -f obj/*

rmstate:
	rm $(wildcard bin/process_*log) 
	rm -r bin/State
	rm -r bin/Last_State


rmdata:
	rm -r  data/lineout
	mkdir data/lineout

	rm -r  data/particle_detector
	mkdir data/particle_detector

	rm -r  data/particle_tracks
	mkdir data/particle_tracks

	rm -r  data/uniform_output
	mkdir data/uniform_output

	rm -r  data/silo
	mkdir data/silo

	rm -r  data/silo_3D
	mkdir data/silo_3D

	rm -r  data/trajectories
	mkdir data/trajectories

rmdir:
	rm -r bin data obj

mkdir:
	-mkdir bin
	-mkdir data
	-mkdir obj
	-mkdir data/lineout
	-mkdir data/particle_detector
	-mkdir data/particle_tracks
	-mkdir data/uniform_output
	-mkdir data/silo
	-mkdir data/silo_3D
	-mkdir data/trajectories
