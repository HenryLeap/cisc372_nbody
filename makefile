# Parallelised by Samhain Ackerman and Henry Leap

FLAGS= -DDEBUG
LOCAL=
LIBS= -lm
ALWAYS_REBUILD=makefile

nbody: nbody.o compute.o
	nvcc $(FLAGS) $(LOCAL) $^ -o $@ $(LIBS)
nbody.o: nbody.cu planets.h config.h vector.h debug.h $(ALWAYS_REBUILD)
	nvcc $(FLAGS) $(LOCAL) -c $<
compute.o: compute.cu config.h vector.h debug.h $(ALWAYS_REBUILD)
	nvcc $(FLAGS) $(LOCAL) -c $<
clean:
	rm -f *.o nbody
