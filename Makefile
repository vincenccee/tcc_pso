PARAMS=-g -c -Wall -std=c++11
# -fopenmp

all: app

app: individual.o population.o main.o particle_swarm.o problem.o scenario.o ackley.o griewank.o rastring.o rosembrock.o schwefel.o sphere.o moving_peaks.o
	g++ individual.o population.o main.o particle_swarm.o problem.o scenario.o ackley.o griewank.o rastring.o rosembrock.o schwefel.o sphere.o moving_peaks.o -o app

main.o: main.cpp
	g++ $(PARAMS) main.cpp

particle_swarm.o: algorithms/particle_swarm.cpp
	g++ $(PARAMS) algorithms/particle_swarm.cpp

individual.o: individual.cpp
	g++ $(PARAMS) individual.cpp

population.o: population.cpp
	g++ $(PARAMS) population.cpp

problem.o: benchmarcks/problem.cpp
	g++ $(PARAMS) benchmarcks/problem.cpp

ackley.o: benchmarcks/statics/ackley.cpp
	g++ $(PARAMS) benchmarcks/statics/ackley.cpp

griewank.o: benchmarcks/statics/griewank.cpp
	g++ $(PARAMS) benchmarcks/statics/griewank.cpp

rastring.o: benchmarcks/statics/rastring.cpp
	g++ $(PARAMS) benchmarcks/statics/rastring.cpp

rosembrock.o: benchmarcks/statics/rosembrock.cpp
	g++ $(PARAMS) benchmarcks/statics/rosembrock.cpp

schwefel.o: benchmarcks/statics/schwefel.cpp
	g++ $(PARAMS) benchmarcks/statics/schwefel.cpp

sphere.o: benchmarcks/statics/sphere.cpp
	g++ $(PARAMS) benchmarcks/statics/sphere.cpp

moving_peaks.o: benchmarcks/dynamics/moving_peaks.cpp
	g++ $(PARAMS) benchmarcks/dynamics/moving_peaks.cpp

scenario.o: benchmarcks/scenario.cpp
	g++ $(PARAMS) benchmarcks/scenario.cpp

clean:
	rm *~ *.o *.eps *.dataset *.output app; clear;

god:
	make clean; make; ./app

try:
	make clean; make

# plot:
# 	gnuplot -e "file='saida.csv'" plot.plg

# to compile
# g++ -o test main.cpp population.cpp individual.cpp particle_swarm.cpp