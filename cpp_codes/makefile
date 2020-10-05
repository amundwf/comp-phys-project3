all: compile run

compile:
	g++ -o main.out main.cpp utils.cpp -larmadillo #solver.cpp planet.cpp

run:
	./main.out


all_test: compile_test run_test

compile_test:
	g++ -o main.out main.cpp solver.cpp planet.cpp -larmadillo

run_test:
	./main.out
