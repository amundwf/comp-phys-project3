all: compile run

compile:
	g++ -o main.out main.cpp utils.cpp solver.cpp planet.cpp -larmadillo

run:
	./main.out


all_test: compile_test run_test

compile_test:
	g++ -o main.out main.cpp utils.cpp -larmadillo

run_test:
	./main.out
