CC=/opt/local/bin/c++-mp-4.9

integrated: main.o sampler.o csv.o
#	$(CC) -std=c++11 -o integrated main.o sampler.o csv.o -lgsl
	$(CC) -std=c++11 -fopenmp -o integrated main.o sampler.o csv.o -lgsl
	cp integrated ../integrated_model

main.o: main.cpp sampler.hpp csv.hpp
	$(CC) -std=c++11 -c -o main.o main.cpp

sampler.o: sampler.cpp sampler.hpp
#	$(CC) -std=c++11 -c -o sampler.o sampler.cpp
	$(CC) -std=c++11 -fopenmp -c -o sampler.o sampler.cpp

csv.o: csv.cpp csv.hpp
	$(CC) -std=c++11 -c -o csv.o csv.cpp

clean:
	rm csv.o sampler.o main.o integrated