CXX = g++
CXXFLAGS = -std=c++17 -Wall -fopenmp -pthread -O3 -mcx16 -march=native

INCLUDE_PATH = -I../../../parlaylib/include/

all: cc scc

cc:
	$(CXX) $(CXXFLAGS) -o cc cc_main.cpp

scc:
	$(CXX) $(CXXFLAGS) -o scc scc_main.cpp

clean:
	rm cc scc
