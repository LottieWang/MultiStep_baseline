CXX = g++
# CXXFLAGS = -Wall -fopenmp -O3 -fcilkplus -DCILK -march=native
CXXFLAGS = -Wall -fopenmp -O3 -fcilkplus -DCILK -march=native

all: cc scc

cc:
	$(CXX) $(CXXFLAGS) -o cc cc_main.cpp

scc:
	$(CXX) $(CXXFLAGS) -o scc scc_main.cpp

clean:
	rm cc scc
