bin = add_scatter_factors

CXX = g++
CXXFLAGS = -std=c++14

$(bin) : main.o
	$(CXX) $(CXXFLAGS) -o $(bin) main.o

main.o : main.cpp

.PHONY : clean

clean :
	-rm -v $(bin) main.o
