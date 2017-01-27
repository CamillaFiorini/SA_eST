CXX = g++
CXXFLAGS = 

main: main.cpp mesh.o state.o flux.o utilities.o
	$(CXX) $(CXXFLAGS) -o main main.cpp mesh.o state.o flux.o utilities.o

mesh.o: mesh.hpp mesh.cpp
	$(CXX) $(CXXFLAGS) -c mesh.cpp

state.o: state.hpp state.cpp
	$(CXX) $(CXXFLAGS) -c state.cpp

utilities.o: utilities.hpp utilities.cpp
	$(CXX) $(CXXFLAGS) -c utilities.cpp

flux.o: flux.hpp flux.cpp state.o
	$(CXX) $(CXXFLAGS) -c flux.cpp

clean:
	$(RM) *.o

distclean:
	$(MAKE) clean
	$(RM) main

.PHONY: clean distclean