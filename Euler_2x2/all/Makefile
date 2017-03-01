CXX = g++
CXXFLAGS = -std=c++11 -Wall

main: main.cpp mesh.o state.o roe.o roe_I.o roe_II.o godunov.o
	$(CXX) $(CXXFLAGS) -o main main.cpp mesh.o state.o roe.o roe_I.o roe_II.o godunov.o

mesh.o: mesh.hpp mesh.cpp
	$(CXX) $(CXXFLAGS) -c mesh.cpp

godunov.o: godunov.hpp godunov.cpp
	$(CXX) $(CXXFLAGS) -c godunov.cpp

roe.o: roe.hpp roe.cpp
	$(CXX) $(CXXFLAGS) -c roe.cpp

roe_II.o: roe_II.hpp roe_II.cpp roe.o
	$(CXX) $(CXXFLAGS) -c roe_II.cpp

roe_I.o: roe_I.hpp roe_I.cpp roe.o
	$(CXX) $(CXXFLAGS) -c roe_I.cpp

state.o: state.hpp state.cpp
	$(CXX) $(CXXFLAGS) -c state.cpp

clean:
	$(RM) *.o

distclean:
	$(MAKE) clean
	$(RM) main

.PHONY: clean distclean