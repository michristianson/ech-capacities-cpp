CC=gcc
CXX=g++
#Release Flags:
CXXFlags=-O3
#Profiling Flags:
#CXXFLAGS=-pg -Wall
#Debugging Flags:
#CXXFLAGS=-g -rdynamic -Wall -Wextra -Werror
LDLIBS=-lboost_filesystem #-lSegFault (for stack traces)

SRCS=ConvexToricDomains.cpp data.cpp lessthans.cpp
OBJS=ConvexToricDomains.o data.o lessthans.o

ConvexToricDomains: $(OBJS)
	$(CXX) -o ConvexToricDomains $(OBJS) $(LDLIBS)

ConvexToricDomains.o: ConvexToricDomains.cpp lessthans.hpp data.hpp
	$(CXX) $(CXXFLAGS) -c ConvexToricDomains.cpp

data.o: data.hpp data.cpp
	$(CXX) $(CXXFLAGS) -c data.cpp

lessthans.o: lessthans.hpp lessthans.cpp data.hpp
	$(CXX) $(CXXFLAGS) -c lessthans.cpp
