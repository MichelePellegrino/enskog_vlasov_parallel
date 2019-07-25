CXX = g++

STANDARD = -std=c++11
WARNINGS = -Wall
OPTIMIZATION = -g
CXXFLAGS = $(WARNINGS) $(STANDARD) $(OPTIMIZATION)

RESULT = ../libraries/libquadrature.a

SRC = domain.cpp mesh.cpp meshGenerators.cpp numerical_integration.cpp numerical_rule.cpp rk45.cpp
OBJS := $(SRC:%.cpp=%.o)
HEADERS = $(wildcard *.hpp)

domain.o: domain.hpp
mesh.o: domain.hpp meshGenerators.hpp mesh.hpp rk45.hpp
meshGenerators.o: meshGenerators.hpp domain.hpp rk45.hpp
numerical_integration.o: numerical_integration.hpp QuadratureRuleBase.hpp mesh.hpp
numerical_rule.o: numerical_rule.hpp QuadratureRule.hpp QuadratureRuleBase.hpp
rk45.o: rk45.hpp

.DEFAULT_GOAL = all

.PHONY: all clean distclean cleanobj

all: $(RESULT)

$(RESULT): $(OBJS) $(HEADERS)
	ar rs $@ $^

clean:
	$(RM) $(RESULT)

distclean:
	$(RM) $(RESULT)
	$(RM) *.o

cleanobj:
	$(RM) *.o
