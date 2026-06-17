CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

SRCDIR = src
BENCHDIR = bench
BUILDDIR = build
INCLUDEDIR = include

SRCS = matrix.cpp vector.cpp householder.cpp hessenberg.cpp qr.cpp tridiagonal.cpp
LIB_OBJS = $(patsubst %.cpp,$(BUILDDIR)/%.o,$(SRCS))
MAIN_OBJ = $(BUILDDIR)/main.o
BENCH_OBJ = $(BUILDDIR)/benchmark.o
REACTOR_OBJ = $(BUILDDIR)/reactor.o

TARGET = $(BUILDDIR)/compute-eigenvalues
BENCHMARK = $(BUILDDIR)/benchmark
REACTOR = $(BUILDDIR)/reactor-network

all: $(TARGET)

$(TARGET): $(MAIN_OBJ) $(LIB_OBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BENCHMARK): $(BENCH_OBJ) $(LIB_OBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(REACTOR): $(REACTOR_OBJ) $(LIB_OBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR)/benchmark.o: $(BENCHDIR)/benchmark.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

benchmark: $(BENCHMARK)

reactor: $(REACTOR)

clean:
	rm -rf $(BUILDDIR)

.PHONY: clean all benchmark reactor
