CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

# Object files
OBJS = main.o matrix.o vector.o householder.o hessenberg.o qr.o
LIB_OBJS = matrix.o vector.o householder.o hessenberg.o qr.o

# Executables
TARGET = compute-eigenvalues
BENCHMARK = benchmark

# Build
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

benchmark: benchmark.o $(LIB_OBJS)
	$(CXX) $(CXXFLAGS) -o $(BENCHMARK) benchmark.o $(LIB_OBJS)

# Compile .cpp to .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) benchmark.o $(TARGET) $(BENCHMARK)

.PHONY: clean all benchmark
