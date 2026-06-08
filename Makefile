CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

# Object files
OBJS = main.o matrix.o vector.o householder.o hessenberg.o

# Executable
TARGET = compute-eigenvalues

# Build
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Compile .cpp to .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean
