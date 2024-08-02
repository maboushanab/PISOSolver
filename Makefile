# Compiler and linker settings
CC = g++
CFLAGS = -g -I./include/eigen-3.4.0 -Wall -Wextra -std=c++17 $(shell python3-config --cflags) -I/path/to/boost/include -fPIE
LDFLAGS = $(shell python3-config --ldflags) -L/path/to/boost/lib -lboost_python38 -lpython3.8 -pie

# Directories
SRCDIR = src
OBJDIR = obj
BINDIR = bin

# Source and object files
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

# Executable
EXEC = $(BINDIR)/program

# Rule to build the executable
$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $@

# Rule to build object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build artifacts
.PHONY: clean
clean:
	rm -rf $(OBJDIR)/*.o $(EXEC)
