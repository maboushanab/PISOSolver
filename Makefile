# Compiler and linker settings
CC = g++
CFLAGS = -g -O0 -I./include/eigen-3.4.0 -Wall -Wextra -std=c++17 -I/usr/include/python3.8 -I/usr/include/python3.8  -Wno-unused-result -Wsign-compare -g -fdebug-prefix-map=/build/python3.8-YBWzqg/python3.8-3.8.10=. -specs=/usr/share/dpkg/no-pie-compile.specs -fstack-protector -Wformat -Werror=format-security  -DNDEBUG -g -fwrapv -Wall -I/path/to/boost/include -fPIE
LDFLAGS = -O0 -I/usr/include/python3.8 -I/usr/include/python3.8  -Wno-unused-result -Wsign-compare -g -fdebug-prefix-map=/build/python3.8-YBWzqg/python3.8-3.8.10=. -specs=/usr/share/dpkg/no-pie-compile.specs -fstack-protector -Wformat -Werror=format-security  -DNDEBUG -g -fwrapv -Wall -L/path/to/boost/lib -lboost_python38 -lpython3.8 -pie

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
