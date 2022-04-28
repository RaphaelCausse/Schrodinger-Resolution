# Directories set up
OBJDIR =obj/
SRCDIR =src/

# Executable name
TARGET =prog

# Project set up, compiler flags and linker flags
CC =gcc -fsanitize=address
CFLAGS =-g -std=c99 -O3 -Wall -Iinclude
LFLAGS =-lgsl -lgslcblas -lm

# Files set up
SRC := $(wildcard $(SRCDIR)*.c)
OBJ := $(SRC:$(SRCDIR)%.c=$(OBJDIR)%.o)

# Compile binary and object files and build target
all: $(TARGET) objclean
	
$(TARGET): $(OBJ)
	$(CC) $(LFLAGS) $^ -o $@ 
	@echo "Compilation and linking completed !"

$(OBJDIR)%.o: $(SRCDIR)%.c
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean project directory
.PHONY: run clean
run:
	./$(TARGET)

clean:
	@rm -rf $(BINDIR) $(OBJDIR)
	@echo "Cleanup completed !"

objclean:
	@rm -rf $(OBJDIR)