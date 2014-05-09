CC := g++
SRCDIR := src
BUILDDIR := build
CFLAGS := -Wall -Wextra -Wl,--no-as-needed -pedantic -march=native -mtune=native -O3 -std=c++11 -pthread
TARGET := tRayce
 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
 
$(TARGET): $(OBJECTS)
	@echo " Linking..."; $(CC) $(OBJECTS) -o $(TARGET)
 
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " CC $<"; $(CC) $(CFLAGS) -c -o $@ $<
 
clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET)
 
.PHONY: clean
