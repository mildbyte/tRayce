CC := g++
SRCDIR := src
BUILDDIR := build
CFLAGS := -Wall -Wextra -pedantic -march=native -mtune=native -O3 -std=c++11
TARGET := tRayce
LFLAGS := -Wl,--no-as-needed -lpthread
 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
 
$(TARGET): $(OBJECTS)
	@echo " Linking..."; $(CC) $(OBJECTS) -o $(TARGET) $(LFLAGS)
 
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " CC $<"; $(CC) $(CFLAGS) -c -o $@ $<
 
clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET)
 
.PHONY: clean
