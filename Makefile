
#Compiler and compiler flags
CXX := mpic++
CXXFLAGS := -Wall -std=c++11  -O2 
LIB := -lconfig++
INC := -I include

#Directories
SRCDIR := src
BUILDDIR := build
TARGETDIR := bin

#Targets
TARGETS = $(TARGETDIR)/DaMaSCUS-Simulator $(TARGETDIR)/DaMaSCUS-Analyzer

#Source files
SRCEXT := cpp
COMMONSRC := $(shell find $(SRCDIR) -maxdepth 1 -type f -name '*.$(SRCEXT)')
SIMSRC :=$(COMMONSRC) $(shell find $(SRCDIR)/simulation -type f -name *.$(SRCEXT))
ANASRC :=$(COMMONSRC) $(shell find $(SRCDIR)/analysis -type f -name *.$(SRCEXT))

#Object files
SIMOBJ := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SIMSRC:.$(SRCEXT)=.o))
ANAOBJ := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(ANASRC:.$(SRCEXT)=.o))

.PHONY: all simulator analyzer clean

all: $(TARGETS)

simulator: $(TARGETDIR)/DaMaSCUS-Simulator

analyzer: $(TARGETDIR)/DaMaSCUS-Analyzer

$(TARGETDIR)/DaMaSCUS-Simulator: $(SIMOBJ)
	$(CXX) $(CXXFLAGS) $(INC) $(LIB) -o $@ $^

$(TARGETDIR)/DaMaSCUS-Analyzer: $(ANAOBJ)
	$(CXX) $(CXXFLAGS) $(INC) $(LIB) -o $@ $^

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)/simulation/
	@mkdir -p $(BUILDDIR)/analysis/
	$(CXX) $(CXXFLAGS) $(INC) $(LIB) -o $@ -c $<

clean:
	rm -f $(SIMOBJ) $(ANAOBJ) $(TARGETS)

