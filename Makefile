

#Compiler and compiler flags
CXX := mpic++
CXXFLAGS := -Wall -std=c++11  -O2 
LIB := -lconfig++
INC := -I /include/

#Directories
SRCDIR := src/simulation
BUILDDIR := build/simulation
TARGETDIR := bin
TESTDIR := test

#Target
TARGET = $(TARGETDIR)/DaMaSCUS-Simulator 

#Source files
SRCEXT := cpp
# COMMONSRC := $(shell find $(SRCDIR) -type f -name ./*.$(SRCEXT))
SIMSRC :=$(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
#$(shell find $(SRCDIR)/simulation -type f -name *.$(SRCEXT))

#Object files
SIMOBJ := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SIMSRC:.$(SRCEXT)=.o))

all: $(TARGET)

$(TARGET): $(SIMOBJ)
	$(CXX) $(CXXFLAGS) $(INC) $(LIB) -o $(TARGET) $^

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<



# 	# $(CXXFLAGS) $(INC) $(LDFLAGS) -o $@ $(SIMOBJ) $(LDLIBS)


# depend: .depend

# .depend: $(SIMOBJ)
# 	rm -f ./.depend
# 	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(SIMOBJ) $(TARGET)

.PHONY: clean

# dist-clean: clean
# 	rm -f *~ .depend

# include .depend