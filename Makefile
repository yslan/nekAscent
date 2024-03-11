CXX ?= mpic++
FC ?= mpif77
CXXFLAGS ?=
LDFLAGS ?=
DEBUG ?= 0
UNDERSCORE ?= 1
ASCENT_DIR ?= 

########################## Don't touch what follows ###########################
ifeq ($(ASCENT_DIR),)
  $(error Specify ASCENT_DIR=<path to Ascent build>)
endif

MKFILEPATH := $(abspath $(lastword $(MAKEFILE_LIST)))
SRCROOT := $(realpath $(patsubst %/,%,$(dir $(MKFILEPATH))))
SRCDIR = $(SRCROOT)/src
FEXAMPLEDIR = $(SRCROOT)/examples/fortran
CEXAMPLEDIR = $(SRCROOT)/examples/cpp
BUILDROOT = $(SRCROOT)/build
INSTALLROOT = $(BUILDROOT)/install
ifneq ($(strip $(DESTDIR)),)
  INSTALLROOT = $(realpath $(DESTDIR))
endif

SRCS = $(wildcard $(SRCDIR)/*.cpp)
SRCOBJS = $(patsubst $(SRCROOT)/%.cpp,$(BUILDROOT)/%.o,$(SRCS))
FEXAMPLES = $(wildcard $(FEXAMPLEDIR)/*.f)
FEXAMPLEBINS = $(patsubst $(SRCROOT)/%.f,$(BUILDROOT)/%,$(FEXAMPLES))
CEXAMPLES = $(wildcard $(CEXAMPLEDIR)/*.cpp)
CEXAMPLEBINS = $(patsubst $(SRCROOT)/%.cpp,$(BUILDROOT)/%,$(CEXAMPLES))

LIB = $(BUILDROOT)/lib/libnekAscent.a

# See $(ASCENT_DIR)/share/ascent/ascent_config.mk for detailed linking info 
include $(ASCENT_DIR)/share/ascent/ascent_config.mk

# make sure to enable c++11 support (conduit's interface now requires it)
ASCENT_CXXFLAGS = -std=c++11
ASCENT_INCFLAGS = $(ASCENT_INCLUDE_FLAGS)
ASCENT_LNKFLAGS = $(ASCENT_LINK_RPATH) $(ASCENT_MPI_LIB_FLAGS)

ifneq ($(DEBUG),0)
  PP += -DDEBUG
  CXXFLAGS += -g
  FFLAGS += -g
else
  CXXFLAGS += -O2
  FFLAGS += -O2
endif
FFLAGS += -fdefault-real-8 -fdefault-double-8

INCFLAGS = -I$(SRCDIR)

CXXLDFLAGS = $(LDFLAGS)
FLDFLAGS = $(LDFLAGS)
CXXLDFLAGS += -lm 
FLDFLAGS += -lstdc++

ifneq ($(UNDERSCORE),0)
  PP += -DUNDERSCORE
endif

CCMD = $(CXX) $(CXXFLAGS) $(ASCENT_CXXFLAGS) $(INCFLAGS) $(ASCENT_INCFLAGS) $(PP)
FCMD = $(FC) $(FFLAGS) $(INCFLAGS)


.PHONY: all lib install examples link.txt clean

all: lib install examples link.txt

lib: $(SRCOBJS)
	@mkdir -p $(BUILDROOT)/lib
	@$(AR) cr $(LIB) $?
	@ranlib $(LIB)

install: lib
	@mkdir -p $(INSTALLROOT)/lib 2>/dev/null
	@cp -v $(LIB) $(INSTALLROOT)/lib 2>/dev/null
	@mkdir -p $(INSTALLROOT)/include 2>/dev/null
	@cp $(SRCDIR)/*.hpp $(INSTALLROOT)/include 2>/dev/null

examples: lib install $(CEXAMPLEBINS) $(FEXAMPLEBINS)

link.txt:
	@echo $(ASCENT_LNKFLAGS) > $(BUILDROOT)/link.txt
	@cp $(BUILDROOT)/link.txt $(INSTALLROOT) 2>/dev/null

clean:
	@$(RM) -rf $(BUILDROOT)

print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info)
	@true

$(BUILDROOT)/%.o: $(SRCROOT)/%.cpp
	$(CCMD) -c $< -o $@ 

$(BUILDROOT)/%.o: $(SRCROOT)/%.f
	$(FCMD) -c $< -o $@ 

# executable (examples)
$(BUILDROOT)/examples/cpp/%: $(BUILDROOT)/examples/cpp/%.o | lib install
	$(CCMD) $< -o $@ -L$(INSTALLROOT)/lib -lnekAscent $(ASCENT_LNKFLAGS) $(CXXLDFLAGS)

$(BUILDROOT)/examples/fortran/%: $(BUILDROOT)/examples/fortran/%.o | lib install
	$(FCMD) -o $@ $< -L$(INSTALLROOT)/lib -Wl,-Bstatic -lnekAscent -Wl,-Bdynamic $(ASCENT_LNKFLAGS) $(FLDFLAGS)

$(shell mkdir -p $(BUILDROOT)/examples)
$(shell mkdir -p $(BUILDROOT)/examples/cpp)
$(shell mkdir -p $(BUILDROOT)/examples/fortran)
$(shell mkdir -p $(BUILDROOT)/src)
