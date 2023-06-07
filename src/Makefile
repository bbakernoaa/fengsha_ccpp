# Compiler
FC ?= gfortran
ifeq ($(FC), f77)  # override possible Make default
  FC := gfortran
endif
$(info FC setting: '$(FC)')

# Compile flags
DEBUG ?= 1
$(info DEBUG setting: '$(DEBUG)')
ifeq ($(DEBUG), 1)
  FCFLAGS := -g -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5
else ifeq ($(DEBUG), 0)
  FCFLAGS := -O3
else
  $(error invalid setting for DEBUG, should be 0 or 1 but is '$(DEBUG)')
endif

# Source objects
OBJS := \
  catchem_constants.o \
  catchem_config.o \
  dust_data_mod.o \
  dust_fengsha_mod.o

PROGRAM := dust_fengsha_mod.mod

# Targets
.PHONY: all clean
all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) $(FCFLAGS) $(LIBS) $(INC) -c $^

%.o: %.F90
	$(FC) $(FCFLAGS) $(INC) -c $<

clean:
	rm -f *.o *.mod $(PROGRAM)
