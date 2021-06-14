NAME = main
DIR = ./MAIN
OBJFILEDIR = ../build/OFILES
MODFILEDIR = ../build/MODFILES
CPU = exe
NAMEFUL=../bin/GRPYlinked.$(CPU)
FC = gfortran
FFLAGS  = -O2 -c -ffpe-summary='none'
LD = gfortran
LDFLAGS = -O2 -ffpe-summary='none' 
LIBRARIES =   
LIBS = 
INCLUDES = -J$(MODFILEDIR) -I$(OBJFILEDIR)
FSRC =  TENSORS/modules_tensors.f \
        HYDRO/modules_hydro.f \
        CONFIGURATION/modules_config.f \
        FORCES/modules_force.f \
        STEPPER/modules_stepper.f \
		LAPACK/lapblas_double_excerpts.f \
        CONFIGURATION/initial_conf.f \
        TENSORS/tensors.f \
        HYDRO/GRPY_matrices.f \
        FORCES/forces.f \
        COLLISIONS/collisions.f \
        IO/IO.f \
        STEPPER/stepper_euler.f

OBJECTS := $(notdir $(FSRC:.f=.obj))
OBJECTSFULL := $(patsubst %,$(OBJFILEDIR)/%,$(notdir $(FSRC:.f=.obj)))
VPATH := $(subst $() $(),:,$(subst /,,$(dir $(FSRC))))

.PHONY: all clean

%.obj: %.f
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $(OBJFILEDIR)/$@

all: folders $(OBJECTS) $(DIR)/$(NAME).f
	$(LD) -o $(NAMEFUL) $(LDFLAGS) $(OBJECTSFULL) $(INCLUDES) $(LIBS) $(DIR)/$(NAME).f $(LIBRARIES)

clean:
	@rm -f $(NAMEFUL) $(OBJFILEDIR)/*.obj $(MODFILEDIR)/*.mod

folders:
	mkdir -p $(OBJFILEDIR) $(MODFILEDIR) ../bin

print-%  :
	@echo $* = $($*)

