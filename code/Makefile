COMPILER   = gfortran
COMPFLAGS  = -c #-mcmodel=medium -i-dynamic
OBJECTS    = core_model.o core_evolution.o parameters.o parameters_core.o splat.o drspln.o mantle_evolution.o #stable_layer.o
EXECUTABLE = thermal_history


all : $(OBJECTS) $(EXECUTABLE).f90
	$(COMPILER) $(COMPFLAGS) $(EXECUTABLE).f90
	$(COMPILER) -o $(EXECUTABLE) $(EXECUTABLE).o $(OBJECTS)

#stable_layer.o : stable_layer.f90
#	$(COMPILER) $(COMPFLAGS) stable_layer.f90

mantle_evolution.o : mantle_evolution.f90 parameters_mantle.o
	$(COMPILER) $(COMPFLAGS) mantle_evolution.f90

core_evolution.o : core_evolution.f90 core_model.o parameters.o
	$(COMPILER) $(COMPFLAGS) core_evolution.f90

core_model.o : core_model.f90 parameters_core.o parameters.o 
	$(COMPILER) $(COMPFLAGS) core_model.f90

parameters_mantle.o : parameters_mantle.f90 parameters.o
	$(COMPILER) $(COMPFLAGS) parameters_mantle.f90

parameters_core.o : parameters_core.f90 
	$(COMPILER) $(COMPFLAGS) parameters_core.f90

parameters.o : parameters.f90
	$(COMPILER) $(COMPFLAGS) parameters.f90

splat.o : splat.f
	$(COMPILER) $(COMPFLAGS) splat.f
drspln.o : drspln.f
	$(COMPILER) $(COMPFLAGS) drspln.f

clean        :
	rm -rf *.o *.mod $(EXECUTABLE)

