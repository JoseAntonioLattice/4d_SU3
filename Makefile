
FC=gfortran

TARGET = d=4_SU3_gauge.exe

SRC=src
BIN=bin

SOURCE= number2string_mod.f90 check_files_directories_mod.f90 create_files.f90 parameters.f90 data_types_observables.f90 matrix_operations.f90 pbc.f90 get_index_mod.f90 local_update_algorithms.f90 dynamics.f90 arrays.f90 starts.f90 statistics.f90 main.f90
OBJECT= $(patsubst %,$(BIN)/%, $(notdir $(SOURCE:.f90=.o)))


FFLAGS= -Wall -Wextra -fcheck=all -O0 -J$(BIN) -I$(BIN)

$(BIN)/$(TARGET): $(OBJECT)
	$(FC) -o $@ $^ -llapack

$(BIN)/%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: help run clean


run :
	@echo input_parameters.par | $(BIN)/$(TARGET)

help :
	@echo "src: $(SOURCE)"
	@echo "bin: $(OBJECT)"

clean:
	rm -f $(OBJECT) $(BIN)/$(TARGET)
