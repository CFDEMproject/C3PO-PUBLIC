# C3PO CSV_interface multiple-machine Makefile

SHELL = /bin/bash

# Definitions
ROOT =	c3po_csv
EXE =	$(ROOT)
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)
COREDIR = $(C3PO_SRC_DIR)/interface_CSV
# Package variables
# List of all targets

help:
	@echo ''
	@echo '***C3PO CSV_interface MAKE SYSTEM***'
	@echo 'make clean-all           delete all object files'
	@echo 'make makelib             create Makefile.lib for static library build'
	@echo 'make makeshlib           create Makefile.shlib for shared library build'
	@echo 'make -f Makefile.lib machine      build C3PO CSV_interface as static library for machine'
	@echo 'make -f Makefile.shlib machine    build C3PO CSV_interface as shared library for machine'
	@echo 'make machine             build C3PO CSV_interface where machine is one of:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''

# Build the code

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
#	@$(SHELL) Make.sh style
#	@$(SHELL) Make.sh models
	@cp MAKE/Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" "SHFLAGS =" \
	  "EXE = ../$(EXE)" ../$(EXE)


# build thirdParty
#thirdParty:

# build test programs
#test:

# Remove all object files, static libraries, and libs of thirdparty
clean-all:
	rm -rf Obj_*  lib$(ROOT)*.a c3po_csv

# NOTIMPLEMENTED
makelib:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.lib

# NOTIMPLEMENTED
makeshlib:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.shlib

# NOTIMPLEMENTED - install C3PO shared lib and Python wrapper for Python usage
install-python:
	@python ../python/install.py

