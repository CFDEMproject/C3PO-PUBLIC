paths for additional libraries
CFDEM_ADD_LIB_PATHS = \
-L$(C3PO_QT5_LIB) \
-L$(C3PO_HDF5_LIB) \
-L$(CFDEM_POEMSLIB_PATH) \

# additional libraries to be linked to solvers
CFDEM_ADD_LIBS = \
-lhdf5_cpp \
-lQt5Core \
-lhdf5 \
-lc3po \


# additional static libraries to be linked to lagrangian library
CFDEM_ADD_STATICLIBS = \
-lmpi_cxx \

#################################################################
## SETTINGS FOR 3.0.x                                          ##
#################################################################
#----------------------------------------------------------------
# incompressible turbulence model settings
#----------------------------------------------------------------
# paths for incompressible turbulence models to use
CFDEM_ADD_INCOMPTURBMOD_PATHS = \
-I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
-I$(LIB_SRC)/TurbulenceModels/incompressible/lnInclude \

# libs for turbulence models to use
CFDEM_ADD_INCOMPTURBMOD_LIBS = \
-lturbulenceModels \
-lincompressibleTurbulenceModels \

#----------------------------------------------------------------
# compressible turbulence model settings
#----------------------------------------------------------------
# paths for compressible turbulence models to use
CFDEM_ADD_COMPTURBMOD_PATHS = \
-I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
-I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
-I$(LIB_SRC)/transportModels/compressible/lnInclude \
-I$(LIB_SRC)/thermophysicalModels/radiation/lnInclude \

# libs for turbulence models to use
CFDEM_ADD_COMPTURBMOD_LIBS = \
-lturbulenceModels \
-lcompressibleTurbulenceModels \
#################################################################
