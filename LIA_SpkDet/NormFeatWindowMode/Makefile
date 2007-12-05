-include ../../config.txt
###############################
# Include your dependances here
###############################
DEPS=
###############################	

#***************************************
# DO NOT EDIT BELOW THIS LINES
#***************************************

###############################
#  Sources 
###############################
SRC_DIR=./src
SRC= $(wildcard $(SRC_DIR)/*.cpp)
OBJ= $(SRC:$(SRC_DIR)/%.cpp=$(SRC_DIR)/%.o)

###############################
# LIBS & INC DIR
###############################
LIBSDIR = -L$(ALIZE_DIR)/src -L$(LIA_SpkTools_DIR) -L$(USER_LIB_DIR)
INCLDIR = -I$(SRC_DIR) -I$(ALIZE_DIR)/src -I$(LIA_SpkTools_DIR)/src
DEPS_DIR=$(dir $(DEPS))
DEPS_INC=$(DEPS_DIR:%=-I%)
INCLDIR:=$(INCLDIR) $(DEPS_INC) 
LIB_LIA=libliatools.a
LIBS = -lliatools -lalize -lm 
DEPLIB=$(wildcard $(LIA_SpkTools_DIR)/src/*.cpp)
DEPLIB:=libliatools.a $(DEPLIB)

###############################
# Find executables
###############################
MAIN=$(wildcard $(SRC_DIR)/*Main.cpp)
EXEC=$(MAIN:$(SRC_DIR)/%Main.cpp=%.exe)

###############################
# Local Mode
###############################
all: deps $(EXEC)

###############################
#  Recursive targets
###############################
distrib: 
	@(cd ../../ && $(MAKE))
all-rec: $(EXEC) $(OBJ)
install: install-rec
install-rec: deps $(EXEC) 
	@echo "----- Installing $^-------------";
	cp $(EXEC) $(USER_BIN_DIR)
clean: clean-rec
clean-rec:
	@echo "----- Cleaning -------------";
	@rm -rf $(SRC_DIR)/*.o *.exe

###############################
# Objs, Exec, Lib targets
###############################
%.o:%.cpp
	@echo "----- Compiling $^ -------------";
	$(CPP) $(INCLDIR) $(CXXFLAGS) -o $@ -c $< 

%.exe:$(DEPS) $(OBJ) $(LIA_SpkTools_DIR)/$(LIB_LIA)
	@echo "----- Linking $@ -------------";
	$(CPP) $(LDFLAGS) $(LIBSDIR) -o $@  $^ $(LIBS)
	@echo "---------- Done ----------------------";

###############################
# Dependancy on LIA_SpkTools
###############################
deps: $(DEPLIB)
$(DEPLIB):
	@$(MAKE) -C $(LIA_SpkTools_DIR)/ $(LIB_LIA)

.PHONY: clean-rec
# Protect deps from cleaning
.PRECIOUS: $(DEPS) $(DEPLIB) $(LIB_LIA)
