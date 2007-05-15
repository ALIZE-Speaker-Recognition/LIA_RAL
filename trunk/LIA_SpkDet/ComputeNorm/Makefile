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
LIB_LIA=$(LIA_SpkTools_DIR)/libliatools.a

###############################
# LIBS & INC DIR
###############################
DEPS_DIR=$(dir $(DEPS))
DEPS_INC=$(DEPS_DIR:%=-I%)
LIBSDIR = -L$(ALIZE_LIB) -L$(LIA_SpkTools_DIR) -L$(USER_LIB_DIR)
INCLDIR = -I$(SRC_DIR) -I$(ALIZE_DIR)/src -I$(LIA_SpkTools_DIR)/src $(DEPS_INC) 

LIBS = -lliatools -lalize -lm 
ifeq ($(THREAD),true)
LIBS +=-lpthread
CXXFLAGS+=-DTHREAD
endif

###############################
# Find executables
###############################
MAIN=$(wildcard $(SRC_DIR)/*Main.cpp)
EXEC=$(MAIN:$(SRC_DIR)/%Main.cpp=%)$(EXE_EXT)

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
install-rec: $(EXEC) 
	@echo "----- Installing $^-------------";
	cp $(EXEC) $(USER_BIN_DIR)
clean: clean-rec
clean-rec:
	@echo "----- Cleaning -------------";
	@rm -rf $(SRC_DIR)/*.o $(EXEC)
cleanlib:
	@echo "----- Clean LIA_SpkTools ----------";
	@$(MAKE) -C $(LIA_SpkTools_DIR) clean

###############################
# Objs, Exec, Lib targets
###############################
%.o:%.cpp
	@echo "----- Compiling $< -------------";
	$(CPP) $(INCLDIR) $(CXXFLAGS) -o $@ -c $<

$(EXEC): $(LIB_LIA) $(DEPS) $(OBJ)
	@echo "----- Linking $@ -------------";
	$(CPP) $(LDFLAGS) $(LIBSDIR) -o $@ $^ $(LIBS)
	@echo "---------- Done ----------------------";

###############################
# Dependancy on LIA_SpkTools
###############################
deps:
	$(MAKE) -C $(LIA_SpkTools_DIR)

.PHONY: clean-rec
# Protect deps from cleaning
.PRECIOUS: $(DEPS) $(DEPLIB) $(LIB_LIA)
