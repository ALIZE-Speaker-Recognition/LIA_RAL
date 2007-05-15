-include ../config.txt
###############################
# LIBS & INC DIR
###############################
LIBSDIR = -L$(ALIZE_LIB) -L$(LIA_SpkTools_DIR) -L$(USER_LIB_DIR)
INCLDIR = -I$(SRC_DIR) -I$(ALIZE_DIR)/src -I$(LIA_SpkTools_DIR)/src
LIB_LIA=libliatools.a
LIBS = -lliatools -lalize -lm 
DEPLIB=$(wildcard $(ALIZE_DIR)/src/*.cpp)
DEPLIB:=$(ALIZE_LIB)/libalize.a $(DEPLIB)

###############################
#  Sources 
###############################
SRC_DIR=./src
SRC= $(wildcard $(SRC_DIR)/*.cpp)
INC= $(wildcard $(SRC_DIR)/*.h)
OBJ= $(SRC:$(SRC_DIR)/%.cpp=$(SRC_DIR)/%.o)

ifeq ($(THREAD),true)
CXXFLAGS+=-DTHREAD
endif

###############################
# Local Mode
###############################
all: $(LIB_LIA)

###############################
#  Recursive targets
###############################
distrib: 
	@(cd ../ && $(MAKE))
all-rec: $(LIB_LIA) 
install: install-rec
install-rec: $(LIB_LIA)
	@echo "----- Installing  -------------";
	cp $(LIB_LIA) $(USER_LIB_DIR)
clean: clean-rec
clean-rec:
	@echo "----- Cleaning -------------";
	@rm -f $(SRC_DIR)/*.o $(LIB_LIA)

###############################
# Objs, Exec, Lib targets
###############################
$(LIB_LIA):$(OBJ) 
	@echo "---------- Building Lib ----------------------";
	ar cru $@ $^
	ranlib $@
	@echo "---------- Done ----------------------";

%.o:%.cpp
	@echo "----- Compiling $< -------------";
	$(CPP) $(INCLDIR) $(CXXFLAGS) -o $@ -c $< 
	
###############################
# Dependancy on ALIZE
###############################	
deps: $(DEPLIB)

$(DEPLIB):
	@$(MAKE) -C $(ALIZE_DIR) all	
	
.PHONY: clean-rec
# Protect deps from cleaning
.PRECIOUS: $(DEPS) $(DEPLIB) $(LIB_LIA)
