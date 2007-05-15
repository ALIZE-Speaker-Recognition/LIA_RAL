-include ../config.txt
SUBDIRS=$(shell find . -mindepth 2 -maxdepth 2 -name Makefile)
PROGS=$(SUBDIRS:./%/Makefile=%)
export 

all: 	
	@($(MAKE) -C ../ utils)
clean:
	@($(MAKE) -C ../ utils_clean)
deps: 
	@($(MAKE) -C ../ deps)
install: 
	@($(MAKE) -C ../ utils_inst)	

all-rec: $(PROGS)
clean-rec: $(PROGS)
deps-rec: $(PROGS)
install-rec: $(PROGS)

.PHONY:$(PROGS)

$(PROGS):
	@$(MAKE) -C $@ $(MAKECMDGOALS) -$(MAKEFLAGS)
