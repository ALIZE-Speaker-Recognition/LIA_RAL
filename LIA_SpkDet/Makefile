-include ../config.txt
SUBDIRS=$(shell find . -mindepth 2 -maxdepth 2 -name Makefile)
PROGS=$(SUBDIRS:./%/Makefile=%)
export 

all: 	
	@($(MAKE) -C ../ spkdet)
clean:
	@($(MAKE) -C ../ spkdet_clean)
deps: 
	@($(MAKE) -C ../ deps)
install: 
	@($(MAKE) -C ../ spkdet_inst)	

all-rec: $(PROGS)
clean-rec: $(PROGS)
deps-rec: $(PROGS)
install-rec: $(PROGS)

.PHONY:$(PROGS)

$(PROGS):
	@$(MAKE) -C $@ $(MAKECMDGOALS) -$(MAKEFLAGS)
