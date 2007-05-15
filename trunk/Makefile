-include config.txt

ToDeleteDistrib = $(shell find . -iname .depend && find . -iname *.exe && find . -iname *.a && find . -type d -name objs && find . -type d -name CVS| xargs)
ToDeleteMrProper = $(shell find . -iname .depend && find . -iname *.exe && find . -iname *.a && find . -type d -name objs| xargs)
export

PROGS=$(LIA_SpkTools_DIR) $(LIA_SpkDet_DIR)  $(LIA_Utils_DIR)

all: tools spkdet utils 

tools:
	@echo "***********************************"
	@echo "**	 LIA SpkTools Package	**"
	@echo "***********************************"
	@$(MAKE) -C $(LIA_SpkTools_DIR) -$(MAKEFLAGS) all-rec

spkdet:
	@echo "***********************************"
	@echo "**	 LIA_SpkDet Package	**"
	@echo "***********************************"
	@$(MAKE) -C $(LIA_SpkDet_DIR) -$(MAKEFLAGS) all-rec

utils:	
	@echo "***********************************"
	@echo "**	 LIA_Utils Package	**"
	@echo "***********************************"
	@$(MAKE) -C $(LIA_Utils_DIR) -$(MAKEFLAGS) all-rec
	@echo "***********************************"

clean: spkdet_clean utils_clean tools_clean

spkdet_clean: 
	@echo "***********************************"
	@echo "Cleaning LIA_SpkDet"
	@echo "***********************************"
	@$(MAKE) -C $(LIA_SpkDet_DIR) -$(MAKEFLAGS) clean-rec
	@echo "***********************************"

utils_clean: 
	@echo "***********************************"
	@echo "Cleaning LIA_Utils"
	@echo "***********************************"
	@$(MAKE) -C $(LIA_Utils_DIR) -$(MAKEFLAGS) clean-rec
	@echo "***********************************"

tools_clean: 
	@echo "***********************************"
	@echo "Cleaning LIA_SpkTools"
	@echo "***********************************"
	@$(MAKE) -C $(LIA_SpkTools_DIR) -$(MAKEFLAGS) clean-rec
	@echo "***********************************"

deps:
	@echo "***********************************"
	@echo "Generating dependencies"
	@echo "***********************************"
	@$(MAKE) -C $(LIA_SpkTools_DIR) deps-rec
	@$(MAKE) -C $(LIA_SpkDet_DIR) deps-rec
	@$(MAKE) -C $(LIA_Utils_DIR) deps-rec
	@echo "***********************************"

install: spkdet_inst utils_inst tools_inst
	@echo "***********************************"
	@echo "Copying to your BIN dir"
	@echo "***********************************"
tools_inst:
	@$(MAKE) -C $(LIA_SpkTools_DIR) -$(MAKEFLAGS) install-rec
spkdet_inst:
	@$(MAKE) -C $(LIA_SpkDet_DIR) -$(MAKEFLAGS) install-rec
utils_inst:
	@$(MAKE) -C $(LIA_Utils_DIR) -$(MAKEFLAGS) install-rec

# Have to remove the man/ dir if it's ok with its location
distrib:
	@echo "***********************************"
	@echo "Let's make a Package: Generating $(PACK).tar.gz"
	@echo "***********************************"
	@cd .. && cp -rf LIA_RAL $(PACK) && cd $(PACK) && rm -rf $(ToDeleteDistrib) config.txt
	@cd .. && tar -zcf LIA_RAL/$(PACK).tar.gz $(PACK) && rm -rf $(PACK)
	@echo "Done"
	@echo "***********************************"

dist-clean:
	@echo "***********************************"
	@echo "Mr Proper Clean: Removing dependencies, objs, exec..."	
	@echo "***********************************"
	@rm -r $(ToDeleteMrProper)
	@echo "***********************************"
	
tags:	
	@echo "***********************************"
	@echo "Generating tags in $(TAGS_DIR)tags"
	@echo "***********************************"
	@$(TAGS) --excmd=number --c-types=pcdgstu `find . -name *.h | xargs`
	@echo "If you use Scite Editor, do a python tags2api.py tags > liaral.api"
	@echo "***********************************"
	
doxygen: 
	@echo "***********************************"
	@echo "Generating documentation with Doxygen in ./doc/doxygen"
	@echo "***********************************"
	@doxygen doxygen.cfg
	@echo "***********************************"
	
docs: 
	@echo "***********************************"
	@echo "**	 LIA ASR Package	**"
	@echo "***********************************"
	@echo "Generating documentation with pdflatex in ./doc/"
	@cd doc && $(MAKE) docs-rec

.PHONY: clean
