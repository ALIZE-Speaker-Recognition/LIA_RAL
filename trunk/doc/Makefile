#TEX= $(wildcard *.tex)
TEX=DOC_LIA_SpkDet.tex
TEX+=doc_SpkSeg.tex
TEX+=DOC_LIA_Utils.tex
AUX=$(TEX:%.tex=%.aux)
LOG=$(TEX:%.tex=%.log)
OUT=$(TEX:%.tex=%.out)
PDF= $(TEX:%.tex=%.pdf)

all: 
	@cd .. && $(MAKE) docs

docs-rec: $(PDF) clean

%.pdf: %.tex
	pdflatex $<
	
clean:
	rm -f $(AUX) $(LOG) $(OUT)
