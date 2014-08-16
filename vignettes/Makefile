
all:
	make Pbase-data.html
	make mapping.html

.PHONY: data2 mapping2

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

Pbase-data.html: Pbase-data.Rmd
	## "$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("Pbase-data.Rmd"); render("Pbase-data.md")'
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); knit2html("Pbase-data.Rmd")'

mapping.html: mapping.Rmd
	## "$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("mapping.Rmd"); render("mapping.md")'
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); knit2html("mapping.Rmd");'

data2: 
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("Pbase-data.Rmd"); render("Pbase-data.md")'

mapping2: 
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("mapping.Rmd"); render("mapping.md")'

clean: 
	rm -f *~
	rm -f .build.timestamp
	rm -f .Rhistory

allclean: clean
	rm -f *html
	rm -f *.md

