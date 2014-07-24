
all:
	make data
	make mapping

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

data: 
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("Pbase-data.Rmd"); render("Pbase-data.md")'

mapping: 
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("mapping.Rmd"); render("mapping.md")'

clean: 
	rm -f *~
	rm -f .build.timestamp
	rm -f .Rhistory

allclean: clean
	rm -f *html
	rm -f *.md

