
all:
	make data

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

data: 
	"$(R_HOME)/bin/Rscript" -e 'require("knitr"); require("rmarkdown"); knit("Pbase-data.Rmd"); render("Pbase-data.md")'

clean: 
	rm -f *~
	rm -f .build.timestamp
	rm -f .Rhistory

allclean: clean
	rm -f *html
	rm -f *.md

