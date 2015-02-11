
all:
	make Pbase-data.html
	make mapping.html

R_HOME = `R RHOME`

Pbase-data.html: Pbase-data.Rmd
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("Pbase-data.Rmd"); render("Pbase-data.md")'

mapping.html: mapping.Rmd
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("mapping.Rmd"); render("mapping.md")'

clean:
	rm -f *~
	rm -f .build.timestamp
	rm -f .Rhistory

allclean: clean
	rm -rf cache figure mapping_cache
	rm -f Pbase-data.R mapping.R
	rm -f *html
	rm -f *.md

