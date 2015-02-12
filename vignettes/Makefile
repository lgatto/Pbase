all:
	make Pbase-data.html
	make mapping.html
	make ensucsc.html

R_HOME = `R RHOME`

Pbase-data.html: Pbase-data.Rmd
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("Pbase-data.Rmd"); render("Pbase-data.md")'

mapping.html: mapping.Rmd
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("mapping.Rmd"); render("mapping.md")'

ensucsc.html: ensucsc.Rmd
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("ensucsc.Rmd"); render("ensucsc.md")'

clean:
	rm -f *~
	rm -f .build.timestamp
	rm -f .Rhistory

allclean: clean
	rm -rf cache figure 
	rm -f Pbase-data.R mapping.R ensucsc.R
	rm -f *html
	rm -f *.md

