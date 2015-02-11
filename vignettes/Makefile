
all:
	make Pbase-data.html
	make mapping.html
	make ensucsc

R_HOME = `R RHOME`

.PHONY = ensucsc

Pbase-data.html: Pbase-data.Rmd
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("Pbase-data.Rmd"); render("Pbase-data.md")'

mapping.html: mapping.Rmd
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("mapping.Rmd"); render("mapping.md")'

ensucsc:
	cp vigsrc/ensucsc.html .
	cp vigsrc/ensucsc.R .

clean:
	rm -f *~
	rm -f .build.timestamp
	rm -f .Rhistory

allclean: clean
	rm -rf cache figure mapping_cache
	rm -f Pbase-data.R mapping.R
	rm -f *html
	rm -f *.md

