
R_HOME = `R RHOME`

ensucsc.html: ensucsc.Rmd
	"$(R_HOME)/bin/Rscript" -e 'library("knitr"); library("rmarkdown"); knit("ensucsc.Rmd"); purl("ensucsc.Rmd"); render("ensucsc.md")'
	rm ensucsc.md



