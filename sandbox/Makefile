%.md: %.Rmd
	Rscript -e 'require(knitr);  knit("$^")'

%.html: %.md
	Rscript -e 'require(rmarkdown); render("$^")'


