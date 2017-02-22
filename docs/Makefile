target : index.html HyDe-notes.html

index.html : HyDe-manual.Rmd
	Rscript -e "rmarkdown::render('HyDe-manual.Rmd')"
	@mv HyDe-manual.html index.html
	
HyDe-notes.html : HyDe-notes.Rmd
	Rscript -e "rmarkdown::render('HyDe-notes.Rmd')"