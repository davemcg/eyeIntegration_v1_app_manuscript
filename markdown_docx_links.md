link | info
---- | ---
https://rmarkdown.rstudio.com/articles_docx.html | overall guide on setting up rmarkdown/docx formatting
https://twitter.com/nicduquette/status/1074337228959014914?lang=en | opinionated guy says use Libertine font, kerning, ligatures, and old style numbers in Word
http://libertine-fonts.org | libertine fonts
https://community.rstudio.com/t/figure-caption-in-r-markdown/6951 | captioner for fig / table legends (bookdown/rmarkdown fig caps too limiting - can't have Fig 1 and Sup Fig 1, for example)
https://github.com/lierdakil/pandoc-crossref/issues/175 | caption for flextable (flextable is the ONLY easy option to get a latex like table in word). Have to use this annoying code bit to get captions for some reason.
https://stackoverflow.com/questions/52918716/authors-and-affiliations-in-the-yaml-of-rmarkdown | Custom lua code to get author affiliations in top of docx
https://github.com/pandoc/lua-filters | repo for above lua code (use  scholarly-metadata and author-info-blocks)
