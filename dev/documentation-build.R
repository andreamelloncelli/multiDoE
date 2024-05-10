#!/usr/bin/env Rscript

# create package PDF document  --------------------------------------------

# update .Rd files (notice the Warning messages!!!)
devtools::document(roclets = c('rd', 'collate', 'namespace'))
# create the manual
devtools::build_manual()
# open the output file (or do it through the file browser)
system(paste0("xdg-open ../multiDoE_0.9.1.pdf"))
# in case it does not work you can debug it with:
devtools::check(manual=TRUE)


# get PDF help ------------------------------------------------------------

# update the docs
devtools::document(roclets = c('rd', 'collate', 'namespace'))
# install also the docs
devtools::install()
# select a function
fun <- "plotPareto"
# create a PDF document for that function only
help(fun, package = "multiDoE", help_type = "pdf")
# open the output file (or do it through the file browser)
system(paste0("xdg-open ",  fun, ".pdf"))
