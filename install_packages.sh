#!/bin/bash

## Default repo
cat >~/.Rprofile  <<EOL
local({r <- getOption("repos")
       r["CRAN"] <- "https://cloud.r-project.org"
       options(repos=r)
})
EOL

Rscript -e 'source("./src/setup_packages.R")'
