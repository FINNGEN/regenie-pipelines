#!/usr/bin/env Rscript

packs <- c("qqman","optparse","data.table","R.utils")

for (p in packs) {
  if( !require(p, character.only = T)) {
    print(p)
    install.packages( p,  repos = c(CRAN = "http://cran.r-project.org") )
  }
}
