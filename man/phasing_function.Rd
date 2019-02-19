% Please edit documentation in R/load.R
\name{phasing_function}
\alias{phasing_function}
\title{Calls phasing_class.py file}
\usage{
  phasing_function(input,min,outpu)
}
\arguments{
  \item{input,min,outpu}{Path to the input file, number of min_snps, output path}
}
\value{
  A matrix of the infile
}
\description{
  This function only calls the phasing_class.py file, so it works like the python file phasing_trios.py
}