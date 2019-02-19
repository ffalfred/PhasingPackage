#' Calls phasing_class.py file
#'
#' This function only calls the phasing_class.py file, so it works like the python file phasing_trios.py
#'
#' @param infile Path to the input file, number of min_snps, output path
#' @return BED files
#'
#' @export
function_phasing<-function(python_path,input,min,output){
  reticulate::use_python(python_path)
  path_phasing<-system.file("python", package = "PhasingPackage")
  print(path_phasing)
  phasing_class<-reticulate::import_from_path("phasing_class",path=path_phasing)
  end<-phasing_class$phasing_class_tool(input,min,output)
  end$phasing() 
}
