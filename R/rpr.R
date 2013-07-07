# RPP Algorithm
rpr <- function(f, resd, data, rprpair='hl-disp') {
  random_terms = str_extract_all(as.character(f)[[3]], "\\([0-9 a-zA-Z.:|]+)")
  random_terms = lapply(random_terms[[1]], function(str) substr(str, 
                        2, nchar(str) - 1))
  
  rprpair = tolower(rprpair)
  location = scale = 2
  if (rprpair == "med-mad") {
      location = scale = 1
  }
  
  if (length(random_terms) == 2) {
    levels = 3
    school_name = tail(strsplit(random_terms[[1]], "\\W")[[1]], 
                       1)
    section_name = tail(strsplit(random_terms[[2]], "\\W")[[1]], 
                        1)
    I = length(unique(factor(data[[school_name]])))
    sec = as.vector(sec_vec(data[[school_name]], data[[section_name]]))
    mat = mat_vec(data[[school_name]], data[[section_name]])
    
    return(rprmeddis(I, sec, mat, resd, location, scale))
  }
  if (length(random_terms) == 1) {
    levels = 2
    school_name = tail(strsplit(random_terms[[1]], "\\W")[[1]], 
                       1)
    I = length(unique(factor(data[[school_name]])))
    mat = mat_vec(data[[school_name]], rep(1, length(data[[school_name]])))
    
    return(rprmeddis2(I, c(), mat, resd, location, scale))
  }
}