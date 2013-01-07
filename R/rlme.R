###
### Test script preparing for main package code
###

#When location and scale are 1, it is Med-MAD in RPP prosedure
#When location and scale are 2, it is HL-Disp in RPP prosedure
location = scale = 2
tol = 10e-24

rlme = function(f, data, method="gr", print = FALSE, na.omit = TRUE, weight='wil', rprpair = 'hl-disp') {
  # Allow users to use upper case method names ("GR") as well as lowercase ("gr)
  method = tolower(method)
  
  # fit is the list that will be returned by rlme.
  # we start building the list here with the inputs to the function
  fit = list(
    formula = f,
    location = location,
    scale = scale
  )
  
  ###
  ### Remove rows with missing values
  ###
  if(na.omit == TRUE) {
    data = na.omit(data)
  }
  
  ###
  ### Extract data column names out of the formula
  ###
  
  # Extract the name of the response variable
  response_name = as.character(as.formula(f)[[2]])
  
  # Extract the random effect sections
  # this is a vector of ["(1 | school), "(1 | school:section)"]
  random_terms = str_extract_all(as.character(f)[[3]], "\\([0-9 a-zA-Z.:|]+)")
  random_terms = lapply(random_terms[[1]], function(str) substr(str, 2, nchar(str) - 1))
  
  covariate_names = setdiff(attributes(terms(f))$term.labels, random_terms)
  
  # Extract the school, section names
  school_name = section_name = ""
  
  levels = 0
  
  if(length(random_terms) == 2) {
    levels = 3
    
    school_name  = tail(strsplit(random_terms[[1]], "\\W")[[1]],1)
    section_name = tail(strsplit(random_terms[[2]], "\\W")[[1]],1)
    
    # Extract numbers of clusters, sections, and sample size in each section
    I = length(unique(factor(data[[school_name]])))   #number of uppest clusters
    sec = as.vector(sec_vec(data[[school_name]], data[[section_name]]))  #number of sections in each school
    ss = rhosch(data[[school_name]], data[[school_name]], data[[section_name]])$nvec  #sample size in each section
    
    fit$num.clusters = I
    fit$num.subclusters = sum(sec)
    
    
    mat = mat_vec(data[[school_name]],data[[section_name]])  #sample size for each section in school
    J = sum(sec)  	# The total number of sections
    n = sum(mat)			# The total sample size
    
    
    one = rep(1, length(data[[response_name]]))
    schsize = aggregate(one, by=list(data[[response_name]]), FUN=sum)[,2]  #this gives school sizes
    school1 = factor(rep(1:length(schsize),schsize))  #this generates school nos as I did alternative
    sss = aggregate(one, by=list(data[[school_name]], data[[section_name]]), FUN=sum)[,3]	#as alternative to my ss in main.begin. try
    section = factor(rep(1:length(ss),ss))
  }
  
  if(length(random_terms) == 1) {
    # This is the two-level case, which Dr. Bilgic is working on
    levels = 2
    
    school_name  = tail(strsplit(random_terms[[1]], "\\W")[[1]],1)
    
    I = length(unique(factor(data[[school_name]])))   #number of uppest clusters
    mat = mat_vec(data[[school_name]], rep(1, length(data[[school_name]])))
    n = sum(mat)
  }

  # Normalize covariates around their mean
  x = as.matrix(data[covariate_names])
  x = apply(x, 2, function(x){ x - mean(x) }) # Center the data around the mean
  
  # Build the response vector
  y = data[[response_name]]
  
  fit$num.obs = length(y)
  fit$y = y
  
  if(length(random_terms) == 0) {
    cat("You have entered an independent linear model.\n")
    cat("The packages Rfit or ww can be used to fit such models.\n")
    cat("Continuing using ww package.\n")
    
    fitw = wwest(x,y,print.tbl=T)
    
    return()
  }

  
  
  if(method == "reml" || method == "ml") {
    if(levels == 3) {
      REML=LM_est(x, y, dat=data.frame(y = y, x, school = data[[school_name]], section = data[[section_name]]), method= toupper(method))  #The function that collects estimates from lme()
    }
    else if(levels == 2) {
      REML=LM_est2(x, y, dat=data.frame(y = y, x, school = data[[school_name]]), method = toupper(method))  #The function that collects estimates from lme()
    }
    REMLb=REML$theta		#fixed effect estimates
    REMLs=REML$sigma		#var. cov estimates sch-sect-err
    REMLe=REML$ses		#se for the fixed estimates
    REML$ehat
    
    tvalue=REMLb/REMLe #this is a vector in the same order as lmer
    pvalue=2*pnorm(-abs(tvalue)) #this is a vector in the same order as lmer
    
    intracoeffs=c(REMLs[1]/sum(REMLs),(REMLs[1]+REMLs[2])/sum(REMLs)) #intra class correlations: intra-section and inter-section. this part is not reported in lmer. 
    
    
    fit$method = "REML"
    
    fit$ehat = REML$ehat
    
    fit$effect.err = REML$effect_err
    fit$effect.cluster = REML$effect_sch
    
    if(levels == 3) {
      fit$effect.subcluster = REML$effect_sec
    }
    
    fit$fixed.effects = data.frame(RowNames = c("(Intercept)", covariate_names), Estimate = REMLb, StdError = REMLe, tvalue = tvalue, pvalue = pvalue)
    fit$var.b = REML$varb
    fit$intra.class.correlations = intracoeffs
    fit$t.value = tvalue
    fit$p.value = pvalue
        
    fit$standr.lme = REML$standr.lme
  }
  else if(method == "jr") {
    if(levels == 3) {
      JR = JR_est(x,y,I,sec,mat,data[[school_name]], data[[section_name]], rprpair = rprpair)
    }
    else if(levels == 2) {
      JR = JR_est2(x, y, I, 1, mat, data[[school_name]], rep(1, length(data[[school_name]])), rprpair = rprpair)
    }
    
    JRb=JR$theta	#fixed effect estimates
    JRe=JR$ses		#se for the fixed estimates
    JRs=JR$sigma	#var. cov estimates using RPP (Groggel and Dubnicka's prosedure)
    
    tvalue=JRb/JRe #this is a vector in the same order as lmer
    pvalue=2*pnorm(-abs(tvalue)) #this is a vector in the same order as lmer
    
    intracoeffs=c(JRs[1]/sum(JRs),(JRs[1]+JRs[2])/sum(JRs)) #intra class correlations: intra-section and inter-section. this part is not reported in lmer. 
    
    fit$method = "JR"
    
    fit$ehat = JR$ehat
    
    fit$fixed.effects = data.frame(RowNames = c("(Intercept)", covariate_names), Estimate = JRb, StdError = JRe, tvalue = tvalue, pvalue = pvalue)
    
    fit$effect.err = JR$effect_err
    fit$effect.cluster = JR$effect_sch
    
    if(levels == 3) {
      fit$random.effects = data.frame(Groups = c(paste(school_name, ':', section_name, sep=''), school_name, "Residual"), Name = c("(Intercept)", "(Intercept)", ""), Variance = JRs)
      fit$effect.subcluster = JR$effect_sec
    }
    else if(levels == 2) {
      fit$random.effects = data.frame(Groups = c(school_name, "Residual"), Name = c("(Intercept)", ""), Variance = JRs)      
    }
    
    fit$intra.class.correlations = intracoeffs
    fit$var.b = JR$varb
    fit$t.value = tvalue
    fit$p.value = pvalue;
  }
  else if(method == "gr") {
    
    if(levels == 3) {
      GR=GR_est(x,y,I,sec,mat,data[[school_name]], data[[section_name]], rprpair = rprpair)
    }
    else if(levels == 2) {
      GR=GR_est2(x,y,I,rep(1, length(data[[school_name]])),mat,data[[school_name]], rep(1, length(data[[school_name]])), rprpair = rprpair)      
    }
    
    GRb=GR$theta  #this is fixed effects estimates (beta): intercept, sex, age 
    GRe=GR$ses #Std. Error for Fixed effects. So 't value' is calculated from GRb/GRe
    tvalue=GRb/GRe #this is a vector in the same order as lmer
    pvalue=2*pnorm(-abs(tvalue)) #this is a vector in the same order as lmer
    varb = GR$varb
    GRs=GR$sigma #variances of school, school:section, residual
    intracoeffs=c(GRs[1]/sum(GRs),(GRs[1]+GRs[2])/sum(GRs)) #intra class correlations: intra-section and inter-section. this part is not reported in lmer. 
    
    #store the followings in a vector when called from users. 
    GR$ehat #raw error 
    GR$ehats #independence error from the last weighted step
    GR$effect_err  #epsilon error
    GR$effect_sch  #random errors school
    GR$effect_sec  #random errors section
    #GR$ehat, effect_sch, effect_err, effect_sec and ehats will be stored too
    GR$xstar
    GR$ystar
    coll.stres <- stanresidgr(GR$xstar,GR$ystar,resid=GR$ehats,delta=.80,param=2,conf=.95) #standardized residual is calculated from this function. then following gets it.
    stresgr <- coll.stres$stanr
    
    fit$method = "GR"
    
    fit$ehat = GR$ehat
    fit$ehats = GR$ehats
    
    fit$xstar = GR$xstar
    fit$ystar = GR$ystar
    
    fit$fixed.effects = data.frame(RowNames = c("(Intercept)", covariate_names), Estimate = GRb, StdError = GRe, tvalue = tvalue, pvalue = pvalue)
    
    fit$effect.err = GR$effect_err
    fit$effect.cluster = GR$effect_sch
    
    if(levels == 3) {
      fit$random.effects = data.frame(Groups = c(paste(school_name, ':', section_name, sep=''), school_name, "Residual"), Name = c("(Intercept)", "(Intercept)", ""), Variance = GRs)
      fit$effect.subcluster = GR$effect_sec
    }
    else if(levels == 2) {
      fit$random.effects = data.frame(Groups = c(school_name, "Residual"), Name = c("(Intercept)", ""), Variance = GRs)    
    }
    
    fit$standard.residual = stresgr
    fit$intra.class.correlations = intracoeffs
    fit$var.b = varb
    fit$t.value = tvalue
    fit$p.value = pvalue;
  }
  else if(method == "geer") {
    tol=tol=1.1e-25  #need tolerance for nontrivial numerical results
    
    if(levels == 3) {
      GEER = GEER_est(x,y,I,sec,mat,data[[school_name]], data[[section_name]], weight = weight, rprpair = rprpair)
    }
    else if(levels == 2) {
      GEER = GEER_est2(x,y,I,rep(1, length(data[[school_name]])),mat,data[[school_name]], rep(1, length(data[[school_name]])), weight = weight, rprpair = rprpair)
    }
    
    GEERb=GEER$theta
    GEERs=GEER$sigma
    GEER_e2=GEER$ses_AP
    #GEER_e3=GEER$ses_CS
    
    intracoeffs=c(GEERs[1]/sum(GEERs),(GEERs[1]+GEERs[2])/sum(GEERs)) #intra class correlations: intra-section and inter-section. this part is not reported in lmer. 
    
    GEERe=GEER$ses_AP #Std. Error for Fixed effects. So 't value' is calculated from GRb/GRe
    tvalue=GEERb/GEERe #this is a vector in the same order as lmer
    pvalue=2*pnorm(-abs(tvalue)) #this is a vector in the same order as lmer
    
    fit$method = "GEER"
    
    fit$ehat = GEER$ehat
    
    fit$fixed.effects = data.frame(RowNames = c("(Intercept)", covariate_names), Estimate = GEERb, StdError = GEERe, tvalue = tvalue, pvalue = pvalue)
    
    fit$effect.err = GEER$effect_err
    fit$effect.cluster = GEER$effect_sch
    
    if(levels == 3) {
      fit$random.effects = data.frame(Groups = c(paste(school_name, ':', section_name, sep=''), school_name, "Residual"), Name = c("(Intercept)", "(Intercept)", ""), Variance = GEERs)
      fit$effect.subcluster = GEER$effect_sec
    }
    else if(levels == 2) {
      fit$random.effects = data.frame(Groups = c(school_name, "Residual"), Name = c("(Intercept)", ""), Variance = GEERs)      
    }
    
    fit$intra.class.correlations = intracoeffs
    
    fit$var.b = GEER$varb
    fit$t.value = tvalue
    fit$p.value = pvalue;
  }
  
  
  
  class(fit) = "rlme"
  
  if(print == TRUE) {
    summary(fit)
  }
  
  return(fit)
  
}

summary.rlme <- function(object, ...) {
  fit = object
  cat("Linear mixed model fit by ", fit$method, "\n")
  
  cat("Formula: ", deparse(fit$formula), "\n")
  
  ###
  ### Print out random effects table
  ###
  
  if('random.effects' %in% attributes(fit)$names) {
    cat("Random effects:\n")
    random.effects = fit$random.effects
    names(random.effects) = c("Groups", "Name", "Variance")
    print(random.effects, row.names = FALSE, right=FALSE)
  }
  
  ###
  ### Print out number of observations, clusters, and subclusters (if available)
  ###
  cat("\nNumber of obs:\n")
  cat(fit$num.obs)
  
  if('num.clusters' %in% attributes(fit)$names && 'num.subclusters' %in% attributes(fit)$names) {
    cat(" observations,", fit$num.clusters, "clusters,", fit$num.subclusters, "subclusters")
  }
  
  cat("\n")
  
  ###
  ### Print out fixed effects table
  ###
  cat("\nFixed effects:\n")
  
  fixed.effects = fit$fixed.effects
  names(fixed.effects) = c("", "Estimate", "Std. Error", "t value", "p value")
  
  print(fixed.effects, row.names = FALSE, right=FALSE)
  
  ###
  ### print out intra-class correlation coefficients
  ###
  
  if("intra.class.correlations" %in% attributes(fit)$names) {
    cat("\nIntra-class correlation coefficients\n")
    
    intra.coeffs = data.frame(names = c("intra-cluster", "intra-subcluster"), Estimates = fit$intra.class.correlations)
    names(intra.coeffs) = c("", "Esimates")
    
    print(intra.coeffs, row.names = FALSE, right = FALSE)
  }
  
  
  ###
  ### print out correlation of fixed effects table
  ### 
  
  cat("\ncov-var (fixed effects)\n")
  
  print(fit$var.b)
}

###
### RLME S3 Class Plotting Function
###
plot.rlme = function(x, ...) {
  if(x$method == "GR") {
    getgrstplot(x)
  }
  else if(x$method == "REML") {
    getlmestplot(x)
  }
  else {
    cat("rlme only supports plotting fits from GR and REML methods\n")
  }
}
