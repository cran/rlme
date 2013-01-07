sourcer = function(){
    source("sourcer.r")

    #3-level data generator
    source("errdist.r")
    source("schdist.r")
    source("secdist.r")
    source("twonestgen.r")

    #Rank-based functions
    source("ww1.7.r")
    source("jmake.r")
    source("bdmake.r")
    source("sigma12i.r")
    source("sigma12.r")
    source("wilstep.r")
    source("projcsp.r")
    source("dispvar.r")

#3level functions
    source("rprmeddis.r")
    source("rhoestimators.r")  #many functions for JR and GR
    source("sigymake.r")

    source("gee_func.r") 	#some functions for gee
    source("three_methods.r") 	#three big functions jr gr geer +reml

    source("stanresidgr.r") 	#standardized resid for gr

#2level functions
    source("rprmeddis2.r") #rpr
    source("three_methods2.r")  #jr gr geer
    source("sigymake2.r")       #one function
    source("rhoestimators2.r")  #many functions
    source("wilstep2.r")

#from locationR, McKean
     source("onesampwil.r")
     source("pairup.r")

    source("fitdvcov.r")  #gee TDBETAS. 1/1/2013 by Joe
    source("hbrwts_gr.r") #hbr weights for gee. this is modified hbr in ww.
    source("getgrstplot.r") #gr stand resid plot and qq plot
}
