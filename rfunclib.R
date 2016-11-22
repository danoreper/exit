## Initialization of properties, loads all functions that will be used
## 
## Author: doreper
###############################################################################



library(R453Plus1Toolbox)
library(Biostrings)
library(yaml)
library(Rsamtools)
library(ggplot2)
library(plyr)
library(doMC)
library(ggplot2)
library(gtools)
library(epitools)
registerDoMC()
library("ShortRead")
library("qrqc")
library("data.table")
library("seqLogo")
library("data.table")


## loads params from command line
getAllPropsFromCL <-  function()
{
    yamldefault  = "../config/config.txt" 
    prop          = yaml::yaml.load_file(yamldefault)
        
    args= commandArgs(T)
    print(paste0("override arg is ",args))
    if(length(args)>0)
    {
        for(arg in args)
        {
            overrideYaml = arg
            print(arg)
            propOver = yaml::yaml.load_file(overrideYaml)
            prop = loadParamRecurseHelp(prop,propOver)
        }
    }
    prop = list2env(prop)
    return(prop)
}

loadParamRecurseHelp <- function(origlist, overlist)
{
    for(aname in ls(overlist))
    {
        if(aname %in% ls(origlist))
        {
            if(!(is.list(origlist[[aname]])))
            {
                origlist[[aname]] = overlist[[aname]]
            } else {
                origlist[[aname]] = loadParamRecurseHelp(origlist[[aname]], overlist[[aname]])
            }
        }
        else
        {
            warning(paste0("prop in override file doesn't exist in default: ", substitute(overlist),":", aname))
        }
    }
    return(origlist)
}


##intentionally in global scope for ease of typing as its everywhere
dat = function(x)
{
    file.path(prop$data, x)
}

fp = file.path



theargs = commandArgs(trailingOnly=TRUE)
if(length(theargs)==0)
{
    theargs = "../config/config.txt"
}
allProps = getAllPropsFromCL()
dir.create(file.path(allProps$outputDir), showWarnings = FALSE)

print("done loading")
for (aname in ls(allProps))
{
    print("*********")
    print(aname)
    print(allProps[[aname]])
}

source("postProcessFuncs.R")
source("alignFuncs.R")
source("genfuncs.R")
source("genomeRef.R")

