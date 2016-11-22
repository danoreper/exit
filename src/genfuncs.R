## General functions which may be used in multiple contexts
## 
## Author: doreper
###############################################################################


getPairOrientation <- function(allProps)
{
    alg = allProps$alg
    print(alg)
    orientation = allProps$pairOrientation[[alg]]
}

convertStrandToInt <- function(thestrand) 
{
    if(thestrand=="+")
    {
        i=2
    } else 
    {
        i=1
    }
    return(i)
}

getRevComplementPos <- function(refGenomeLen, pos)
{
    return(refGenomeLen - pos + 1)
}

reverseComplementStrings <- function(strings) 
{
    return(as.character(reverseComplement(DNAStringSet(strings))))
}

DNA2AA <- function(DNASeq) 
{

    s1 =  1+ nchar(DNASeq)%%3
    DNASeq = subseq(DNASeq, s1,nchar(DNASeq))
    nc = nchar(DNASeq)
    Dtriples = substring(DNASeq, seq(1, nc, by=3), seq(3, nc, 3))
    aanew = GENETIC_CODE[Dtriples];
    aanew[which(is.na(aanew))] = '#'
    dna = paste(aanew, collapse="")
}	

getAlgOptsString <- function(alg, mode)
{
    astring = ""
    optsForMode = allProps[[alg]][[mode]]
    keys = names(optsForMode$valued)
    for(key in keys)
    {
        value = optsForMode$valued[[key]]
        
        key = gsub("^prop_", replacement = "", x = key)

        astring = paste(astring, paste(" -", key, " ", value, sep=""), sep="")
    }
    for (binarySwitch in optsForMode$binary)
    {
        binarySwitch = gsub("^prop_", replacement = "", x = binarySwitch)

        astring = paste("-",binarySwitch, " ", astring, sep="")
    }
    print(astring)
    return(astring)
}

## improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
        fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

## shorthand
lsos <- function(..., n=10)
{
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}


seqLogo2 <- function(apwm) 
{
    seqLogo(t(t(apwm) * 1/colSums(apwm)))
}

writePostProcessedCSV <- function(aframe, fname) 
{
    write.table(aframe, 
                file.path(allProps$outputDir, allProps$postProcessDir, fname),
                sep="\t", quote=T,row.names=F)
}

##nearest neighbor median, for purposes of data.table, which cant deal with plain old median
##http://stackoverflow.com/questions/12125364/why-does-median-trip-up-data-table-integer-versus-double
getMedian <- function(numFusions) 
{
    unname(quantile(numFusions,.5, na.rm=T, type=3))
}

addRecord <- function(key, value, d)
{
    d[["statType"]]= append(d[["statType"]], key)
    d[["value"]]= append(d[["value"]], value)
    return(d)
}

getReadFile <- function(experimentID, run, readIndex=1)
{

    if(run=="")
    {
        return(NA)
    }
    return(file.path(allProps$dataDir, allProps$seqDir, 
                     paste("EXIT_",experimentID,"_",readIndex,".reads", sep="")))
}


buildAllJackpotTags <- function()
{
    x=(permutations(4, 7, v = c("A","C","G","T"), repeats.allowed=TRUE))
    jackpotTags=apply(x,1,paste,collapse="")
    return(jackpotTags)
}

allJackpotTags = factor(buildAllJackpotTags())
jackpotHash = data.table(jackpotTag = allJackpotTags, index=seq(1,length(allJackpotTags)))
setkey(jackpotHash, jackpotTag)

getValidExperimentalIds <- function(experimentData)
{
    validExperimentalIds = allProps$includedExperiments
    print(validExperimentalIds)
    if(-1 %in% validExperimentalIds)
    {
        print("inside valid experimental ids")
        validExperimentalIds = experimentData$experimentID
    }
    return(validExperimentalIds)
}

getExperimentData <- function()  
{
    mediumToInitial = data.table(medium = c("neither", "sucrose", "betalactam"), initial=c("h","s","c"))
    setkey(mediumToInitial, medium)
    data = (data.table(read.table(file.path(allProps$dataDir, allProps$seqDir, allProps$experimentDataFile), header=T, sep="\t")))
    validExperimentalIds = getValidExperimentalIds(data)
    
    setkey(data,"experimentID")
    data = data[J(experimentalID=validExperimentalIds),]
    
    data = data[data$run!=""]
    data$run = as.factor(data$run)
    data$fullname = as.factor(paste(data$replicate, data$medium, data$organ,sep="_"))
    data$experimentID = as.factor(data$experimentID)
    data$lane         = as.factor(data$lane)
    data$medium       = as.factor(data$medium)
    data$organ        = as.factor(data$organ)
    data$replicate    = as.factor(data$replicate)
    setkey(data,"experimentID")
    return(data)
}

