## Functions called by postProcessAlignment for postprocessing and assembling 
## data frames of the aligned read results from .bam diles
## 
###############################################################################

getReadFrameForPooledFile <- function(bamFile, trimFile, experimentID, startId)
{	
    print(paste("started loading frame:", Sys.time()))
    
    x=scanBam(bamFile)
    y = x[[1]]
    z=do.call("DataFrame", y)
    z$strand=factor(as.factor(z$strand))
    
    z = as.data.frame(z)
    ##remove all this stuff since we dont use it to save memory
    z$mrnm  = NULL
    z$rname = NULL
    z$mapq  = NULL
    z$flag  = NULL	
    z$mpos  = NULL
    z$qual  = NULL
    z$qname = as.character(z$qname)
    {
        samtools = file.path(allProps$resourcesDir, allProps$samtoolsDir, "samtools")
        if(is.null(allProps$samtoolsDir))
        {
            samtools = "samtools"
        }
        samFilename = paste(bamFile,".sam", sep="")
        command = paste(samtools," view ", bamFile, " -h -o ", samFilename,sep="")
        system(command)
        cigarStrings = read.table(samFilename, skip=5)$V6
        z$cigar = as.character(cigarStrings)
        z$numMM=convertToErrorCounts(z$cigar,"X")
        z$numI=convertToErrorCounts(z$cigar,"I")
        z$numD=convertToErrorCounts(z$cigar,"D")
    }
    
    if(allProps$turnOffCigar)
    {
        z$cigar = NULL
    }
    if(allProps$turnOffSeq)
    {
        z$seq   = NULL
        
    } else
    {
        z$seq = as.character(z$seq)
    }
    
    z$ID = rep(-1,nrow(z))
    z$readEnd = rep(NA, nrow(z))
    z$alignCounter = rep(-1,nrow(z))
    readSides = c("l","r")
    for (readSide in readSides)
    {
        splittedNames = strsplit(z$qname, readSide)
        ids = sapply(splittedNames, function (x) x[1])
        indexes= ids!=z$qname
        alignCounters = splittedNames[indexes]
        alignCounters = sapply(alignCounters, function (x) x[2])
        
        z$ID[indexes]=as.integer(ids[indexes])
        z$alignCounter[indexes]=as.integer(alignCounters)
        z$readEnd[indexes]=readSide
    }	
    z$qname = NULL
    z$readEnd = as.factor(z$readEnd)
    
    d1 = z[z$readEnd=="r",]
    d2 = z[z$readEnd=="l",]
    
    z = merge(d2,d1, by=c("ID","alignCounter"), all.x=T, all.y=T)
    z$alignCounter=NULL
    
    trimData = read.table(trimFile, header=T)
    if(nrow(trimData)!=(max(z$ID)-min(z$ID)+1))
    {
        stop(paste(nrow(trimData), max(z$ID)-min(z$ID)+1))
    }
    trimData$ID = seq(min(z$ID), max(z$ID))
    z = merge(z, trimData, by="ID")

    colnames(z) = sub(pattern=".x", replacement = ".3", x=colnames(z), fixed=T)
    colnames(z) = sub(pattern=".y", replacement = ".1", x=colnames(z), fixed=T)

    z$isize = z$isize.3 #using the 3 end to ensure we get NA's rather than 0 when we have an unpaired set
    z$isize.1 = NULL
    z$isize.3 = NULL
    z$readEnd.1=NULL
    z$readEnd.3=NULL
    z$strand = z$strand.1
    z$strand.1=NULL
    z$strand.3=NULL
    z$pos.3=NULL #THis is already taken into account by isize
    
    z$ID = z$ID + startId
    z$experimentID = experimentID
    z$jackpotIndex = allJackpotTags[z$jackpotIndex]
    colnoindel = (is.na(z$numD.1)|(z$numD.1+z$numI.1)==0) & (is.na(z$numD.3)|(z$numD.3+z$numI.3)==0)
    noindelIDs = z[colnoindel,]$ID
    remainingIDs =  setdiff(z$ID, noindelIDs)
    z = data.table(z)
    setkey(z, ID)
    neededIndel = z[J(remainingIDs),]
    z = rbind(z[colnoindel,],neededIndel)
    z$strand.1=NULL
    return(z)	
    
}

getName <- function(experimentRecord)
{
    return(experimentRecord$experimentID)
}

getBlocksForExperiment <- function(experimentID, experimentData)
{
    mask=experimentData$experimentID==experimentID
    experimentRecord = experimentData[mask,]
    read1file = getReadFile(experimentID, experimentRecord$run, readIndex=1)
    if (is.na(read1file))
    { 
        return(c(0))
    }
    if (allProps$maxNumReads!=-1) 
    {
        numReads = allProps$maxNumReads
    }
    else
    {
        numReadsCommand = paste("wc -l ", file.path(read1file), sep="")
        
        print(numReadsCommand)
        numReads = as.integer(strsplit(system(numReadsCommand, intern=T), " ")[[1]][1])/2
    }
    
    
    print(paste("num reads for ", getName(experimentRecord), ":" , numReads), sep="")
    
    
    defaultBlockSize = allProps$readProcessingBlockSize
    blocks = seq(1, numReads, defaultBlockSize)
    blocks = append(blocks, numReads + 1)
    return(blocks)
}


getMergedReadFrameFromBlocks <- function(allProps)
{
    print(paste("started loading:", Sys.time()))
    i = 1
    allReadFrames = list()
    experimentData = getExperimentData()
    
    totalRows = 0

    validExperimentalIds=getValidExperimentalIds( experimentData = experimentData)
    for (experimentID in validExperimentalIds)
    {

        print(experimentID)
        blocks = try(getBlocksForExperiment(experimentID, experimentData))
        
        if(class(blocks)!="try-error" & blocks[1]!=0)
        {
            for(b in 1:(length(blocks)-1))
            { 
                startId  = blocks[b]
                blockLen = format(blocks[b+1]-blocks[b], scientific=F)
                afolder = file.path(allProps$outputDir, "blocks", paste(experimentID,"_",startId, "_",blockLen,sep=""))
                bamFile  = file.path(afolder, "javaNew.bam")
                trimFile = file.path(afolder, "jackpotIndex.txt")
                readFrame = 
                    tryCatch(getReadFrameForPooledFile(bamFile, trimFile, experimentID, startId),	
                             error=function(e) {print(paste("Failure!!!", bamFile, trimFile, experimentID, startId, e)); return(NULL)}); #keep going if not enough data
                print(paste("loading:",bamFile))
                print(paste("frame size:",object.size(readFrame)))
                if(is.null(readFrame))
                {
                    next;
                }
                allReadFrames[[i]] = readFrame
                totalRows = totalRows + nrow(readFrame)
                i = i+1
            }
            print(totalRows)
            
            gc()
        }
    }
    print(paste("started rbind:", Sys.time()))
    print(paste("mem used on readFrames list",object.size(allReadFrames)))

##    print(str(allReadFrames[[1]]))
    print(lsos())
    save(file=file.path(allProps$outputDir, allProps$postProcessDir, "savedMergedList"), list=ls())
    gc()
    newReadFrame=do.call(rbind, allReadFrames)
    ## print(str(newReadFrame))
    ## print("mem used on newReadFrame")
    ## print(object.size(newReadFrame))
    
    gc()
    rm(allReadFrames)
    gc()
    
    print(paste("finished loading into frame:", Sys.time()))
    print(paste("nrow after binding", nrow(newReadFrame)))
    
    return(newReadFrame)
}

getFreqNormalizingTable <- function(readFrame) 
{
    normalizing=data.table(data.frame(table(readFrame$experimentID)))
    normalizing$experimentID =normalizing$Var1
    normalizing$Var1=NULL
    setkey(normalizing, "experimentID")
    return(normalizing)
}

filterReadFrame <- function(groupCols, readFrame, removedRows, measuredValue, removeRows=T) 
{
    for(groupCol in groupCols)
    {
        atab=table(readFrame[[groupCol]], removedRows,  dnn=c(groupCol, measuredValue))
        print(atab)
        printPercentage(atab)
        ratios = (atab[,1]/(atab[,1]+atab[,2]))
        minval = min(ratios)
        try(print(t(t(ratios/minval))))
    }
    x=gc()
    print(measuredValue)
    if(removeRows)
    {
        readFrame = readFrame[which(removedRows),]
        print(paste("rows remaining after removal:",nrow(readFrame)))
    }
    x=gc()
    return(readFrame)
}

printPercentage <- function(tabl)
{
    nrml = (format(tabl/rowSums(tabl), scientific=F))
    print(nrml)
    return(nrml)
}

getProcessedReadFrame <- function(tbref, readFrame=NULL)
{
    experimentData = getExperimentData()
    experimentData = experimentData[,j=list(experimentID,fullname)]
    setkey(experimentData, "experimentID")
    
    print("getting processed readframe")
    if(is.null(readFrame))
    {
        print("loading up readframe from bam file")
        readFrame = getMergedReadFrameFromBlocks(allProps = allProps)
        
        x=gc()
        
        save(file=file.path(allProps$outputDir, allProps$postProcessDir, "savedMerged"), list=ls())
    }
    if (!(TRUE %in%(class(readFrame)=="data.table")))
    {
        print("converting to data.table")
        readFrame=data.table(readFrame)
    }

    if(allProps$removeQualityMetrics)
    {
        readFrame$numD.1 = NULL
        readFrame$numMM.1 = NULL
        readFrame$numI.1 = NULL
        readFrame$numD.3 = NULL
        readFrame$numMM.3 = NULL
        readFrame$numI.3 = NULL
        readFrame$a_insdel.1 = NULL
        readFrame$a_mm.1 = NULL
        readFrame$a_insdel.3 = NULL
        readFrame$a_mm.3 = NULL
    }
    
    groupCols = "experimentID"
    
    removedRows <- !is.na(readFrame$strand)
    measuredValue = "read1_aligned"
    readFrame = filterReadFrame(groupCols = groupCols, readFrame = readFrame, removedRows = removedRows, measuredValue = measuredValue)

    removedRows <- !is.na(readFrame$qwidth.3)
    measuredValue = "read_3_aligned"
    readFrame = filterReadFrame(groupCols = groupCols, readFrame = readFrame, removedRows = removedRows, measuredValue = measuredValue)

    readFrame$qwidth.3=NULL
    
    x=gc()
    
    revStrand = which(readFrame$strand=="-")
    posStrand = which(readFrame$strand=="+")
    revStrandPos = readFrame$pos.1[revStrand]
    forStrandPos = readFrame$pos.1[posStrand]
    x=gc()
    
    orientation = "FF"#getPairOrientation(allProps)
    if(orientation =="FR")
    {
        seqlen     = readFrame$qwidth.1[revStrand]
        fuspos     = (revStrandPos + seqlen - 1)

        readFrame$endPos.1=rep(-1, length(readFrame$pos.1))
        readFrame$endPos.1[revStrand] = fuspos
        
        seqlen     = readFrame$qwidth.1[posStrand]
        fuspos     = forStrandPos
        readFrame$endPos.1[posStrand] = getRevComplementPos(tbref$length, pos = fuspos)
        readFrame$pos.1[posStrand]    = readFrame$endPos.1[posStrand] - seqlen + 1
        
        print("flipping read1 alignment")
        readFrame$strand[revStrand] = "+"
        readFrame$strand[posStrand] = "-"
    }
    if(orientation=="FF")
    {
        revStrandPos = readFrame$pos.1[revStrand]
        seqlen     = readFrame$qwidth.1[revStrand]
        fuspos     = (revStrandPos + seqlen - 1)
        fuspos     = getRevComplementPos(tbref$length, pos = fuspos)
        readFrame$pos.1[revStrand] = fuspos
        readFrame$endPos.1 = readFrame$pos.1 + readFrame$qwidth.1-1
    }

    readFrame$pos.1=NULL

    x=gc()
    rm(revStrand)
    rm(posStrand)
    rm(revStrandPos)		
    rm(forStrandPos)
    rm(seqlen)
    rm(fuspos)
    
    readFrame$experimentID = as.factor(readFrame$experimentID)
    setkey(readFrame, "experimentID")
    readFrame = experimentData[readFrame]
    x=gc()
    
    print(table(readFrame$fullname))
    print(nrow(readFrame))

    readFrame$qwidth.1=NULL
    setkey(readFrame, "strand", "endPos.1")
    readFrame$cutDist = tbref$positionData[readFrame]$cutDist
    
    print(".95,.99,.995 quantiles for distance from nearest known cut site")
    print(quantile(abs(readFrame$cutDist),c(.95,.99,.995)))
    
    fixDists=allProps$fixDists
    
    setkey(readFrame, "strand", "endPos.1")
    removedRows = abs(readFrame$cutDist)<5
    measuredValue = "lt_5bp_to_cut_site"
    readFrame = filterReadFrame(groupCols = groupCols, readFrame = readFrame, removedRows = removedRows, measuredValue = measuredValue,removeRows=fixDists)
    x=gc()
    
    print("neg cut dist bias?")
    atab=(table(readFrame$experimentID,
                cut(readFrame$cutDist, 
                    c(-999,-.001,.001,999),
                    labels=c("neg","zero","pos")), dnn=c("experiment","neg cut dist")))
    print(atab)
    printPercentage(atab)
    atab=atab[,-2]
    try(print(riskratio(atab,rev="columns")))
    x=gc()
    
    if(fixDists)
    {
        readFrame$endPos.1=readFrame$endPos.1+readFrame$cutDist
        readFrame$cutDist=NULL
    }
    
    setkey(readFrame, "strand","endPos.1")
    readFrame$geneIndex =  tbref$positionData[readFrame]$geneIndex
    readFrame$geneStatus = tbref$positionData[readFrame]$geneStatus
    gc()
    
    print("gene status table")
    atab = table(readFrame$fullname, readFrame$geneStatus)
    print(atab)
    printPercentage(atab)
    
    removedRows <- !is.na(readFrame$jackpotIndex)
    measuredValue = "parsed_jackpot_index"
    readFrame = filterReadFrame(groupCols = groupCols, readFrame = readFrame, removedRows = removedRows, measuredValue = measuredValue,removeRows=allProps$removeReadsWithoutJackpot)
    
    readFrame$geneStatus=factor(readFrame$geneStatus, levels =c("inFrame","outOfFrame","downstreamOfGene","upstreamOfGene"))
    
    readFrame$experimentID=factor(readFrame$experimentID)
    x=gc()
    
    if(!allProps$removeQualityMetrics)
    {
	print("a_mm.1")
	atab = table(readFrame$fullname, (readFrame$a_mm.1))
	printPercentage(atab)
	gc()
	
	print("a_mm.3")
	atab = table(readFrame$fullname, cut(readFrame$a_mm.3, c(-10,0,3,50)))
	printPercentage(atab)
	gc()
	
	print("a_insdel.1")
	atab = table(readFrame$fullname, readFrame$a_insdel.1)
	printPercentage(atab)
	gc()
	
	print("a_insdel.3")
	atab = table(readFrame$fullname, cut(readFrame$a_insdel.3, c(-10,0,1,2,3,50)))
	printPercentage(atab)
	gc()
	
	print("numD.1")
	atab = table(readFrame$fullname, readFrame$numD.1)
	printPercentage(atab)
	gc()
	
	print("numD.3")
	atab = table(readFrame$fullname, readFrame$numD.3)
	printPercentage(atab)
	gc()
	
	print("numI.1")
	atab = table(readFrame$fullname, readFrame$numI.1)
	printPercentage(atab)
	gc()
	
	print("numI.3")
	atab = table(readFrame$fullname, readFrame$numI.3)
	printPercentage(atab)
	gc()
	
	print("numMM.1")
	atab = table(readFrame$fullname, readFrame$numMM.1)
	printPercentage(atab)
	gc()
	
	print("numMM.3")
	atab = table(readFrame$fullname, readFrame$numMM.3)
	printPercentage(atab)
	gc()
    }
    
    return(readFrame)
}

getUniqueJackpotInfo <- function(jackpotIndexes,cutDists)
{
    uniqueIndexes =unique(jackpotIndexes)
    numFusions = length(uniqueIndexes)
    numRemovedFusions = length(jackpotIndexes)-numFusions
    maxSingleJackpotRemoved = max(table(jackpotIndexes))-1
    return(list(numFusions = numFusions, 
                maxSingleJackpotRemoved = maxSingleJackpotRemoved, 
                numRemovedFusions = numRemovedFusions, 
                cutDist=cutDists[1]))
}

getEventsPerGenePerDist <- function(eventsPerPosPerJackpot, normalizing, uniqueGenes) 
{
    eventsPerGenePerDist=
        eventsPerPosPerJackpot[,j=list(numFusionsWithJackpot = sum(numFusions), 
                                       numWholeFusions = sum(numWholeFusions>0),
                                       numPartialFusions = sum(numPartialFusions),
                                       numFusions=sum(pmin(numFusions,1))), 
                               by=c("experimentID","strand","endPos.1")]
    print("first")
    setkey(eventsPerGenePerDist,experimentID)
    
    eventsPerGenePerDist[ ,freqFusions:=eventsPerGenePerDist[normalizing, numFusions/Freq]]
    setkey(eventsPerGenePerDist, "strand","endPos.1")
    
    eventsPerGenePerDist$geneIndex = tbref$positionData[eventsPerGenePerDist,  j=list(geneIndex=geneIndex)]$geneIndex
    eventsPerGenePerDist$geneStatus = tbref$positionData[eventsPerGenePerDist, j=list(geneStatus=geneStatus)]$geneStatus
    
    setkey(eventsPerGenePerDist, geneIndex)
    setkey(uniqueGenes, geneIndex)
    eventsPerGenePerDist = uniqueGenes[eventsPerGenePerDist]
    
    x=gc()
    setkey(eventsPerGenePerDist,experimentID)
    expData = getExperimentData()
    expData$read1file = NULL
    expData$read3file = NULL
    eventsPerGenePerDist = expData[eventsPerGenePerDist]
    setkey(eventsPerGenePerDist,experimentID,strand, endPos.1)
    
    return(eventsPerGenePerDist)
}

getEventsPerGenePerJackpot <- function(readFrame, normalizing, uniqueGenes) 
{
    eventsPerGenePerDistPerJackpot=readFrame[,j=list(numFusions = sum(readFraction),
                                                     numPartialFusions = sum(readFraction<1&readFraction>0), 
                                                     numWholeFusions = sum(readFraction==1)),
                                             by=c("experimentID","geneIndex","jackpotIndex","geneStatus")]
    setkey(eventsPerGenePerDistPerJackpot, experimentID)
    
    eventsPerGenePerDistPerJackpot[ ,freqFusions:=eventsPerGenePerDistPerJackpot[normalizing, numFusions/Freq]]
    setkey(eventsPerGenePerDistPerJackpot, geneIndex)
    setkey(uniqueGenes, geneIndex)
    eventsPerGenePerDistPerJackpot = uniqueGenes[eventsPerGenePerDistPerJackpot]
    
    x=gc()
    setkey(eventsPerGenePerDistPerJackpot,experimentID,geneIndex, jackpotIndex,geneStatus)
    
    expData = getExperimentData()
    expData$read1file = NULL
    expData$read3file = NULL
    eventsPerGenePerDistPerJackpot = expData[eventsPerGenePerDistPerJackpot]
    setkey(eventsPerGenePerDistPerJackpot,experimentID,geneIndex, jackpotIndex,geneStatus)
    
    
    return(eventsPerGenePerDistPerJackpot)
}

getEventsPerGenePerPerDistPerJackpot <- function(readFrame, normalizing, tbref) 
{
    uniqueGenes = tbref$uniqueGenes
    positionData = tbref$positionData
    readFrame$distFromStart = positionData[readFrame, j=distFromGeneStart]$distFromGeneStart
    eventsPerGenePerDistPerJackpot=readFrame[,j=list(numFusions = sum(readFraction), numPartialFusions = sum(readFraction<1), numWholeFusions = sum(readFraction==1)),			
                                             by=c("experimentID","geneIndex","jackpotIndex","distFromStart", "geneStatus")]
    setkey(eventsPerGenePerDistPerJackpot, experimentID)
    
    eventsPerGenePerDistPerJackpot[ ,freqFusions:=eventsPerGenePerDistPerJackpot[normalizing, numFusions/Freq]]
    setkey(eventsPerGenePerDistPerJackpot, geneIndex)
    eventsPerGenePerDistPerJackpot = eventsPerGenePerDistPerJackpot[uniqueGenes]
    eventsPerGenePerDistPerJackpot = eventsPerGenePerDistPerJackpot[!is.na(eventsPerGenePerDistPerJackpot$experimentID),]
    x=gc()

    setkey(eventsPerGenePerDistPerJackpot,experimentID,geneIndex, jackpotIndex,geneStatus)
    expData = getExperimentData()
    expData$read1file = NULL
    expData$read3file = NULL
    eventsPerGenePerDistPerJackpot = expData[eventsPerGenePerDistPerJackpot]
    setkey(eventsPerGenePerDistPerJackpot,experimentID,geneIndex, jackpotIndex,geneStatus)
    
    return(eventsPerGenePerDistPerJackpot)
}

getEventsPerPosPerJackpot <- function(readFrame, normalizing, uniqueGenes, mergeToGeneInfo=F) 
{
    eventsPerGenePerDistPerJackpot=readFrame[,j=list(numFusions = sum(readFraction), numPartialFusions = sum(readFraction<1), numWholeFusions = sum(readFraction==1)),			
                                             by=c("experimentID","strand","endPos.1","jackpotIndex")]
    setkey(eventsPerGenePerDistPerJackpot, experimentID)
        
    eventsPerGenePerDistPerJackpot[ ,freqFusions:=eventsPerGenePerDistPerJackpot[normalizing, numFusions/Freq]]
    setkey(eventsPerGenePerDistPerJackpot, "strand", "endPos.1")
    eventsPerGenePerDistPerJackpot$geneIndex = tbref$positionData[eventsPerGenePerDistPerJackpot, j=list(geneIndex=geneIndex)]$geneIndex
    eventsPerGenePerDistPerJackpot$geneStatus = tbref$positionData[eventsPerGenePerDistPerJackpot, j=list(geneStatus=geneStatus)]$geneStatus
    
    setkey(eventsPerGenePerDistPerJackpot, geneIndex)
    if(mergeToGeneInfo)
    {
        eventsPerGenePerDistPerJackpot = eventsPerGenePerDistPerJackpot[uniqueGenes]
        eventsPerGenePerDistPerJackpot = eventsPerGenePerDistPerJackpot[!is.na(eventsPerGenePerDistPerJackpot$experimentID),]
    }
    x=gc()
    setkey(eventsPerGenePerDistPerJackpot,experimentID,geneIndex, jackpotIndex,geneStatus)
    
    expData = getExperimentData()
    expData$read1file = NULL
    expData$read3file = NULL
    eventsPerGenePerDistPerJackpot = expData[eventsPerGenePerDistPerJackpot]
    setkey(eventsPerGenePerDistPerJackpot,experimentID,geneIndex, jackpotIndex,geneStatus)
    
    return(eventsPerGenePerDistPerJackpot)
}

getEventsPerGene <- function(eventsPerGenePerDist, normalizing, tbref)
{
    uniqueGenes= tbref$uniqueGenes
    eventsPerGenePerDist$distFromStart = tbref$positionData$distFromGeneStart[eventsPerGenePerDist$endPos.1]
    eventsPerGenePerJackpot=eventsPerGenePerDist[,list(
        uniqueFusionPoints=length(numFusions),
        numFusions=sum(numFusions),
        numWholeFusions = sum(numWholeFusions),
        numPartialFusions = sum(numPartialFusions),
        maxFusionsPerUnique=max(numFusions), 
        medFusionsPerUnique=getMedian(numFusions),
        maxDistFromStart = distFromStart[which.max(numFusions)],
        avgWeightedDistFromStart = sum(distFromStart*numFusions)/sum(numFusions)),
        by=c("experimentID","geneIndex","geneStatus")]
    
    setkey(eventsPerGenePerJackpot, experimentID)
    freqFusions = eventsPerGenePerJackpot[normalizing, list(freqFusions=numFusions/Freq)]$freqFusions
    eventsPerGenePerJackpot$freqFusions=freqFusions
    setkey(eventsPerGenePerJackpot, geneIndex)
    eventsPerGenePerJackpot = uniqueGenes[eventsPerGenePerJackpot]
    eventsPerGenePerJackpot=eventsPerGenePerJackpot[!is.na(eventsPerGenePerJackpot$experimentID),]
    x=gc()
    theOrder = order(eventsPerGenePerJackpot$experimentID, eventsPerGenePerJackpot$geneStatus, eventsPerGenePerJackpot$numFusions, decreasing=T)
    eventsPerGenePerJackpot=eventsPerGenePerJackpot[theOrder,]
    setkey(eventsPerGenePerJackpot,experimentID,geneIndex, geneStatus)
    x=gc()
    
    expData = getExperimentData()
    expData$read1file = NULL
    expData$read3file = NULL
    eventsPerGenePerJackpot = expData[eventsPerGenePerJackpot]
    setkey(eventsPerGenePerJackpot,experimentID,geneIndex, geneStatus)
    
    return(eventsPerGenePerJackpot)
}

buildMiriamEventsByColumn <- function(eventsPerGene, uniqueGenes, normalizing)
{
    allGenesRv = data.table(data.frame(locus = tbref$uniqueGenes$locus))
    allGenesRvPlus = data.table(data.frame(locus = tbref$uniqueGenes$locus, symbol=tbref$uniqueGenes$symbol, name=tbref$uniqueGenes$name))
    setkey(allGenesRvPlus, "locus","symbol","name")
    setkey(allGenesRv, "locus")
    setkey(eventsPerGene, experimentID,geneStatus,locus)
    
    frameStatuses = c("inFrame", "outOfFrame")
    experimentData = getExperimentData()
    
    maxReads = max(normalizing$Freq)
    inOutFrameData = list()
    for(expID in experimentData$experimentID)
    {
        
        inOutFrameData[[expID]] =list()
        prefix = experimentData[i=expID,]$fullname
        print(prefix)
        for(frameStatus in frameStatuses)
        {
            uniqueFusionsColName = paste(prefix, "UniqueFusions", sep="")
            forMediumAndStatus  =  eventsPerGene[J(experimentID=expID, geneStatus=frameStatus),    
                                                 j=list(uniqueFusionPoints, numFusions, freqFusions=freqFusions*maxReads,
							numWholeFusions,
							numPartialFusions), by="locus"]
            setkey(forMediumAndStatus, "locus")
            forMediumAndStatus = forMediumAndStatus[allGenesRv]
            
            colnamesToRename = c("uniqueFusionPoints","numFusions", "numWholeFusions", "numPartialFusions", "freqFusions")
            for(name in colnamesToRename)
            {
                setnames(forMediumAndStatus, name, paste(prefix,"_",name,sep=""))
            }
            
            inOutFrameData[[expID]][[frameStatus]]= forMediumAndStatus
            
        }
        inOutRatioCol = paste(prefix, "InOutRatio", sep="")
        
        inFrameData                   = inOutFrameData[[expID]][["inFrame"]]
        outFrameData                  = inOutFrameData[[expID]][["outOfFrame"]]
        joinedStatuses                = inFrameData[outFrameData]
        num = joinedStatuses[[paste(prefix, "_freqFusions",sep="")]]
        denom = joinedStatuses[[paste0("i.",prefix, "_freqFusions")]]
        inOutFrameData[[expID]][["inFrame"]][[inOutRatioCol]]  = num/denom
    }
    
    experimentalIds = experimentData$experimentID
    firstID = experimentalIds[1]
    
    joined = inOutFrameData[[firstID]]$inFrame[allGenesRv]
    for (i in 2:length(experimentalIds))
    {
        expID = experimentalIds[i]
        joined = inOutFrameData[[expID]]$inFrame[joined]
    }
    joined = allGenesRvPlus[joined]
    gc()
    return(joined)
}

getJackpotPlots <- function(readFrame)
{
    amat = table(readFrame$jackpotIndex)
    amat = sort(amat)
    
    x=unlist(lapply(gregexpr("G", names(amat)),length))
    
    pdf(file.path(allProps$outputDir, allProps$postProcessDir, "numGs.pdf"))
    plot(x,amat, ylab="numReads", xlab="number of G in jackpot")
    dev.off()
    
    apwm= PWM(na.omit(as.character(readFrame$jackpotIndex)))
    seqLogo(t(t(apwm) * 1/colSums(apwm)))

    return(amat)
}

convertToErrorCounts <- function(cigarCol,errType)
{
    numMax = 2
    numMM = rep(0,length(cigarCol))
    for (i in 1:numMax)
    {
        qstring= paste(i,errType, sep="")
        a = gregexpr(qstring,cigarCol)
        mask = (grepl(qstring,cigarCol))
        numMM= numMM+unlist(lapply(a,length))*mask*i #weight by the number i of mismatches in the mismatch type
        
    }
    return(as.integer(numMM))
}

frameForOutcome <- function(outcome, eventsPerGenePerDist, strandVal) 
{
    
    formulaString = as.formula(paste(outcome, " ~ ", "endPos.1 + fullname", sep=""))
    fullMat = (xtabs(formulaString, eventsPerGenePerDist[eventsPerGenePerDist$strand==strandVal]))
    
    fullFrame = data.table(data.frame(fullMat[,]))
    genes = data.table(endPos.1 = as.numeric(rownames(fullMat)))
    rm(fullMat)
    fullFrame= cbind(genes, fullFrame)
    oldnames = base::setdiff(colnames(fullFrame), c("strand","endPos.1"))
    setnames(fullFrame, old = oldnames, new = paste(oldnames, "_",outcome, sep=""))

    setkey(fullFrame, "endPos.1")
    return(fullFrame)
}

formFullFrame <- function(eventsPerGenePerDist, tbref, normalizing=NULL) 
{
    thestrands = c("+", "-")
    crossFrames = list()
    counter = 1
    for(thestrand in thestrands)
    {
        outcome = "numFusions"
        fullFrame = frameForOutcome(outcome = outcome, eventsPerGenePerDist = eventsPerGenePerDist,thestrand)
        infoPerPos = tbref$positionData[J(strand = rep(thestrand, nrow(fullFrame)), position = fullFrame$endPos.1), j= list(strand, geneIndex, geneStatus, distFromGeneStart, position)]
        infoPerPos = infoPerPos[infoPerPos$strand==thestrand]
        infoPerPos = unique(infoPerPos)
        setkey(infoPerPos, "position")

        infoPerPos=infoPerPos[fullFrame, j=list(locus = geneIndex, name=geneIndex, symbol = geneIndex, geneIndex = geneIndex, geneStatus, distFromGeneStart)]
        infoPerPos$locus = tbref$uniqueGenes$locus[as.integer(infoPerPos$locus)]
        infoPerPos$name  = tbref$uniqueGenes$name[as.integer(infoPerPos$name)]
        infoPerPos$symbol  = tbref$uniqueGenes$symbol[as.integer(infoPerPos$symbol)]


        strands = infoPerPos$strand

        fullFrame = cbind(infoPerPos, fullFrame)
        setkey(fullFrame, "endPos.1", "fullname_numFusions")
        
        outcome = "numWholeFusions"
        partFrame = frameForOutcome(outcome = outcome, eventsPerGenePerDist = eventsPerGenePerDist,thestrand)
        setkey(partFrame, "endPos.1", "fullname_numWholeFusions")
        
        fullFrame = fullFrame[partFrame] 
        
        outcome = "numPartialFusions"
        partFrame = frameForOutcome(outcome = outcome, eventsPerGenePerDist = eventsPerGenePerDist,thestrand)
        setkey(partFrame, "endPos.1", "fullname_numPartialFusions")
        fullFrame = fullFrame[partFrame] 
        
        outcome = "numFusionsWithJackpot"
        partFrame = frameForOutcome(outcome = outcome, eventsPerGenePerDist = eventsPerGenePerDist,thestrand)
        setkey(partFrame, "endPos.1", "fullname_numFusionsWithJackpot")
        fullFrame = fullFrame[partFrame]
        
        fullFrame= cbind(data.table(strand = rep(thestrand, nrow(fullFrame))), fullFrame)
        
        crossFrames[[counter]]= fullFrame
        counter = counter +1
    }
    fullFrame = rbindlist(crossFrames)
    fullFrame$endPos.1=NULL
    setkey(fullFrame, geneIndex)
    fullFrame =fullFrame[tbref$geneMetadata]
    fullFrame$geneIndex = NULL
    
    if(!is.null(normalizing))
    {
        setkey(normalizing, "fullname")
        fullFrame = data.frame(fullFrame)
        maxval = max(normalizing$Freq)
        for(aname in normalizing$fullname)
        {
            print(aname)
            thecount = normalizing[J(fullname=aname)]$Freq
            n_const = maxval/thecount
            
            print(n_const)
            thecols = colnames(fullFrame)
            validCols = grepl(aname, thecols)
            inds = which(validCols)
            fullFrame[,inds] = fullFrame[,inds]*n_const
        }
        setkey(normalizing, "experimentID")
    }
    return(fullFrame)
}
