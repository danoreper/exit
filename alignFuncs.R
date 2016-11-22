## Functions for performing alignments
## 
## Author: doreper
###############################################################################

bwa_samtools2 <- function() 
{
    outputs = list()
    bwa= allProps$bwaCom
    samtools = allProps$samtoolsCom
    
    fastaDir = allProps$outputDir
    
    errFile           = file.path(fastaDir, 'bwaIndexOut.txt')
    tbReferenceFile = file.path(allProps$dataDir, allProps$tbReferenceGenomeFile)
    acommand = (paste(bwa, 'index', tbReferenceFile, paste('2>', errFile, sep="")))

    print(getwd())
    print(acommand)
    x = system(acommand, intern=TRUE)
    print("got index!")
    
    alignPrefix = "trimmedRead1.fasta"
    fastaTargetsFile1 = file.path(fastaDir, alignPrefix)
    bwaOptsString = getAlgOptsString("bwa","aln")
    
    fastaAlnFile1      = file.path(fastaDir, paste(alignPrefix, '.sai', sep=""))
    errFile           = file.path(fastaDir, paste(alignPrefix, 'bwaAlnOut.txt',sep=""))
    acommand = paste(bwa, 'aln ', bwaOptsString, tbReferenceFile, fastaTargetsFile1, '>', fastaAlnFile1,
                     paste('2>', errFile, sep=""), sep=" ")
    print(acommand)
    x = system(acommand, intern=TRUE)
    
    
    alignPrefix = "trimmedRead3.fasta"
    fastaTargetsFile3 = file.path(fastaDir, alignPrefix)
    
    print(Sys.time());
    bwaOptsString = getAlgOptsString("bwa","samse")
    alignPrefix = "merged"
    errFile         = file.path(fastaDir, paste(alignPrefix, 'samseOut.txt',sep=""))
    outputs$samFile = file.path(fastaDir, paste('bwa.sam', sep=""))
    acommand = paste(bwa, 'samse', bwaOptsString, tbReferenceFile, fastaAlnFile1, fastaTargetsFile1, '>', outputs$samFile, paste('2>', errFile, sep=""), sep=" ")
    x = system(acommand, intern=TRUE)
    print(acommand)
    return(outputs)
}

jalign <- function(jversion) 
{
    print("calling jalign")
    java = "java -jar";
    outputDir = allProps$outputDir
    
    errFile           = file.path(outputDir, paste(jversion,'.txt', sep=""))
    acommand = paste("nohup", java, paste("../java/lib/",jversion,".jar",sep=""), paste(theargs, collapse = " "), paste(">",errFile, sep=""))
    print(acommand)
    x = system(acommand, intern=TRUE)
}

parsePositions <- function(samFile, numbp) 
{
    print("calling parsePositions")
    java = "java -jar";
    outputDir = allProps$outputDir
    errFile           = file.path(outputDir, 'parse.txt')
    acommand = paste("nohup", java, "../java/lib/ParsePositionsSAMSE.jar",
                     paste(samFile,
                           paste(samFile,".parsed",sep=""),
                           numbp, sep=" "),
                     paste(">",errFile, sep=""))
    print(acommand)
    x = system(acommand, intern=TRUE)
}


buildJackpotCostMatrix <- function() 
{
    jackpotMat = diag(6)
    rnames = c("A", "C", "G", "T", "-", "N")
    ncol = which(rnames=="N") 
    dashcol = which(rnames=="-")
    rownames(jackpotMat) = rnames
    colnames(jackpotMat) = rnames
    jackpotMat[ncol,(1:length(rnames))[-dashcol]]=.25
    jackpotMat[(1:length(rnames))[-dashcol],ncol]=.25
    return(jackpotMat)
}

removeAdaptersRead3 <- function(areads)
{
    print(paste("started remove adapters3:", Sys.time()))
    expectedAdapter   = DNAString("ATGATGGCCGGTGGATTTGTGNNANNANNNTGGTCGTGGTAT")
    bestAdapterRegion = DNAString("ATGATGGCCGGTGGATTTGTGAAAAAAAAATGGTCGTGGTAT")
    ##substitution matrix for the step identifying the whole of the adapter region
    adapterMat = nucleotideSubstitutionMatrix()
     
    maxNumIndel = 3;
    ncharToMatch = nchar(toString(expectedAdapter))+maxNumIndel
    numReads = length(areads)
    expectedAdapterRegions = subseq(areads, 1, ncharToMatch)
    readLength = length(areads[[1]])
    
    ##substitution matrix for the step identifying the jackpot tag region
    jackpotMat = buildJackpotCostMatrix()
    numFailedReads = 0
    jackpotPattern=DNAString("NNANNANNN")
    
    nindex = as.vector(gregexpr("N", jackpotPattern)[[1]])
    
    print(paste("started align adapters3:", Sys.time()))
    adapterAlignment = pairwiseAlignment(
        pattern=expectedAdapterRegions,
        subject=expectedAdapter, 
        substitutionMatrix = adapterMat,
        type='local', 
        gapOpening=0,
        gapExtension=-1)
    
    alignedAdapterString=as.character(subject(adapterAlignment));	
    
    insdel = (deletion(nindel(adapterAlignment)))[,"WidthSum"]+(insertion(nindel(adapterAlignment)))[,"WidthSum"]
    print(paste("started count mismatches 1:", Sys.time()))
    adapterAlignmentPattern <- pattern(adapterAlignment)
    
    nooccurences = length(nindex)
    
    mm = nmismatch(adapterAlignment)- nooccurences #noccurences$V1 # at every N in the pattern there will be a mismatch.
    rm(adapterAlignment)
    gc()
        
    nonAdapterStartIndex = end(adapterAlignmentPattern)+1
    trimmedReads = subseq(areads, nonAdapterStartIndex, readLength)
    rm(nonAdapterStartIndex)
    gc()
    
    print("start aligned adpter 2b")
    alignment2 = pairwiseAlignment(pattern=alignedAdapterString, 
                                   subject=jackpotPattern,
                                   substitutionMatrix = jackpotMat,
                                   type='local-global', 
                                   gapOpening=0,
                                   gapExtension=-1)
    
    ##must have exact match for jackpot index to be valid
    invalidJackpots = as.character(pattern(alignment2))!=toString(jackpotPattern)
    
    offsets = start(pattern(alignment2))
    print(paste("started loop adapters3:", Sys.time()))
    rm(alignment2)
    
    jackpotRegion = substring(adapterAlignmentPattern, offsets,offsets+length(jackpotPattern)-1)
    
    rm(offsets)
    gc()
    jackpotChars = do.call(rbind, strsplit(jackpotRegion,""))

    rm(jackpotRegion)
    gc()
    jackpotChars = jackpotChars[,nindex]

    callList = list()
    for(index in 1:ncol(jackpotChars))
    {
        callList[[index]]=jackpotChars[,index]
    }
    rm(jackpotChars)
    gc()
    callList$sep=""
    jackpotStrings = do.call(paste, callList)
    jackpotStrings[invalidJackpots]=NA
    
    print(paste("end loop adapters3:", Sys.time()))
    print(paste("ended remove adapters3:", Sys.time()))
    
    return(list(trimmedReads=trimmedReads, jackpotStrings=jackpotStrings, a_mm=mm, a_insdel=insdel))
}


removeAdaptersRead1 <- function(areads)
{
    expectedAdapter   = DNAString("ACCAGCGTTTCTGGATCG")
    bestAdapterRegion = DNAString("ACCAGCGTTTCTGGATCG")
    
    expectedAdapter   = DNAString("GTTTCTGGATCG")
    bestAdapterRegion = DNAString("GTTTCTGGATCG")
    
    adapterMat = nucleotideSubstitutionMatrix()
    maxScore = score(pairwiseAlignment(subject=bestAdapterRegion, 
                                       pattern=expectedAdapter,
                                       substitutionMatrix = adapterMat,
                                       type='local-global', 
                                       gapOpening=0,
                                       gapExtension=-1))
    
    maxNumIndel = 26;
    ncharToMatch = nchar(toString(expectedAdapter))+maxNumIndel
    numReads = length(areads)
    expectedAdapterRegions = subseq(areads, 1, ncharToMatch)
    readLength = length(areads[[1]])
    
    ##substitution matrix for the step identifying th
    adapterAlignments = pairwiseAlignment(pattern=expectedAdapterRegions, 
                                          subject=expectedAdapter,
                                          substitutionMatrix = adapterMat,
                                          type='local-global', 
                                          gapOpening=0,
                                          gapExtension=-1)

    insdel = (deletion(nindel(adapterAlignments)))[,"WidthSum"]+(insertion(nindel(adapterAlignments)))[,"WidthSum"]
    adapterAlignmentPattern <- pattern(adapterAlignments)
     
    subjectNs = gregexpr("N", expectedAdapter)[[1]]
    subjectNs=length(subjectNs[[1]])*grepl("N", expectedAdapter)
    
    nooccurences = subjectNs #+ patternNs #subject+pattern n's. yes we could double up, but not really a big deal. 
    mm = nmismatch(adapterAlignments)- nooccurences #noccurences$V1 # at every N in the pattern there will be a mismatch.
    
    scoreDifferentials = maxScore - score(adapterAlignments)
    nonAdapterStartIndex = end(adapterAlignmentPattern)+1
    trimmedReads = subseq(areads, nonAdapterStartIndex, readLength)
    
    orientation = getPairOrientation(allProps)
    ##	if(orientation =="FF")
    {
        trimmedReads = reverseComplement(trimmedReads)
    }
    return(list(trimmedReads=trimmedReads, a_mm = mm, a_insdel=insdel, scoreDifferentials = scoreDifferentials))
}

removeAdaptersWrapper <- function(experimentRecord, nrec, startRecord)
{
    
    print(paste("started start remove adapters wrapper:", Sys.time()))
    dir.create(allProps$outputDir, showWarnings = FALSE)
    print(allProps$outputDir)
    print(commandArgs(trailingOnly=TRUE))
    print(nrec)
    print(startRecord)
    
    print("processing reads1")
    fastaFilename=getReadFile(experimentRecord$experimentID, experimentRecord$run, readIndex=1)
    print(fastaFilename)
    areads = readDNAStringSet(fastaFilename, format="fasta",nrec=nrec,skip=startRecord-1)
    
    trimmedReadInfo1 = removeAdaptersRead1(areads)
    gc() 
    
    print("processing reads3")
    fastaFilename=getReadFile(experimentRecord$experimentID, experimentRecord$run, readIndex=3)
    print(fastaFilename)
    areads = readDNAStringSet(fastaFilename, format="fasta",nrec=nrec,skip=startRecord-1)
    trimmedReadInfo3 = removeAdaptersRead3(areads)
    
    trimmedReads1File = file.path(allProps$outputDir, "trimmedRead1.fasta")
    writeXStringSet(trimmedReadInfo1$trimmedReads, trimmedReads1File, append=FALSE, format="FASTA")
    
    ##write everything about the reads, but the reads themselves to a separate data frame	
    trimmedReads3File = file.path(allProps$outputDir, "trimmedRead3.fasta")
    writeXStringSet(trimmedReadInfo3$trimmedReads, trimmedReads3File, append=FALSE, format="FASTA")
    
    jackpotIndex = 	jackpotHash[J(jackpotTag=trimmedReadInfo3$jackpotStrings)]$index
    
    noReads = data.frame(
        a_mm.1            = trimmedReadInfo1$a_mm,
        a_insdel.1        = trimmedReadInfo1$a_insdel, 
        a_mm.3            = trimmedReadInfo3$a_mm,
        a_insdel.3        = trimmedReadInfo3$a_insdel,
        jackpotIndex      = jackpotIndex)
    
    
    write.table(noReads, 
                file.path(allProps$outputDir, paste(allProps$jackpotIndexFile, '.txt', sep="")),
                col.names=T,
                row.names =F,
                quote=F,
                sep="\t")
    gc()
                                        
    print(paste("ended start remove adapters wrapper:", Sys.time()))
}
