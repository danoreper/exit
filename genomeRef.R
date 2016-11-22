## Functions used to build up infromation about the TB genomic reference
## 
## Author: doreper
###############################################################################

getGenomeAnnotationRefseq <- function(dataPerStrand) 
{

    tbrefLength = dataPerStrand$bpLength[1]
    
    genomeSummaryFile = file.path(allProps$dataDir, allProps$genomeSummaryFile)
    
    fullAnnot = data.table(read.csv(file.path(allProps$dataDir, allProps$tbAnnotationsFile), header=T, sep="\t"))
    fullAnnot$RefSeq.Locus.Tag= fullAnnot$Locus.Tag
    fullAnnot[fullAnnot$RefSeq.Locus.Tag=="Rv2098c"]$Gene.Symbol="PE_PGRS36"
    fullAnnot$X=NULL

    genomeAnnotation = fullAnnot[fullAnnot$Feature.Type=="CDS"|fullAnnot$Feature.Type=="pseudogene"|fullAnnot$RefSeq.Locus.Tag=="Rv2098c"]
    handAnnotation = data.table(read.csv(file.path(allProps$dataDir, allProps$handAnnotationsFile), header=T, sep="\t"))
    setkey(handAnnotation,"RefSeq.Locus.Tag")
    setkey(genomeAnnotation,"RefSeq.Locus.Tag")
    genomeAnnotation[J(handAnnotation$RefSeq.Locus.Tag),j="Gene.Symbol"]=handAnnotation$Gene.Symbol
    
    tempstart = getRevComplementPos(tbrefLength, genomeAnnotation[as.character(genomeAnnotation$Strand)=="-",]$End)
    tempstop = getRevComplementPos(tbrefLength, genomeAnnotation[as.character(genomeAnnotation$Strand)=="-",]$Start)
    
    genomeAnnotation[as.character(genomeAnnotation$Strand)==as.character("-"),]$Start = tempstart
    genomeAnnotation[as.character(genomeAnnotation$Strand)==as.character("-"),]$End  = tempstop
    
    uniqueGenes = data.table(unique(genomeAnnotation))
    setkey(uniqueGenes, "RefSeq.Locus.Tag")
    numUniqueGenes = nrow(uniqueGenes)
    uniqueGenes$geneIndex = 1:numUniqueGenes
    
    geneMetadata = data.table(read.csv(file.path(allProps$dataDir, allProps$ellenAnnotationsFile), header=T, sep="\t",  fileEncoding = "UTF-8"))
    setkey(geneMetadata, "RefSeq.Locus.Tag")
    
    geneMetadata$geneIndex = uniqueGenes[geneMetadata]$geneIndex
    setkey(geneMetadata, geneIndex)

##    browser()
    if(allProps$includeMetaData)
    {
        setkey(geneMetadata, "RefSeq.Locus.Tag")
        uniqueGenes = merge(uniqueGenes,geneMetadata,by=c("RefSeq.Locus.Tag","geneIndex"))
        setkey(geneMetadata, geneIndex)
    }
    
    print("pre set names")	
    setnames(x = uniqueGenes,
             old = c("RefSeq.Locus.Tag","Strand", "Start", "End", "Gene.Symbol", "Product"),
             new = c("locus", "strand","start","stop","symbol","name")) 
   print("post set names")
    
    uniqueGenes$Locus.Tag=NULL
    uniqueGenes$Accession=NULL
    uniqueGenes$Genome = NULL
    uniqueGenes$Annotation=NULL
    uniqueGenes$Protein.ID=NULL
    uniqueGenes$AntiCodon=NULL
    uniqueGenes$Pseudo.Gene=NULL


    return(list(uniqueGenes=uniqueGenes, geneMetadata=geneMetadata))
}

getVectorsPerStrand <- function(dataPerStrand, defaultFill=NA)
{
    byStrand = list()
    for(i in 1:nrow(dataPerStrand))
    {
        strand = dataPerStrand$strand[i]
        strand = as.character(strand)
        reflength = dataPerStrand$bpLength[i]
        byStrand[[strand]] = rep(defaultFill, reflength)
    }
    return(byStrand)
}

convertListsToDataTable <- function(byStrand)
{
    alldts = list()
    counter = 1
    for(theStrand in names(byStrand))
    {
        theStrand = as.character(theStrand)
        dataForStrand = byStrand[[theStrand]]
        dt = data.table(dataForStrand = dataForStrand,
                        position = 1:length(dataForStrand),
                        strand   = rep(theStrand, length(dataForStrand)))
        
        alldts[[counter]] = dt
        counter = counter +1
    }
    alldts = rbindlist(alldts)
    alldts$strand = factor(alldts$strand)
    setkey(alldts, "strand","position")
    return(alldts)
}

buildGeneIndexPerPosition <- function(uniqueGenes, dataPerStrand) 
{
    indexPerPosition = getVectorsPerStrand(dataPerStrand = dataPerStrand)
    for(j in (1:nrow(uniqueGenes)))
    {
        theStart = uniqueGenes$start[j]
        theEnd   = uniqueGenes$stop[j]
        theStrand = uniqueGenes$strand[j]
        indexPerPosition[[as.character(theStrand)]][theStart:theEnd] = j
    }

    indexPerPosition2 = convertListsToDataTable(byStrand = indexPerPosition)
    setnames(indexPerPosition2, old="dataForStrand", new="geneIndex")
    return(indexPerPosition2)
}

getNearestUpstreamDownstreamGenes <- function(uniqueGenes, dataPerStrand)
{

    print("started up down")
    ordering = order(uniqueGenes$strand, uniqueGenes$start)
    orderedUnique = uniqueGenes[ordering,]
    
    nearestUpstreamGene   = getVectorsPerStrand(dataPerStrand)
    nearestDownstreamGene = getVectorsPerStrand(dataPerStrand)
    for (strandType in dataPerStrand$strand)
    {
        strandType = as.character(strandType)
        orderedGeneSubset = orderedUnique[orderedUnique$strand==strandType,]
        tbrefLength = dataPerStrand$bpLength[dataPerStrand$strand==strandType]
        
        startGeneRecord = orderedGeneSubset[1,]
        stopNonGene = startGeneRecord$start-1
        if(stopNonGene!=0)
        {
            nearestDownstreamGene[[strandType]][1:stopNonGene]=startGeneRecord$geneIndex
        }
        print(tbrefLength)
        for(i in 1:nrow(orderedGeneSubset))
        {
            upstreamGene= orderedGeneSubset[i,]
            startNonGene = upstreamGene$stop+1
            
            if(i<nrow(orderedGeneSubset))
            {
                downstreamGene = orderedGeneSubset[i+1,]
                stopNonGene = downstreamGene$start-1

                nearestDownstreamGene[[strandType]][startNonGene:stopNonGene]=downstreamGene$geneIndex
            }
            else
            {
                stopNonGene = tbrefLength
            }
            if(stopNonGene>=startNonGene)
            {
                nearestUpstreamGene[[strandType]][startNonGene:stopNonGene]=upstreamGene$geneIndex
            }
        }
    }
    print("done looping")
    nearestUpstreamGene = convertListsToDataTable(nearestUpstreamGene)
    nearestDownstreamGene = convertListsToDataTable(nearestDownstreamGene)
    
    setnames(nearestUpstreamGene, old ="dataForStrand", new = "nearestUpstreamGene")
    setnames(nearestDownstreamGene, old ="dataForStrand", new = "nearestDownstreamGene")
    nearestUpAndDownstreamGenes = nearestUpstreamGene[nearestDownstreamGene]
    
    return(nearestUpAndDownstreamGenes)
}

getGeneStatuses <- function(uniqueGenes, dataPerStrand)
{
    geneStatusLabelsPerModulo = factor(c("inFrame","outOfFrame", "outOfFrame", "downstreamOfGene"))
    geneStatuses = getVectorsPerStrand(dataPerStrand, geneStatusLabelsPerModulo[4])
    distFromGeneStart = getVectorsPerStrand(dataPerStrand)

    for(i in 1:nrow(uniqueGenes))
    {
        geneRecord = uniqueGenes[i,]
        start = geneRecord$start
        stop = geneRecord$stop
        strand = as.character(geneRecord$strand)
        distFromGeneStart[[strand]][start:stop] = 0:(stop-start)
        ##+1 to convert distance to a spacn of nucleiotides, +1 to deal with single extra base inthe primer that shifts things out of frame, 
        ##frameShifts = ((distFromGeneStart[[strand]][start:stop] + 2)%%3)
        ##+1 to index into modulo array
        ##geneStatuses[[strand]][start:stop] = geneStatusLabelsPerModulo[frameShifts+1]   
    }
    
    for(strand in dataPerStrand$strand)
    {
        strand = as.character(strand)
        frameShifts = (distFromGeneStart[[strand]] + 2)%%3
        geneStatuses[[strand]] = geneStatusLabelsPerModulo[frameShifts+1]   
        geneStatuses[[strand]][is.na(frameShifts)]=geneStatusLabelsPerModulo[4]
    }
    
    geneStatuses = convertListsToDataTable(geneStatuses)
    distFromGeneStart = convertListsToDataTable(distFromGeneStart)
    
    setnames(geneStatuses,      old ="dataForStrand", new = "geneStatus")
    setnames(distFromGeneStart, old ="dataForStrand", new = "distFromGeneStart")
    geneStatusesAndDists = geneStatuses[distFromGeneStart]
    
    return(geneStatusesAndDists)
}

buildCutSitesFrame <- function(sstr, sstrRev, recSites, recOffsets, dataPerStrand) 
{
    cutSites = getVectorsPerStrand(dataPerStrand)
    cutSites[["-"]] = getCutSites(sstrRev, recSites, recOffsets)
    cutSites[["+"]] = getCutSites(sstr, recSites, recOffsets)
    
    cutSites = convertListsToDataTable(cutSites)
    setnames(cutSites, old ="dataForStrand", new = "cutSites")
    cutSites$cutSites = as.factor(cutSites$cutSites)
    return(cutSites)
}

getCutSites <- function(sstr, recSites, recOffsets) 
{
    cutNames  = rep(NA, nchar(sstr))
    for (name in names(recSites))
    {
        recSite = recSites[[name]]
        recOffsetForSite = recOffsets[[name]]
        
        cutStartsForRecSite = as.integer(gregexpr(recSite, sstr)[[1]] + recOffsetForSite)-1
        cutNames[cutStartsForRecSite]=name
        
        revRecSite=toString(reverseComplement(DNAString(recSite)))
        if(recSite != revRecSite )
        {
            cutStartsForRecSite = gregexpr(revRecSite, sstr)[[1]] + recOffsetForSite -1
            cutNames[cutStartsForRecSite]=name
        }
    }
    return(cutNames)
}

getDistanceToNearestCutSite <- function(positionData, dataPerStrand)
{
    distanceToNearestCutSite = getVectorsPerStrand(dataPerStrand)
    for(strandType in dataPerStrand$strand)
    {
        strandType = as.character(strandType)
        cutPos=which(!is.na(positionData$cutSites[positionData$strand==strandType]))
        prevCut = 1
        for(i in 1:length(cutPos))
        {
            curCut = cutPos[i]
            if(i==1)
            {
                cutDists = (curCut - 1):0	
            }else if(curCut==prevCut)
            {
                cutDists = 0
            }else if(curCut ==(prevCut+1)) 
            {
                cutDists = c(0,0)
            } else
            {
                midpoint = floor((curCut + prevCut)*.5)
                cutDists = c((0:-(midpoint-prevCut)),((curCut-midpoint-1):0))
            }
            
            distanceToNearestCutSite[[strandType]][prevCut:curCut]=cutDists
            prevCut = curCut
        }
        tbrefLength =dataPerStrand[dataPerStrand$strand==strandType]$bpLength
        distanceToNearestCutSite[[strandType]][prevCut:tbrefLength]=(0:(prevCut-tbrefLength))
        
    }	
    distanceToNearestCutSite= convertListsToDataTable(distanceToNearestCutSite)
    setnames(distanceToNearestCutSite, old="dataForStrand", new="cutDist")
    return(distanceToNearestCutSite)
}

getTbBareBones <- function() 
{
    tbReferenceFile = file.path(allProps$dataDir, allProps$tbReferenceGenomeFile)
    tbref = readDNAStringSet(tbReferenceFile, format="fasta")
    return(tbref)
}

getTbRef <- function() 
{
    tbref = getTbBareBones()
    tbrefString    = toString(tbref)
    tbrefStringRev = toString(reverseComplement((tbref)))
    tbrefLength = nchar(tbrefString)
    dataPerStrand = data.table(strand = c("+","-"), bpLength= c(tbrefLength, tbrefLength))
    gcContent = table(strsplit(tbrefString,split=""))
    uniqueGenesAndMeta = getGenomeAnnotationRefseq(dataPerStrand)
    uniqueGenes = uniqueGenesAndMeta$uniqueGenes
    geneMetadata    = uniqueGenesAndMeta$geneMetadata
    positionData = buildGeneIndexPerPosition(uniqueGenes, dataPerStrand)

    moreData = getNearestUpstreamDownstreamGenes(uniqueGenes = uniqueGenes, dataPerStrand)
    positionData = positionData[moreData]
    
    moreData = getGeneStatuses(uniqueGenes = uniqueGenes, dataPerStrand)
    positionData = positionData[moreData]
    
    positionData$geneIndex[is.na(positionData$geneIndex)]=positionData$nearestUpstreamGene[is.na(positionData$geneIndex)]
    positionData$geneIndex[is.na(positionData$geneIndex)]=positionData$nearestDownstreamGene[is.na(positionData$geneIndex)]

    
    print("getting cut sites")
    recSites = list("AcI1"="CCGC", "HpaII"="CCGG")
    recOffsets = list("AcI1"=1, "HpaII"=1)
    moreData = buildCutSitesFrame(sstr = tbrefString,
                                  sstrRev = tbrefStringRev,
                                  recSites = recSites, 
                                  recOffsets = recOffsets,
                                  dataPerStrand = dataPerStrand)
    positionData = positionData[moreData]
    
    print("started distance to cut sites")
    moreData = getDistanceToNearestCutSite(positionData,dataPerStrand)
    positionData = positionData[moreData]
    
    print("loaded cut sites")
    
    tbref = list(tbref=tbref, 
                 string=tbrefString, 
                 stringRev = tbrefStringRev, 
                 length=tbrefLength, 
                 gcContent=gcContent,
                 geneMetadata    = geneMetadata,
                 uniqueGenes=uniqueGenes,
                 numUniqueGenes = nrow(uniqueGenes),
                 positionData=positionData)
    
    return(tbref)
}


plotCutDiffs <- function(positionData) 
{
    positionData = positionData[!is.na(positionData$cutSites)]
    cutDiffs = c(diff(sort(positionData[positionData$strand == "+",]$position)))

    
    pdf(file.path(allProps$outputDir, allProps$postProcessDir, "cutDistHist.pdf"))
    hist(cutDiffs, xlab = "Distance between restriction sites", main= "Hist of distance between restriction sites")
    dev.off()
}

