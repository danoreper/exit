##
## 
## Author: doreper
##
## After clusterDistribute has finished, call this script to assemble all the generated
## bam files into a dataframe, and then write the data frame out in various forms to file:
## 
## Among the output files:
## output/eventsPerPosPerJackpot.csv
## output/eventsPerGenePerDist.csv
## output/eventsPerGene.csv
## output/miriamEvents.csv
## output/miriamDiscreteEvents.csv


#####################################################
readFrame=NULL
x=gc()

source("rfunclib.R")
experimentData = getExperimentData()

x=gc()
dir.create(file.path(allProps$outputDir, allProps$postProcessDir), showWarnings = FALSE)
if(allProps$useProfiler)
{
	Rprof(file.path(allProps$outputDir, allProps$postProcessDir, "Rprof.out"));
}
allStats = list("statType"=c(), "value"=c())
tbref = getTbRef()

partCuts = data.table(strand=tbref$positionData$strand, position=tbref$positionData$position)
partCuts = partCuts[!is.na(tbref$positionData$cutSites)]
setkey(partCuts, "position")
write.table(partCuts, "cutSites.txt", quote=F, sep="\t", col.names=F)

readFrame = getProcessedReadFrame(tbref,readFrame)
print(paste("old nrow:", nrow(readFrame)))
print("loaded genome")


setkey(readFrame, "experimentID","ID", "endPos.1")

readCounts = data.table(data.frame(table(readFrame$ID, readFrame$experimentID, dnn=c("ID","experimentID"))))
setkey(readCounts, "experimentID","ID")
readCounts$ID = as.integer(as.character(readCounts$ID))
setnames(readCounts, "Freq", "readFraction")
readCounts$readFraction[readCounts$readFraction!=0] = 1.0/readCounts$readFraction[readCounts$readFraction!=0]
setkey(readCounts, "experimentID","ID")
readFrame= readCounts[readFrame]

rm(readCounts)
readFrame$strand.1 ==NULL

plotCutDiffs(tbref$positionData)
amat = getJackpotPlots(readFrame = readFrame)
jackpotDistribution = table(readFrame$geneIndex, readFrame$jackpotIndex)

normalizing = getExperimentData()
normalizing$Freq = rep(0, nrow(normalizing))
for(expname in normalizing$fullname)
{
    mask=normalizing$fullname==expname
    experimentRecord = normalizing[mask,]
    numReadsCommand = paste("wc -l ",  
                            getReadFile(experimentRecord$experimentID, experimentRecord$run, readIndex=1),
                            sep="")
    print(numReadsCommand)
    numReads = try(as.integer(strsplit(system(numReadsCommand, intern=T), " ")[[1]][1])/2)
    if(class(numReads)=="try-error")
    {
        next
    }
    normalizing[mask,]$Freq = numReads
}
setkey(normalizing,"experimentID")


x=gc()
print("got fully processed read frame")
save(file=file.path(allProps$outputDir, allProps$postProcessDir, "fullyProcessReads.RData"), list=c("readFrame"))
gc()

lsos()

eventsPerPosPerJackpot = getEventsPerPosPerJackpot(readFrame,normalizing,tbref$uniqueGenes, mergeToGeneInfo=F)
save(file=file.path(allProps$outputDir, allProps$postProcessDir, "eventsPerPosPerJackpot.RData"), list=c("eventsPerPosPerJackpot"))


eventsPerGenePerJackpot = getEventsPerGenePerJackpot(readFrame, normalizing, tbref$uniqueGenes)
save(file=file.path(allProps$outputDir, allProps$postProcessDir, "eventsPerGenePerJackpot.RData"), list=c("eventsPerGenePerJackpot"))
x=gc()

if(allProps$removeQualityMetrics)
{
	rm(readFrame)
	x = gc()
}

eventsPerGenePerDist = getEventsPerGenePerDist(eventsPerPosPerJackpot=eventsPerPosPerJackpot, normalizing, tbref$uniqueGenes)
save(file=file.path(allProps$outputDir, allProps$postProcessDir, "eventsPerGenePerDist.RData"), list=c("eventsPerGenePerDist"))
x=gc()

eventsPerGene        = getEventsPerGene(eventsPerGenePerDist, normalizing, tbref)
save(file=file.path(allProps$outputDir, allProps$postProcessDir, "eventsPerGene.RData"), list=c("eventsPerGene"))
x=gc()
lsos()

save(file=file.path(allProps$outputDir, allProps$postProcessDir, "eventsPerGeneAll.RData"), list=ls())

readsPerExperiment = eventsPerPosPerJackpot[,list(numReads=sum(numFusions)),by="experimentID,geneStatus"]
print("following dejackpotting")
print(readsPerExperiment)
print("total num reads in frame:")
print(sum(readsPerExperiment[readsPerExperiment$geneStatus=="inFrame"]$numReads))

miriamEventsByColumn = buildMiriamEventsByColumn(eventsPerGene, tbref$uniqueGenes, normalizing)

	
setkey(eventsPerGenePerDist, "strand", "endPos.1")

fullFrame = formFullFrame(eventsPerGenePerDist = eventsPerGenePerDist, tbref = tbref, normalizing)


writePostProcessedCSV(eventsPerPosPerJackpot, "eventsPerPosPerJackpot.csv")
writePostProcessedCSV(eventsPerGenePerDist, "eventsPerGenePerDist.csv")
writePostProcessedCSV(eventsPerGene, "eventsPerGene.csv")
writePostProcessedCSV(miriamEventsByColumn,"miriamEvents.csv")
writePostProcessedCSV(fullFrame,"miriamDiscreteEvents.csv")
