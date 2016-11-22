## Script that is called by cluster distrubute to remove adapters and align reads 
## This is repeatedly called, once per experiment/block, untill all reads from all experiments are processed.        
##
## Bwa_samtools performs a rough alignment of the two ends of each read, while jalign
## then refines the alignment of the short left end read.
##
## Important results of calling this script is that at the end, it generates javanew.bam,
## a file describing the alignment of all reads in a block of reads within an experiment.
##
## Author: doreper
###############################################################################

print(paste("started loading rfuncs:", Sys.time()))
source("rfunclib.R")
print(paste("ended loading rfuncs:", Sys.time()))
dir.create(file.path(allProps$outputDir, allProps$postProcessDir), showWarnings = FALSE)
                                     
print("processIllumina:")

if(allProps$readProcessingBlockSize!=0)
{
    print(paste("started remove adapters overall:", Sys.time()))
    
    experimentData =getExperimentData()

    experimentData = experimentData[experimentID == as.character(allProps$experimentId)]
    {
        experimentRecord = experimentData[1,]
        
        removeAdaptersWrapper(experimentRecord=experimentRecord,
                              nrec = allProps$readProcessingBlockSize,
                              startRecord = allProps$blockStart)
        print(paste("ended remove adapters overall:", Sys.time()))
        if(allProps$alg=="bwa")
        {
            print(paste("started bwa overall:", Sys.time()))
            outs=bwa_samtools2()
            print(paste("ended bwa overall:", Sys.time()))
            
            print(paste("started post processing:", Sys.time()))
            gc()
            print(outs$samFile)
            parsePositions(outs$samFile, width(getTbBareBones()))
            print(paste("ended postprocessing:", Sys.time()))
            
            print(paste("started jalign5 :", Sys.time()))
            rm(outs)

            gc()
            x = lsos()
            print(x)
            jalign("jalign5")
            print(paste("ended jalign5 :", Sys.time()))
            
        }
    }
}
