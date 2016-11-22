'''
Functions for splitting up read alignment into a bunch of blocks of reads, and then 
repeatedly calling the read aligner (i.e. processIllumina)

This can be called on the cluster or locally, though it will take on the order of 100 hours
locally. To change whether its called cluster or locally, see the clusterRunCommand property
Created on Apr 20, 2012
@author: doreper

'''
import loadConfig;
import os;
import subprocess;
import sys;
import random;
import logging;
import time;
import yaml;
import csv;
import time;

log = logging.getLogger()
def hello():
    print("hello")

def getName(experimentRecord):
    return experimentRecord['experimentID']

def getBlocksForExperiment(experimentRecord):
    allProps = loadConfig.allProps
    if(allProps.includedExperiments[0]==-1 or int(experimentRecord['experimentID']) in allProps.includedExperiments):
        numReads = -1
        try:
            if not ('run' in experimentRecord.keys()) or experimentRecord['run']=="":
                raise("no read file defined")
            read1file = getReadFile(experimentRecord['experimentID'], 1, ".reads")
            numReadsCommand = "wc -l " + read1file
            log.info(numReadsCommand)
            numReads = (subprocess.check_output((numReadsCommand), shell=True))
            numReads = int(numReads.split()[0])/2
        except:
            return [0]
    
        if allProps.maxNumReads!=-1:
            numReads = allProps.maxNumReads
        log.info("num reads for " + getName(experimentRecord) + ":" + str(numReads))
        defaultBlockSize = allProps.readProcessingBlockSize
        blocks = range(1, numReads + 1, defaultBlockSize)
        blocks = blocks + [numReads + 1]
    else:
        return [0]
    return blocks



def getReadFile(experimentID, readIndex, suffix=".reads"):
    if(readIndex==""):
        inner = ""
    else:
        inner = "_" +str(readIndex)
    thename = str(os.path.join(loadConfig.allProps.dataDir, loadConfig.allProps.seqDir, "EXIT_"+ str(experimentID) + inner + suffix))
    return(thename)


def getRecordsReader(filename):
    """
    Get a reader for the tab delimited standard output file
    """
    fle = open(filename, 'rb')
    reader = csv.DictReader(fle,dialect = 'excel-tab')
    return reader, fle;

def main(args):
    if len(args) < 2:
        raise Exception("a config filename should have been provided as the first argument to this script.")
    
    configFullFilenames = args[1:len(args)]
    loadConfig.buildConfigMaps(configFullFilenames)
    allProps = loadConfig.allProps
    logDir = os.path.join(loadConfig.allProps.outputDir, loadConfig.allProps.clusterLogDir)
    buildIfDoesNotExist(loadConfig.allProps.outputDir)
    buildIfDoesNotExist(logDir)
    blocksDir = os.path.join(loadConfig.allProps.outputDir, "blocks")
    buildIfDoesNotExist(blocksDir)
    command = loadConfig.allProps.clusterRunCommand
    
    onCluster = 'bsub' in command
    log.info("scripts are running on cluster: " + str(onCluster))

    experimentRecords, fle = getRecordsReader(os.path.join(loadConfig.allProps.dataDir, 
                                                           loadConfig.allProps.seqDir, 
                                                           loadConfig.allProps.experimentDataFile))
    totalNumBlocks = 0
    for experimentRecord in experimentRecords:
        blocks = getBlocksForExperiment(experimentRecord)
        log.info(getName(experimentRecord))
        if blocks == [0]:
            continue
        log.info(blocks)
        for i in range(0,len(blocks)-1):
            blockStart = blocks[i]
            blockEnd   = blocks[i+1]
            blockSize = blockEnd - blockStart 
            randConfigFilename =  _writeConfigForBlock(configFullFilenames, experimentRecord, blockStart, blockSize)
            allConfigs = configFullFilenames + [randConfigFilename]
           
            command = loadConfig.allProps.clusterRunCommand.replace("--args ", "--args " + " ".join(allConfigs))
            name = getName(experimentRecord)
            rLogName = getBlockName(blockStart, blockSize, name)+".rlog"
            command = command + " " 
            command = command + os.path.join(logDir, rLogName)
                
            if onCluster:
                jobname = rLogName+"_job.txt"
                command = specifyLogAndJobName(command, logDir, jobname)
                
            log.info(command)
            try:
                theout = subprocess.check_output(command, shell=True)
            except Exception, err:
                print Exception, err
                #log.info("subproc msg:" +theout)
            
    print("finished cluster distribute!!!")
    
    ##restore original config state.
    loadConfig.buildConfigMaps(configFullFilenames)
    

def getBlockName(blockStart, blockSize, name):
    return "_".join((name, str(blockStart), str(blockSize)))

def _writeConfigForBlock(configFullFilenames, experimentRecord, blockStart, blockSize):
    """
    Generates a config file for a random run. It is basically just setting a config file with a random seed and a distinct output directory, and writing out that config file.
    """
    
    configWrap = loadConfig.yamlReadConfigs(configFullFilenames)
    outputDir = loadConfig.allProps.outputDir
    buildIfDoesNotExist(outputDir)
    name = getName(experimentRecord)
    
    configWrap['outputDir'] = os.path.join(outputDir, "blocks", getBlockName(blockStart, blockSize, name))
    configWrap['blockStart'] = blockStart
    configWrap['readProcessingBlockSize'] = blockSize
    configWrap['experimentID'] = experimentRecord['experimentID']

    buildIfDoesNotExist(os.path.join(outputDir, "configs"))
    randConfigFilename = os.path.join(loadConfig.allProps.outputDir, "configs", "_".join([name, str(blockStart), str(blockSize),'config.txt']))
    randConfigHandle = open(randConfigFilename, 'wb')
    yaml.dump(configWrap, randConfigHandle)
    randConfigHandle.close()
    return randConfigFilename

def specifyLogAndJobName(command, logDir, jobname):
    return command.replace('bsub', 'bsub -o ' + logDir + jobname + '.txt' + ' -J illumina_' + jobname)
def buildIfDoesNotExist(adir):

    if not os.path.exists(adir):
        os.mkdir(adir)

if __name__ == '__main__':    
    args = sys.argv
    main(args)
    

        

