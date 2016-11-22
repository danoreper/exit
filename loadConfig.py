'''
Created on Mar 19, 2012

@author: doreper

Loads up configuration files for python code(i.e. cluster distribute), sets up global properties based on those configuration files (see config folder for default config checked into source).
These will be shared throughout a given run, but can be reset dynamically by calling buildConfigMaps at runtime repeatedly, pointing to a different file each time.
Also reads in the column header files specifying which columns to use from our data files, sets up some of the necessary output directories for analysis,
intializes logging system.

'''
import logging.config;
import os;
import yaml;

class configwrap(object):
    """puls items out of configParser into a data structure with fields named by the keys in he config file """
    def __init__(self, configFnames):
        yamlMap = yamlReadConfigs(configFnames);
        for k,v in yamlMap.items(): 
            setattr(self, k, v)
            
def yamlReadConfigs(configFnames):
    endMap = {}
    for configFname in configFnames:
        configFle = open(configFname)
        yamlVals = yaml.load(configFle);
        endMap.update(yamlVals)
        configFle.close()
    return(endMap)

def str2bool(v):
    return v.lower() in ("true", "t", "1", "y", "yes")

allProps = None
headersByFile = None
outputDirs = None
configFileSource = None

def buildConfigMaps(configFilenames):
    """ Reads in config files, in order of configFileNamames. 
    The ability to parse mutlpile files lends itself to the ability to override selected properties.
    Last file in wins overwrites. 
    As an example,  this could be helpful when there is a base configuration, and a subset of configuration changes 
    which only apply to running locally
    
    """
    ##these are global so that we can set the scope outside of this method.
#    global headersByFile
    global allProps
    global configFileSource
    
    configFileSource = configFilenames
    allProps         = configwrap(configFilenames)
    logging.config.fileConfig(allProps.pythonLogConfigFile)
       
#intentionally does not call dataStructureWriting- we dont want circular dependency 
def buildIfDoesNotExist(self, adir):
    if not os.path.exists(adir):
        os.mkdir(adir)
