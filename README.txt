Author: Daniel Oreper, UNC, Valdar Lab, 2016-06-13
Contents: Description of code for aligning EXIT experiment reads and for postprocessing to
create summary information about the read counts.
License: The LGPL License, Version 3.0. See LICENSE.txt
Third party licenses: See NOTICE.txt

Disclaimer: This is academic code, written merely to run correctly-- not to be reusable!

======================================================
Prerequisites that are not included in this bundle.
=======================================================
BWA 0.7.12-r1039 or greater, installed on path
SAMTOOLS 1.3 or greater, installed on path
JRE 1.8.0 or greater, installed on path


R  prerequisites:
R 3.2 or greater
R453Plus1Toolbox
Biostrings
yaml
Rsamtools
ggplot2
plyr
doMC
ggplot2
gtools
epitools
registerDoMC
ShortRead
qrqc
data.table
seqLogo
data.table

Python prerequsites: (uncertain if libs come with python)
Python 2.7.11 or greater
yaml
loadConfig
logging

============================================================
Run Code 
=============================================================
To run alignment, using PROJECT_LOCATION/src as working directory,

> python clusterDistribute.py ../config/config.txt

To run alignment on the the cluster (with parameters tuned for UNC kure cluster)

> python clusterDistribute.py ../config/config.txt ../config/configCluster.txt

If cluster params need to be changed see config/configCluster.txt, especially clusterRunCommand property.



When alignment (i.e. clusterDistribute.py) completes, to collate all information, run

> R CMD BATCH '--args ../config/config.txt' postProcessAlignment.R

or

> R CMD BATCH '--args ../config/config.txt ../config/configCluster.txt' postProcessAlignment.R

if on cluster.


After running these scripts, PROJECT_LOCATION/output/postProcess will contain resulting collated information in the form of .csv files.

==============================================================
Debugging logs
=============================================================

If clusterDistribute fails, 
PROJECT_LOCATION/output/clusterLogs contains an R log file for each experiment/block of reads.
Check these files if folders PROJECT_LOCATION/output/blocks/BLOCKNAME fail to contain javaNew.bam files.

If postProcessAlignment fails,
PROJECT_LOCATION/src/postProcessAlignment.ROut will contain the error message. 


==============================================================
Adding experiments
=============================================================
If additional samples need to be processed, 3 steps must be taken.
1) Add the EXIT_****_1.reads, EXIT_****_3.reads, and EXIT_****.qvals files for the sample to PROJECT_LOCATION/dataset/seq
2) Edit PROJECT_LOCATION/dataset/seq/experimentData.csv to include all experiment metadata
3) Edit PROJECT_LOCATION/config/config.txt and PROJECT_LOCATION/config/configCluster.txt
so that the includedExperiments list contains the new sample ids.



================================================================
Overview of Directory Structure
================================================================

PROJECT_LOCATION/src         : R and python source code. Run from this directory

 	 	/dataset     : input files

                /dataset/seq :raw sequencing files,
			      of the form EXIT_****_1.reads, EXIT_****_3.reads, EXIT_****.qvals, experimentData.csv
			      where EXIT_****_1.reads is the right end read
			      EXIT_****_3.reads is the left end read,
			      EXIT_****.qvals is the file containing q values for both reads
			      and  **** is a numeric experimentId.
			      experimentData.csv is a tab separated file containing metadata about each experiment id

		/java/src    : java code implementing banded dynamic programming alignment, given a set of candidate regions
		
		/java/lib    : external jars upon which java code depends, as well as runnable jars of compiled source from /java/src; i.e. jalign5.jar and ParsePositionsSAMSE.jar

		/config      : configuration files, including config.txt, the set of default parameters for alignment, etc.

		/output      : location to which output files are written.

                /output/postprocess : location to which final csv summary files are written.
