/**
 * 
 */
package edu.unc.doreper;

import java.io.File;
import java.util.Collections;

import org.apache.log4j.Logger;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.Cigar;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.SequenceUtil;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;

public class AlignmentWriter
{ 
	
	public static final int POS_STRAND = 0;
	public static final int NEG_STRAND = 1;
	public static final int[] strands = {POS_STRAND, NEG_STRAND};
	public static final boolean[] strandsBool = {false,true};
	
	public static final int NONE = 0;
	public static final int LEFT = 1;
	public static final int DIAG = 2;
	public static final int UP   = 3;
	
	static Logger log = Logger.getLogger(AlignmentWriter.class.getCanonicalName());
	static final boolean DEBUG = log.isDebugEnabled();
	static final boolean INFO = log.isInfoEnabled();

	
	public static class Intpair
	{
		static Intpair intPairSingleton = new Intpair();
		static
		{
			intPairSingleton.start = -10000;
		}
		
		
		int start;
		int maxval;
		public Cigar cigar;
		public String cigarString;
		public int numMM;
		public int numInsDel;
		public int end;
	}
	
	public static final class AlignResults
	{
		//TODO reimplement as 3d aray
		public TIntObjectHashMap<Intpair>[][] dataPerRead;
		public AlignResults(int numReads)
		{
			dataPerRead = new TIntObjectHashMap[numReads][2];
			for(int r=0;r<numReads; r++)
			{
				for(int s : strands)
				{
					dataPerRead[r][s] = new TIntObjectHashMap<Intpair>();
				}
			}
		}
	}

	static void formPairAndWriteToBam(final AlignResults alignResultsR,
									  final AlignResults alignResultsL,
									  final byte[] rbases,
											  final byte[] lbases,
											  final int r,
											  final int maxMMpair,
											  final int maxInDelPair,
											  final int maxInsertSize,
											  final int refGenomeLen,
											  final BAMFileWriter outbam)
	{
		if (DEBUG)
		{
			log.debug("assembling read " + r);
			log.debug("length: " + rbases.length);
			log.debug("num pos right read align: "+alignResultsR.dataPerRead[0][Alignment2.POS_STRAND].size());
			log.debug("num neg right read align: "+alignResultsR.dataPerRead[0][Alignment2.NEG_STRAND].size());
		}
		
		int alignCounter =0;
		for(int s: Alignment2.strands)
		{
			final TIntObjectHashMap<Intpair> startsL = alignResultsL.dataPerRead[0][s];
			final TIntObjectHashMap<Intpair> startsR = alignResultsR.dataPerRead[0][s];
			final int startsL_size = startsL.size();
			final int startsR_size = startsR.size();
			if(startsL_size > 0 && startsR_size>0)
			{
				alignCounter=findPairingForReadAndStrand(refGenomeLen,
						startsR,
						startsL, 
						maxMMpair, maxInDelPair, maxInsertSize,
						r,
						s, 
						rbases,
						lbases,
						outbam, 
						alignCounter);
			}
		}
		
		if(alignCounter==0)
		{
			boolean foundAlignment = false; 
			for(int s :Alignment2.strands)
			{
				final TIntObjectHashMap<Intpair> startsR = alignResultsR.dataPerRead[0][s];
				final int startsR_size = startsR.size();
				final TIntObjectIterator<Intpair> riter = startsR.iterator();
				for(int i=0; i<startsR_size; i++)
				{
					foundAlignment=true;
					SAMRecord samrecR = buildUnmappedSamRec(outbam, r, s, "r", rbases);
					riter.advance();
					final Intpair intpair = riter.value();
					samrecR.setAlignmentStart(intpair.start+1); 
					samrecR.setReadUnmappedFlag(false);
					samrecR.setCigarString(intpair.cigarString);
					samrecR.setReadName(samrecR.getReadName()+alignCounter);
//						samrecR.setInferredInsertSize(-1);

					if(s==Alignment2.NEG_STRAND)
					{
						invertStartPosition(intpair, samrecR, refGenomeLen);
					}
					outbam.addAlignment(samrecR);
					alignCounter++;
				}
			}
			if(!foundAlignment)
			{
				outbam.addAlignment(buildUnmappedSamRec(outbam, r, Alignment2.POS_STRAND, "r", rbases));
			}
		}
	}

//	private static SAMRecord copyUnmappedRead(SAMRecord samRecord, BAMFileWriter outbam)
//	{
//		int strandInt = samRecord.getReadNegativeStrandFlag() ? Alignment2.NEG_STRAND : Alignment2.POS_STRAND;
//		return buildUnmappedSamRec(outbam, strandInt, label, bases)
//	}

	private static int findPairingForReadAndStrand(final int refGenomeLen,
			final TIntObjectHashMap<Intpair> startsR,
			final TIntObjectHashMap<Intpair> startsL,
			final int maxMMpair, final int maxInDelPair,final int maxInsertSize, 
			final int r,
			final int s,
			final byte[] rbases,
			final byte[] lbases,
			BAMFileWriter outbam,
			int counter) 
	{
		final TIntObjectIterator<Intpair> riter = startsR.iterator();
		for(int i=0; i<startsR.size(); i++)
		{
			riter.advance();
			final int startRIndex = riter.key();
			//TODO this is only working because the spacing between sites looking for an alignment is ~1000. Fix this to be more general
			for (int lookup = startRIndex-1; lookup<=startRIndex; lookup++)
			{
				if(!startsL.containsKey(lookup))
				{
					continue;
				}
				final Intpair startRPair = riter.value();
				final Intpair startLPair = startsL.get(lookup);
				final int startL = startLPair.start;
				final int startR = startRPair.start;
				if(startR<startL)
				{
					continue;
				}
				if(!((startR - startL)<maxInsertSize))
				{
					continue;
				}
				if((startRPair.numMM+startLPair.numMM>maxMMpair) || (startRPair.numInsDel+startLPair.numInsDel>maxInDelPair))
				{
					continue;
				}

				SAMRecord samrecL = buildUnmappedSamRec(outbam, r, s, "l", lbases);
				samrecL.setFirstOfPairFlag(true);
				SAMRecord samrecR = buildUnmappedSamRec(outbam, r, s, "r", rbases);
				
				samrecL.setAlignmentStart(startL+1); 
				samrecR.setAlignmentStart(startR+1);

				samrecL.setInferredInsertSize(startL - startRPair.end  + 1);
				samrecR.setInferredInsertSize(startL - startRPair.end  + 1);

				samrecL.setMateUnmappedFlag(false);
				samrecR.setMateUnmappedFlag(false);

				samrecL.setProperPairFlag(true);
				samrecR.setProperPairFlag(true);

				samrecL.setCigarString(startLPair.cigarString);
				samrecR.setCigarString(startRPair.cigarString); 
				if(s==Alignment2.NEG_STRAND)
				{
					invertStartPosition(startLPair, samrecL, refGenomeLen);
					invertStartPosition(startRPair, samrecR, refGenomeLen);
				}

				samrecR.setMateAlignmentStart(samrecL.getAlignmentStart());
				samrecL.setMateAlignmentStart(samrecR.getAlignmentStart());
				samrecL.setReadUnmappedFlag(false);
				samrecR.setReadUnmappedFlag(false);
				
				samrecL.setReadName(samrecL.getReadName()+counter);
				samrecR.setReadName(samrecR.getReadName()+counter);
				
				outbam.addAlignment(samrecL);
				outbam.addAlignment(samrecR);
				counter++;
			}
		}
		return counter;
	}

	private static void invertStartPosition(final Intpair startRPair,
											final SAMRecord samrecR, 
											final int refGenomeLen) 
	{
		final int endr = samrecR.getAlignmentStart() + startRPair.cigar.getReferenceLength()-1;
		samrecR.setAlignmentStart(getRevStrandPos(endr, refGenomeLen));
	}
	
	
	
	static void formPairsAndOutputAll(final int[][] bothRefs,
									  final AlignResults alignResultsR, 
									  final AlignResults alignResultsL,
									  final String readsFileL,
									  final String readsFileR)
	{
		//assumes number of right reads and left reads is the same
		final int numReads = alignResultsR.dataPerRead.length;
		//The length of the reference genome
		final int refGenomeLen = bothRefs[0].length;
		
		final FastaSequenceFile fastaFileL = new FastaSequenceFile(new File(readsFileL), true);
		final FastaSequenceFile fastaFileR = new FastaSequenceFile(new File(readsFileR), true);
		
		final BAMFileWriter outbam = getBamFileWriter(bothRefs);
		
		final AlignResults[] bothResults = new AlignResults[]{alignResultsL, alignResultsR};
		final String[] labels = new String[]{"l","r"};
		for(int read=0; read<numReads; read++)
		{
			final byte[][] bothBases = {fastaFileL.nextSequence().getBases(), fastaFileR.nextSequence().getBases()};
			
			for(int leftOrRight=0; leftOrRight<bothResults.length; leftOrRight++)
			{
				final AlignResults alignResult = bothResults[leftOrRight];
				final String label = labels[leftOrRight];
				final byte[] bases = bothBases[leftOrRight];
				
				final SAMRecord unmappedReadSamRecord = buildUnmappedSamRec(outbam, read, label, bases);
				boolean foundMapping = false;
				
				for(int strandInt :Alignment2.strands)
				{
					if(strandInt==Alignment2.NEG_STRAND)
					{
						SequenceUtil.reverseComplement(bases);
					}
					
					final boolean strandBool = Alignment2.strandsBool[strandInt];
					final TIntObjectHashMap<Intpair> mappings = alignResult.dataPerRead[read][strandInt];
					final int numMappings = mappings.size();
					final TIntObjectIterator<Intpair> mappingsIter = mappings.iterator();
					for(int mappingIndex=0; mappingIndex<numMappings; mappingIndex++)
					{
						foundMapping=true;
						mappingsIter.advance();
						final Intpair intpair = mappingsIter.value();
						final SAMRecord samrec = buildUnmappedSamRec(outbam, read, label, bases);
						
						samrec.setReadNegativeStrandFlag(strandBool);
						samrec.setAlignmentStart(intpair.start+1); 
						samrec.setReadUnmappedFlag(false);
						samrec.setCigarString(intpair.cigarString);

						if(strandInt==Alignment2.NEG_STRAND)
						{
							invertStartPosition(intpair, samrec, refGenomeLen);
						}
						
						outbam.addAlignment(samrec);
					}
				}
				if(!foundMapping)
				{
					outbam.addAlignment(unmappedReadSamRecord);
				}
			}
		}
		outbam.close();
	}
	


	static BAMFileWriter getBamFileWriter(final int[][] bothRefs) 
	{
		//TODO: fix this
		final File outFile = new File((String) Alignment5.allProps.get("outputDir"), "javaNew.bam");
		
		final BAMFileWriter outbam = new BAMFileWriter(outFile);
		outbam.setSortOrder(SAMFileHeader.SortOrder.unsorted, false);
		final SAMFileHeader samfh = formHeader(bothRefs);
		outbam.setHeader(samfh);
		return outbam;
	}

	private static SAMRecord buildUnmappedSamRec(BAMFileWriter outbam, 
			int read,
			final String label,
			final byte[] bases) 
	{
		final SAMRecord samrec = new SAMRecord(outbam.getFileHeader());
		samrec.setReadName(read+label);
		samrec.setFirstOfPairFlag(true);
		samrec.setReadPairedFlag(true);
		samrec.setReadUnmappedFlag(true);
		samrec.setReadString(new String(bases));
		samrec.setNotPrimaryAlignmentFlag(true);
		return samrec;
	}

	private static SAMRecord buildUnmappedSamRec(BAMFileWriter outbam, 
			int read,
			int s,
			final String label,
			final byte[] bases) 
	{
		String readName = read+label;
		return buildUnmappedSamRec(outbam, s, bases, readName);
	}

	private static SAMRecord buildUnmappedSamRec(BAMFileWriter outbam,
		int s,
		final byte[] bases,
		String readName)
	{
		final SAMRecord samrec = new SAMRecord(outbam.getFileHeader());

		final byte[] lbasesCopy = bases.clone();
		if(s==Alignment2.NEG_STRAND)
		{
			SequenceUtil.reverseComplement(lbasesCopy);
		}
		final boolean strandBool = Alignment2.strandsBool[s];
		samrec.setReadNegativeStrandFlag(strandBool);
	
		samrec.setReadName(readName);
		samrec.setFirstOfPairFlag(true);
		samrec.setReadPairedFlag(true);
		samrec.setReadUnmappedFlag(true);
		samrec.setReadString(new String(lbasesCopy));
		samrec.setNotPrimaryAlignmentFlag(true);
		return samrec;
	}
	
	static void formPairsAndOutput(final int[][] bothRefs,
										   final AlignResults alignResultsR, 
										   final AlignResults alignResultsL,
										   final String readsFileL,
										   final String readsFileR,
										   final int maxMMpair,
										   final int maxInDelPair,
										   final int maxInsertSize) 
	{
		FastaSequenceFile fastaFileL = new FastaSequenceFile(new File(readsFileL), true);
		FastaSequenceFile fastaFileR = new FastaSequenceFile(new File(readsFileR), true);
		
		String outDir = (String) Alignment2.allProps.get("outputDir");
		BAMFileWriter outbam = new BAMFileWriter(new File(outDir, "java.bam"));
		outbam.setSortOrder(SAMFileHeader.SortOrder.unsorted, false);
		SAMFileHeader samfh = formHeader(bothRefs);
		outbam.setHeader(samfh);
		
		for(int r=0; r<alignResultsR.dataPerRead.length; r++)
		{
		
			SAMRecord samrecL = new SAMRecord(samfh);
			samrecL.setReadName("l");
			SAMRecord samrecR = new SAMRecord(samfh);
			samrecR.setReadName("r");
			samrecL.setFirstOfPairFlag(true);
			
			samrecR.setReadPairedFlag(true);
			samrecL.setReadPairedFlag(true);

			final byte[] lbases = fastaFileL.nextSequence().getBases();
			final byte[] rbases = fastaFileR.nextSequence().getBases();
			
			
			
			boolean foundPairing = false;
			samrecL.setReadUnmappedFlag(true);
			samrecR.setReadUnmappedFlag(true);
			if(DEBUG)
			{
				log.debug("assembling read " + r);
				log.debug("length: " + rbases.length);
				log.debug("num pos right read align: "+alignResultsR.dataPerRead[r][Alignment2.POS_STRAND].size());
				log.debug("num neg right read align: "+alignResultsR.dataPerRead[r][Alignment2.NEG_STRAND].size());
			}
			for(int s :Alignment2.strands)
			{
				if(foundPairing)
				{
					break;
				}
				
				final TIntObjectHashMap<Intpair> startsL = alignResultsL.dataPerRead[r][s];
				final TIntObjectHashMap<Intpair> startsR = alignResultsR.dataPerRead[r][s];
				final int startsL_size = startsL.size();
				final int startsR_size = startsR.size();
				if(startsL_size > 0 && startsR_size>0)
				{
					final boolean strandBool = Alignment2.strandsBool[s];
					samrecL.setReadNegativeStrandFlag(strandBool);
					samrecR.setReadNegativeStrandFlag(strandBool);
					
					samrecL.setReadUnmappedFlag(false);
					samrecR.setReadUnmappedFlag(false);
					
					final TIntObjectIterator<Intpair> riter = startsR.iterator();
					for(int i=0; i<startsR_size; i++)
					{
						if(foundPairing)
						{
							break;
						}
						riter.advance();
						final int startRIndex = riter.key();
						for (int lookup = startRIndex-1; lookup<=startRIndex; lookup++)
						{
							if(foundPairing)
							{
								break;
							}
							if(!startsL.containsKey(lookup))
							{
								continue;
							}
							final Intpair startRPair = riter.value();
							final Intpair startLPair = startsL.get(lookup);
							final int startL = startLPair.start;
							final int startR = startRPair.start;
							if(!((startR - startL)<maxInsertSize))
							{
								continue;
							}
							if((startRPair.numMM+startLPair.numMM>maxMMpair) ||
							   (startRPair.numInsDel+startLPair.numInsDel>maxInDelPair))
							{
								continue;
							}
							
							samrecL.setAlignmentStart(startL+1); 
							samrecR.setAlignmentStart(startR+1);
							
							samrecL.setInferredInsertSize(startR - startL );
							samrecR.setInferredInsertSize(startR - startL );
							
							samrecL.setMateUnmappedFlag(false);
							samrecR.setMateUnmappedFlag(false);

							samrecL.setProperPairFlag(true);
							samrecR.setProperPairFlag(true);
							
							samrecL.setCigarString(startLPair.cigarString);
							samrecR.setCigarString(startRPair.cigarString); 
							if(s==Alignment2.NEG_STRAND)
							{
								final int endl = samrecL.getAlignmentStart() + startLPair.cigar.getReferenceLength()-1;
								final int endr = samrecR.getAlignmentStart() + startRPair.cigar.getReferenceLength()-1;
								final int refGenomeLen = bothRefs[0].length-1;
								
								samrecL.setAlignmentStart(getRevStrandPos(endl, refGenomeLen));
								samrecR.setAlignmentStart(getRevStrandPos(endr, refGenomeLen));
								
								SequenceUtil.reverseComplement(lbases);
								SequenceUtil.reverseComplement(rbases);
								
							}
							
							samrecR.setMateAlignmentStart(samrecL.getAlignmentStart());
							samrecL.setMateAlignmentStart(samrecR.getAlignmentStart());
							
							foundPairing = true;
							break;
						}
					}
				}
			}
			if(!foundPairing)
			{
				samrecL.setMateUnmappedFlag(true);
				samrecR.setMateUnmappedFlag(true);
				boolean foundL = false;
				boolean foundR = false;
				
				for(int s :Alignment2.strands)
				{
					final boolean strandBool = Alignment2.strandsBool[s];
					
					if(!foundL && !alignResultsL.dataPerRead[r][s].isEmpty())
					{
						samrecL.setReadNegativeStrandFlag(strandBool);
						final TIntObjectIterator<Intpair> iter = alignResultsL.dataPerRead[r][s].iterator();
						iter.advance();
						final Intpair intpair = iter.value();
						samrecL.setAlignmentStart(intpair.start+1); 
						samrecL.setReadUnmappedFlag(false);
						foundL=true;
						samrecL.setCigarString(intpair.cigarString);

						if(s==Alignment2.NEG_STRAND)
						{
							final int endl = samrecL.getAlignmentStart() + intpair.cigar.getReferenceLength() - 1;
							final int refGenomeLen = bothRefs[0].length;
							samrecL.setAlignmentStart(getRevStrandPos(endl, refGenomeLen));
							SequenceUtil.reverseComplement(lbases);
						}
					}
					if(!foundR &&!alignResultsR.dataPerRead[r][s].isEmpty())
					{
						samrecR.setReadNegativeStrandFlag(strandBool);
						TIntObjectIterator<Intpair> iter = alignResultsR.dataPerRead[r][s].iterator();
						iter.advance();
						final Intpair intpair = iter.value();
						samrecR.setAlignmentStart(intpair.start+1); 
						samrecR.setReadUnmappedFlag(false);
						samrecR.setNotPrimaryAlignmentFlag(true);
						samrecR.setCigarString(intpair.cigarString);
						foundR=true;
						
						if(s==Alignment2.NEG_STRAND)
						{
							final int endr = samrecR.getAlignmentStart() + intpair.cigar.getReferenceLength() - 1;
							final int refGenomeLen = bothRefs[0].length;
							samrecR.setAlignmentStart(getRevStrandPos(endr, refGenomeLen));
							SequenceUtil.reverseComplement(rbases);
						}
					}
					if(foundL && foundR)
					{
						break;
					}
				}
			}
			
			samrecL.setReadString(new String(lbases));
			samrecR.setReadString(new String(rbases));
			outbam.addAlignment(samrecL);
			outbam.addAlignment(samrecR);
		}
		outbam.close();
	}
	

	private static int getRevStrandPos(final int endl, final int refGenomeLen) 
	{
		return refGenomeLen - endl + 1;
	}


	private static SAMFileHeader formHeader(final int[][] bothRefs)
	{
		SAMFileHeader samfh = new SAMFileHeader();
		samfh.setTextHeader("sometext!");
		samfh.addComment("header!");
		samfh.setGroupOrder(SAMFileHeader.GroupOrder.none);
		samfh.setReadGroups(Collections.singletonList(new SAMReadGroupRecord("dummygroup")));
		samfh.setProgramRecords(Collections.singletonList(new SAMProgramRecord("dummyprogramgroup")));
		samfh.setSequenceDictionary(new SAMSequenceDictionary(Collections.singletonList(new SAMSequenceRecord("tbref", bothRefs[0].length))));
		return samfh;
	}
}
