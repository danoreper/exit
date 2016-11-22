/**
 * 
 */
package edu.unc.doreper;

import edu.unc.doreper.Alignment5.Pair;
import edu.unc.doreper.AlignmentWriter.AlignResults;
import edu.unc.doreper.AlignmentWriter.Intpair;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import au.com.bytecode.opencsv.CSVReader;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.CigarUtil;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.util.SequenceUtil;

import org.apache.log4j.Logger;
import org.yaml.snakeyaml.Yaml;

public class Alignment5 
{
	
	static Logger log = Logger.getLogger(Alignment5.class.getCanonicalName());
	static final boolean DEBUG = log.isDebugEnabled();
	static final boolean INFO = log.isInfoEnabled();
	
	static byte negOneByte = (byte)'-';
	
	public static class Pair<R,S>
	{
		public R r;
		public S s;
	}
	
	
	//TODO move into main
	static Map allProps;
	private static TIntIntHashMap smallIntToOrigByteMap;
	
	public static void main(String[] args) throws IOException
	{				
		log.info("running jalign5");
		allProps= loadProps(args);
		
		log.debug(allProps);
		
		String sep = File.pathSeparator;
		
		String referenceFilename  = new File((String)allProps.get("dataDir"), (String)allProps.get("tbReferenceGenomeFile")).getAbsolutePath();
		String leftReadsFilename  = new File((String)allProps.get("outputDir"), "trimmedRead3.fasta").getAbsolutePath();
		String rightReadsFilename = new File((String)allProps.get("outputDir"), "trimmedRead1.fasta").getAbsolutePath();
		String parsedFilename     = new File((String)allProps.get("outputDir"), "bwa.sam.parsed").getAbsolutePath();
		
//		String cutSitesFilename = new File((String)allProps.get("dataDir"), "h37rv/derived/cutSites.txt").getAbsolutePath();
		align(referenceFilename, leftReadsFilename, rightReadsFilename, parsedFilename);
	}


	private static Map loadProps(String[] args) {
		Yaml yaml = new Yaml();
		Map propsMap = new HashMap();
		for(String arg: args)
		{
			InputStream input;
			try 
			{
				input = new FileInputStream(new File(arg));
			}
			catch (FileNotFoundException e)
			{
				throw new RuntimeException(e);
			}
			propsMap.putAll((Map) yaml.load(input));
		}
		return propsMap;
	}

	private static class WorkingMem
	{
		private final AlignResults alignResults = new AlignResults(1);
		private final int[] prevRow;
		private final int[] curRow;
		private final int[][] ptrMatrix;
		private final int[] val1;
		private final int[] val2;
		private final int bandWidth;
		private final int maxReadLength;
		
		public WorkingMem(int bandwidth,
						  int maxReadLength)
		{
			this.bandWidth = bandwidth;
			this.maxReadLength = maxReadLength;
			prevRow = new int[bandwidth+maxReadLength+1];
			curRow  = new int[bandwidth+maxReadLength+1];

			ptrMatrix   = new int[bandwidth+maxReadLength+1][bandwidth+maxReadLength];
			val1 = new int[bandwidth+maxReadLength+1];
			val2 = new int[bandwidth+maxReadLength+1];

			Arrays.fill(val1, negOneByte);
			Arrays.fill(val2, negOneByte);
		}
	}
	
	public static int hashcode(int[] arr, int start, int len)
	{
		int hash =0;
		final int len1 = start + len;
		final int len2 = arr.length;
		final int lenSmallest = len1<len2 ? len1 : len2;
		for(int i =start; i<lenSmallest; i++)
		{
			hash<<=2;
			hash+=arr[i];
			//TODO avoid harcode
		}
		return hash;
	}
	
	public static void align(String referenceFilename, 
						     String leftReadsFilename,
						     String rightReadsFilename,
						     String positions) throws IOException
	{
		BufferedReader positionsReader = new BufferedReader(new FileReader(new File(positions)));
		final int[][] bothRefs = getRefSeq(referenceFilename);
		
		Pair<int[], TIntIntHashMap> pair = convertToSmallInts(bothRefs);
		final int[] origByteToSmallIntMap = pair.r;
		//TODO fix static.
		smallIntToOrigByteMap = new TIntIntHashMap();
		for(int orig : pair.s.keys())
		{
			int k = pair.s.get(orig);
			smallIntToOrigByteMap.put(k,orig);
		}
		
		
		final ArrayList<int[]> readsL = getReads(leftReadsFilename);
		final ArrayList<int[]> readsR = getReads(rightReadsFilename);
		final int refGenomeLen = bothRefs[0].length; 

		final FastaSequenceFile fastaFileL = new FastaSequenceFile(new File(leftReadsFilename), true);
		final FastaSequenceFile fastaFileR = new FastaSequenceFile(new File(rightReadsFilename), true);
		
		//TODO move into properties
		final int mmPenalty  = -1;
		final int gapPenalty = -2;
		final int maxInsertLength = 600;
		final int maxReadLength = 500;
		
		final int maxScorePerBase = 1;
		final int step = 1000;

		final int maxMM_l     = 1;
		final int maxInsDel_l = 1;
		
		final int maxMM_r     = 2;
		final int maxInsDel_r = 1;

		final int maxMM_total = 2;
		final int maxInsDel_total = 1;
		
		final int effectiveBandwidthL = maxInsDel_l;
		final int effectiveBandwidthR = maxInsDel_r;
		
		final int threshR = getScoreDiffThresh(mmPenalty, gapPenalty, maxScorePerBase, maxMM_r, maxInsDel_r); 
		final int threshL = getScoreDiffThresh(mmPenalty, gapPenalty, maxScorePerBase, maxMM_l, maxInsDel_l); 
		
		final int bandwidthL = getBandwidth(effectiveBandwidthL, step);
		final int bandwidthR = getBandwidth(effectiveBandwidthR, step);

		WorkingMem workingMemR = new WorkingMem(bandwidthR, maxReadLength);
		WorkingMem workingMemL = new WorkingMem(bandwidthL, maxReadLength);
	
		//assumes number of right reads and left reads is the same
		final int numReads = readsR.size();

		final BAMFileWriter outbam = AlignmentWriter.getBamFileWriter(bothRefs);
		
		long start = System.currentTimeMillis();
		for(int r=0; r<numReads;r++)
		{
			final int[] readR = readsR.get(r);
			
			log.debug("Aligning read "+r);
			TIntArrayList[][] candidatePositionsPerStrand = getCandidatePositionsForRead(positionsReader, readR.length);
			
			fixRead(readR, origByteToSmallIntMap);
			
			final int[] readL = readsL.get(r);
			fixRead(readL, origByteToSmallIntMap);
			
			
			log.debug("alligning R");
//			try{
			alignToReference(bothRefs, 
					readR, 
					mmPenalty, 
					gapPenalty,
					threshR, 
					maxMM_r,
					maxInsDel_r, 
					refGenomeLen,
					candidatePositionsPerStrand,
					workingMemR);
//			}
//			catch(Exception e)
//			{
//				throw new RuntimeException(e+"failed on read "+r);
//			}
			

			log.debug("alligning L");
			
			alignToMappedRead(bothRefs, 
					readL, 
					mmPenalty, 
					gapPenalty,
					threshL, 
					maxMM_l,
					maxInsDel_l,
					maxInsertLength,
					workingMemR.alignResults,
					refGenomeLen,
					workingMemL);

			log.debug("outputing");

			final byte[] readBytesR = fastaFileR.nextSequence().getBases();
//			SequenceUtil.reverseComplement(readBytesR);
			final byte[] readBytesL = fastaFileL.nextSequence().getBases();
			AlignmentWriter.formPairAndWriteToBam(workingMemR.alignResults, 
												  workingMemL.alignResults, 
												  readBytesR,
												  readBytesL, 
												  r, 
												  maxMM_total, 
												  maxInsDel_total,
												  maxInsertLength, 
												  refGenomeLen, 
												  outbam);
		}
//		printAlignResults(arr);
		outbam.close();
		log.info("Total time ms: " + (System.currentTimeMillis() - start));
	}


	private static TIntArrayList[][] getCandidatePositionsForRead(BufferedReader positionsReader, int readLen) throws IOException
	{
		TIntArrayList[] candidatePositionsPerStrand = getCandPositionsPerStrand(positionsReader);
		if(candidatePositionsPerStrand==null)
		{
			return null;
		}
		
		//first dimensino is strand, second dimension is cluster start or end
		TIntArrayList[][] clustersPerStrand = new TIntArrayList[2][];
		for (int s = 0; s < candidatePositionsPerStrand.length; s++)
		{
			 TIntArrayList storForStrand = candidatePositionsPerStrand[s];
			 storForStrand.sort();
			
			 if(storForStrand.isEmpty())
			 {
				 //store empty cluster starts.
				 
			 }
			 else
			 {
				 //first element is cluster starts, second is cluster ends
				 clustersPerStrand[s] = new TIntArrayList[2];
				 TIntArrayList clusterStartsForStrand = new TIntArrayList();
				 TIntArrayList clusterEndsForStrand = new TIntArrayList();
				 
				 final int firstPos = storForStrand.get(0);
				 clusterStartsForStrand.add(firstPos);
				 clusterEndsForStrand.add(firstPos);
				 int numClusters = 1;

				 for(int j=1;j<storForStrand.size();j++)
				 {
					 int candPosition = storForStrand.get(j);
					 if(candPosition-clusterEndsForStrand.get(numClusters-1)>readLen)
					 {
						 clusterStartsForStrand.add(candPosition);
						 clusterEndsForStrand.add(candPosition);
						 numClusters++;
//						 if(numClusters>1)
//						 {
//							 System.out.println("dubcluster");
//						 }
					 }
					 else
					 {
						 clusterEndsForStrand.set(numClusters-1,candPosition);
					 }
				 }
				 clustersPerStrand[s][0]= clusterStartsForStrand;
				 clustersPerStrand[s][1]= clusterEndsForStrand;
			 }
			 
			 
		}
		return clustersPerStrand;
	}


	private static TIntArrayList[] getCandPositionsPerStrand(BufferedReader positionsReader) throws IOException
	{
		String line = positionsReader.readLine().split(":")[1];
		if(line.equals("NA"))
		{
			return null;
		}
		
		TIntArrayList[] candidatePositionsPerStrand = new TIntArrayList[2];
		candidatePositionsPerStrand[0] = new TIntArrayList();
		candidatePositionsPerStrand[1] = new TIntArrayList();
		String[] candidatePositionsForRead =line.split(",");

		for(String posString: candidatePositionsForRead)
		{
			TIntArrayList storForStrand;
			if(posString.charAt(0)=='+')
			{
				storForStrand = candidatePositionsPerStrand[AlignmentWriter.POS_STRAND];
			}
			else
			{
				storForStrand = candidatePositionsPerStrand[AlignmentWriter.NEG_STRAND];
			}
			storForStrand.add(Integer.parseInt(posString.substring(1)));
		}
		return candidatePositionsPerStrand;
	}


	
	private static void fixRead(final int[] readR, int[] intHash)
	{
		for(int i = 0; i<readR.length; i++)
		{
			readR[i]=intHash[readR[i]];
		}
	}


	private static Pair<int[], TIntIntHashMap> convertToSmallInts(final int[][] bothRefs)
	{
		int[] refseq = bothRefs[Alignment2.POS_STRAND];
		int[] refseq_rc = bothRefs[Alignment2.NEG_STRAND];
		
		final TIntIntHashMap hash = new TIntIntHashMap();
		int counter = -1;
		for(int i=0; i<refseq.length;i++)
		{
			int elem = refseq[i];
			if(!hash.containsKey(elem))
			{
				counter++;
				hash.put(elem, counter);
			}
			
		}
		//TODO dont hardcode
		int[] intHash = new int[100];
		for (int original: hash.keys())
		{
			intHash[original]=hash.get(original);
		}
		 
		
		for(int i=0; i<refseq.length;i++)
		{
			{
				int elem = refseq[i];
				refseq[i] = hash.get(elem);
			}

			{
				int elem = refseq_rc[i];
				refseq_rc[i] = hash.get(elem);
			}
		}
		Pair<int[], TIntIntHashMap> pair = new Pair<int[], TIntIntHashMap>();
		pair.r = intHash;
		pair.s = hash;
		return pair;
	}

	private static int getScoreDiffThresh(final int mmPenalty,
										  final int gapPenalty, 
										  final int maxScorePerBase, 
										  final int maxMM,
										  final int maxInsDel) 
	{
		return (maxScorePerBase-mmPenalty) * maxMM + (maxScorePerBase-gapPenalty)*maxInsDel;
	}


	private static int getBandwidth(int effectiveBandwidth, final int step) 
	{
		return effectiveBandwidth + step - 1;
	}
	
	private static void alignToReference(final int[][] bothRefs,
										 final int[] read,
										 final int mmPenalty,
										 final int gapPenalty,
										 final int coarseThresh, 
										 
										 final int maxMM,
										 final int maxInsDel,
										 final int reflen,
										 final TIntArrayList[][] candidatePositionsPerStrand,
									     final WorkingMem workingMem)
	{
		final AlignResults alignResults = workingMem.alignResults;
		alignResults.dataPerRead[0][0].clear();
		alignResults.dataPerRead[0][1].clear();
		if(candidatePositionsPerStrand==null)
		{
			return;
		}
//		final int bandwidth = workingMem.bandWidth;
		final int bandwidth = maxInsDel;
		final int[] prevRow = workingMem.prevRow;
		final int[] curRow = workingMem.curRow;
		final int[][] ptrMatrix = workingMem.ptrMatrix;
		final int[] val1 = workingMem.val1;
		final int[] val2 = workingMem.val2;
		final int readLen = read.length;
		for (int s=0; s<2;s++)
		{
			final TIntArrayList[] clustersForStrand = candidatePositionsPerStrand[s];
			if(clustersForStrand==null)
			{
				continue;
			}
			TIntArrayList clusterStarts = clustersForStrand[0];
			TIntArrayList clusterEnds = clustersForStrand[1];
			if(clusterStarts==null)
			{
				continue;
			}
			
			final int[] refstrand = bothRefs[s];
			final TIntObjectHashMap<Intpair> storageForReadForStrand = alignResults.dataPerRead[0][s];
	
//				l.System.out.println(matchingPositions.size());
//			Intpair prevAlign = Intpair.intPairSingleton;
			for(int r=0; r<clusterStarts.size();r++)
			{
				final int clusterStart = clusterStarts.get(r);
				final int clusterEnd   = clusterEnds.get(r);
				
				int refstart = Math.max(clusterStart - bandwidth,0); 	


//				if(deltaR<=coarseThresh) //TODO careful here
				{
					final Intpair intpair = alignFullBand(read, 
							0,
							readLen,

							refstrand,
							refstart, 
							Math.min(readLen + (clusterEnd - clusterStart) + 2*bandwidth,  reflen - refstart),

							mmPenalty,
							gapPenalty,
							readLen + (clusterEnd - clusterStart) + 2*bandwidth,

							prevRow,
							curRow,
							ptrMatrix,
							val1, 
							val2);

					if(intpair.numMM<=maxMM &&intpair.numInsDel<=maxInsDel)//&&intpair.start!=prevAlign.start)
					{
						storageForReadForStrand.put(refstart, intpair);
//						prevAlign = intpair;
					}
				}
			}
		}
		return;
	}

	private static void alignToMappedRead(final int[][] bothRefs,
		final int[] read,
		final int mmPenalty,
		final int gapPenalty,
		final int coarseThresh, 

		final int maxMM,
		final int maxInsDel,
		final int maxInsertLength,
		final AlignResults prevResults,
		final int reflen,
		final WorkingMem workingMem)
	{
		final AlignResults alignResults = workingMem.alignResults;
		alignResults.dataPerRead[0][0].clear();
		alignResults.dataPerRead[0][1].clear();
		//final int bandwidth = workingMem.bandWidth;
		final int bandwidth = maxInsDel;
		final int[] prevRow = workingMem.prevRow;
		final int[] curRow = workingMem.curRow;
		final int[][] ptrMatrix = workingMem.ptrMatrix;
		final int[] val1 = workingMem.val1;
		final int[] val2 = workingMem.val2;
		final int readLen = read.length;

		for (int s : AlignmentWriter.strands)
		{
			final int[] refstrand = bothRefs[s];
			final TIntObjectHashMap<Intpair> storageForReadForStrand = alignResults.dataPerRead[0][s];

			final TIntObjectHashMap<Intpair> cutSites = prevResults.dataPerRead[0][s];
			final int numRefStarts = cutSites.size();
			final TIntObjectIterator<Intpair> cutSitesIter = cutSites.iterator();

//			Intpair prevAlign = Intpair.intPairSingleton;
			for(int r=0; r<numRefStarts;r++)
			{
				cutSitesIter.advance();
				final int prevStart = cutSitesIter.key();
				final int refstart = Math.max(prevStart - maxInsertLength,0);

				{
					final Intpair intpair = alignFullBand(read, 
						0,
						readLen,

						refstrand,
						refstart, 
						Math.min(maxInsertLength+readLen+maxInsDel,  reflen - refstart),

						mmPenalty,
						gapPenalty,
						maxInsertLength+readLen+maxInsDel,

						prevRow,
						curRow,
						ptrMatrix,
						val1, 
						val2);

					if(intpair.numMM<=maxMM &&intpair.numInsDel<=maxInsDel)//&&intpair.start!=prevAlign.start)
					{
						storageForReadForStrand.put(prevStart, intpair);
//						prevAlign = intpair;
					}
				}
			}
		}
		return;
	}


	private static int[][] getRefSeq(String referenceFilename)
	{
		final int[] refseq    = getReads(referenceFilename).get(0);
		final int[] refseq_rc = getReads(referenceFilename,true).get(0);
		
		final int[][] bothRefs = new int[2][];
		bothRefs[AlignmentWriter.POS_STRAND] = refseq;
		bothRefs[AlignmentWriter.NEG_STRAND] = refseq_rc;
		
		return bothRefs;
	}

	
	private static String toString(int[] arr, int from, int len)
	{ 
		final byte[] refBases = new byte[len];
		for (int i = from; i<from+len;i++)
		{
			refBases[i-from] = (byte)smallIntToOrigByteMap.get(arr[i]);
		}
		return new String(refBases);
	}
	
	
	private static ArrayList<int[]> getReads(String referenceFilename) 
	{
		return getReads(referenceFilename, false);
	}
	private static ArrayList<int[]> getReads(String referenceFilename, boolean reverseComplement) 
	{
		FastaSequenceFile fastaFile = new FastaSequenceFile(new File(referenceFilename), true);
		ReferenceSequence seq = fastaFile.nextSequence();
		
		
		final ArrayList<int[]> reads = new ArrayList<int[]>();
		while((seq!= null))
		{
			reads.add(getIntSeq(seq, reverseComplement));
			seq = fastaFile.nextSequence();
		}
		fastaFile.close();
		return reads;
	}
	
	private static int[] getIntSeq(ReferenceSequence refseq, boolean reverseComplement) 
	{
		final byte[] refBases = refseq.getBases();
		if(reverseComplement)
		{
			SequenceUtil.reverseComplement(refBases);
		}
		int len = refBases.length;
		final int[] refBasesInt = new int[len];
		refBasesInt[0]=	(int)negOneByte;
		for (int i = 0; i<len;i++)
		{
			refBasesInt[i] = (int)refBases[i];
		}
		return refBasesInt;
	}
	
	//I think this is semi-global
	public static int alignFasterBand(final int[] letters1, 
			final int start1,
			final int len1,

			final int[] letters2,
			final int start2, 
			final int len2,

			final int mmPenalty,
			final int gapPenalty,
			final int bandwidth,

			int[] prevRow,
			int[] curRow)
	{
		if(DEBUG)
		{
			log.debug("aligning approx: "+(toString(letters1, start1, len1)));
			log.debug("aligning to "+(toString(letters2, start2, len2)));
		}
		
		int maxVal = Integer.MIN_VALUE;
		int[] tempRow=null;
		
//		Arrays.fill(prevRow,0,len2+1,0);
		//This is what makes it semiglobal
		Arrays.fill(prevRow,0); //TODO remove cur row... probably not necessary.
		Arrays.fill(curRow,0);
		
		for(int i =0; i< len1; i++)
		{
			final int letter1 = letters1[i+start1];
			final int ulimit = i + bandwidth + 1;
			final int ulimitMax = ulimit<len2 ? ulimit : len2;
			final int llimit = i - bandwidth;
			final int llimitMin = llimit<0 ? 0 : llimit;
			curRow[0] = prevRow[0] + gapPenalty;
			for(int j = llimitMin; j<ulimitMax; j++)
			{
				final int match = letter1 == letters2[j+start2] ? 1 : mmPenalty;
				final int k = j + 1;
				final int diagVal = prevRow[k-1] + match;
				final int upVal   = prevRow[k]   + gapPenalty;
				final int leftVal = curRow[k-1]  + gapPenalty;
				if(diagVal>upVal)
				{
					if(diagVal>leftVal)
					{							
						curRow[k]=diagVal;
					}
					else
					{
						curRow[k]=leftVal;
					}
				}
				else
				{
					if(upVal>leftVal)
					{
						curRow[k]=upVal;
					}
					else
					{
						curRow[k] = leftVal;
					}
				}
			}
//			System.out.println(Arrays.toString(curRow));
			tempRow = prevRow;
			prevRow = curRow;
			curRow = tempRow;
		}
		for(int k=0;k<len2+1;k++)
		{
			final int val = prevRow[k];
			if(val>maxVal)
			{
				maxVal=val;
			}
		}
		//System.out.println(maxi);
		//System.out.println(maxj);
		//System.out.println(maxVal);

//		System.out.println(start2+":"+maxVal);
		return (maxVal);
	}
	
	//I think this is semi-global
		public static int alignFaster(final int[] letters1, 
				final int start1,
				final int len1,

				final int[] letters2,
				final int start2, 
				final int len2,

				final int mmPenalty,
				final int gapPenalty,
				final int bandwidth,

				int[] prevRow,
				int[] curRow)
		{
			if(DEBUG)
			{
				log.debug("aligning approx: "+(toString(letters1, start1, len1)));
//				log.debug("aligning to "+(toString(letters2, start2, len2)));
			}
			
			int maxVal = Integer.MIN_VALUE;
			int[] tempRow=null;
			
//			Arrays.fill(prevRow,0,len2+1,0);
			//This is what makes it semiglobal
			Arrays.fill(prevRow,0);
//			Arrays.fill(curRow,0);
			
			for(int i =0; i< len1; i++)
			{
				final int letter1 = letters1[i+start1];
				curRow[0] = prevRow[0] + gapPenalty;
				for(int j = 1; j<len2; j++)
				{
					final int match = letter1 == letters2[j+start2] ? 1 : mmPenalty;
					final int k = j + 1;
					final int diagVal = prevRow[k-1] + match;
					final int upVal   = prevRow[k]   + gapPenalty;
					final int leftVal = curRow[k-1]  + gapPenalty;
					if(diagVal>upVal)
					{
						if(diagVal>leftVal)
						{							
							curRow[k]=diagVal;
						}
						else
						{
							curRow[k]=leftVal;
						}
					}
					else
					{
						if(upVal>leftVal)
						{
							curRow[k]=upVal;
						}
						else
						{
							curRow[k] = leftVal;
						}
					}
				}
//				System.out.println(Arrays.toString(curRow));
				tempRow = prevRow;
				prevRow = curRow;
				curRow = tempRow;
			}
			for(int k=0;k<len2+1;k++)
			{
				final int val = prevRow[k];
				if(val>maxVal)
				{
					maxVal=val;
				}
			}
			//System.out.println(maxi);
			//System.out.println(maxj);
			//System.out.println(maxVal);

//			System.out.println(start2+":"+maxVal);
			return (maxVal);
		}
	
	public static int alignFasterBandLocal(final int[] letters1, 
			final int start1,
			final int len1,

			final int[] letters2,
			final int start2, 
			final int len2,

			final int mmPenalty,
			final int gapPenalty,
			final int bandwidth,

			int[] prevRow,
			int[] curRow)
	{
		if(DEBUG)
		{
			log.debug("aligning approx: "+(toString(letters1, start1, len1)));
//			log.debug("aligning to "+(toString(letters2, start2, len2)));
		}
		
		int maxVal = 0;
		int[] tempRow=null;
		
//		Arrays.fill(prevRow,0,len2+1,0);
		Arrays.fill(prevRow,0);
		
		for(int i =0; i< len1; i++)
		{
			final int letter1 = letters1[i+start1];
			final int ulimit = i + bandwidth + 1;
			final int ulimitMax = ulimit<len2 ? ulimit : len2;
			final int llimit = i - bandwidth;
			final int llimitMin = llimit<0 ? 0 : llimit;

			for(int j = llimitMin; j<ulimitMax; j++)
			{
				final int k= j +1;
				final int match = letter1 == letters2[j+start2] ? 1 : mmPenalty;
				final int diagVal = prevRow[k-1] + match;
				final int upVal   = prevRow[k]   + gapPenalty;
				final int leftVal = curRow[k-1]  + gapPenalty;
				if(diagVal>0)
				{
					if(diagVal>upVal)
					{
						if(diagVal>leftVal)
						{
							curRow[k]=diagVal;
							if(diagVal>maxVal)
							{
								maxVal=diagVal;
							}
						}
						else
						{
							curRow[k]=leftVal;
							if(leftVal>maxVal)
							{
								maxVal = leftVal;
							}
						}
					}
					else
					{
						if(upVal>leftVal)
						{
							curRow[k]=upVal;
							if(upVal>maxVal)
							{
								maxVal = upVal;
							}
						}
						else
						{
							curRow[k] = leftVal;
							if(leftVal>maxVal)
							{
								maxVal = leftVal;
							}
						}
					}
				}
				else
				{
					if(upVal>0)
					{
						if(upVal>leftVal)
						{
							curRow[k] = upVal;
							if(upVal>maxVal)
							{
								maxVal = upVal;
							}
						}
						else
						{
							curRow[k] = leftVal;
							if(leftVal>maxVal)
							{
								maxVal = leftVal;
							}
						}
					}
					else
					{
						if(leftVal>0)
						{
							curRow[k]=leftVal;
							if(leftVal>maxVal)
							{
								maxVal = leftVal;
							}
						}
						else
						{
							curRow[k] = 0;
						}
						//otherwise do nothing
					}
				}
			}
//			System.out.println(Arrays.toString(curRow));
			tempRow = prevRow;
			prevRow = curRow;
			curRow = tempRow;

		}
		//System.out.println(maxi);
		//System.out.println(maxj);
		//System.out.println(maxVal);

//		System.out.println(start2+":"+maxVal);
		return (maxVal);
	}
	
	public static Intpair alignFullBand(final int[] letters1, 
			final int start1,
			final int len1,

			final int[] letters2,
			final int start2, 
			final int len2,

			final int mmPenalty,
			final int gapPenalty,
			final int bandwidth,

			int[] prevRow,
			int[] curRow,
			final int[][] backptrMat,
			final int[] val1, 
			final int[] val2)
	{
		if(DEBUG)
		{
			log.info("aligning exact"+(toString(letters1, start1, len1)));
			log.info("aligning to "+(toString(letters2, start2, len2)));
		}
		
		Intpair intpair = new Intpair();
		int maxVal = Integer.MIN_VALUE;
		int[] tempRow=null;
		
//		Arrays.fill(prevRow,0,len2+1,0);
		Arrays.fill(prevRow,0);
		Arrays.fill(curRow,0);
		int[] prevPtr = backptrMat[0];
		//TODO get rid of this... maybe. unclear.
		for(int i=0; i<len1; i++)
		{
			Arrays.fill(backptrMat[i],0);
		}
		Arrays.fill(prevPtr,0);
		int maxi = 0;
		int maxj = 0;
		//TODO look into getting rid of this.
		Arrays.fill(val1, 0);
		Arrays.fill(val2, 0);
		
		for(int i =0; i< len1; i++)
		{
			final int[] ptrRow = backptrMat[i+1];
			final int letter1 = letters1[i+start1];
			final int ulimit = i + bandwidth + 1;
			final int ulimitMax = ulimit<len2 ? ulimit : len2;
			final int llimit = i - bandwidth;
			final int llimitMin = llimit<0 ? 0 : llimit;
			
			curRow[0] = prevRow[0] + gapPenalty;
			ptrRow[0] = AlignmentWriter.UP;
			for(int j = llimitMin; j<ulimitMax; j++)
			{
				final int match = letter1 == letters2[j+start2] ? 1 : mmPenalty;
				final int k = j + 1;
				final int diagVal = prevRow[k-1] + match;
				
				final int upVal   = j<ulimitMax-1 ? prevRow[k]   + gapPenalty : Integer.MIN_VALUE;
				final int leftVal = j>llimitMin ? curRow[k-1]  + gapPenalty : Integer.MIN_VALUE;
				if(diagVal>upVal)
				{
					if(diagVal>leftVal)
					{							
						ptrRow[k] = AlignmentWriter.DIAG;
						curRow[k]=diagVal;

					}
					else
					{
						ptrRow[k] = AlignmentWriter.LEFT;
						curRow[k]=leftVal;
					}
				}
				else
				{
					if(upVal>leftVal)
					{
						ptrRow[k] = AlignmentWriter.UP;
						curRow[k]=upVal;
					}
					else
					{
						ptrRow[k] = AlignmentWriter.LEFT;
						curRow[k] = leftVal;
					}
				}
			}
			//swap prev and current row.
			tempRow = prevRow;
			prevRow = curRow;
			curRow = tempRow;
		}
		maxi = start1 + len1-1;
		int k = maxi+1;
		final int ulimit = k + bandwidth + 1;
		final int ulimitMax = ulimit<len2 ? ulimit : len2;
		final int llimit = k - bandwidth;
		final int llimitMin = llimit<0 ? 0 : llimit;
		
		for(int j=llimitMin;j<ulimitMax;j++)
		{
			final int val = prevRow[j];
			if(val>maxVal)
			{
				maxVal=val;
				maxj=j;
			}
		}
		maxj = maxj -1;
//		Cigar cigarString = new Cigar();
//		CigarUtil.cigarStringFromArray(arg0)
//		cigarString.add(new CigarElement(length, operator))
		
		int i = maxi; 
		int j = maxj;
		int backPtr = backptrMat[i+1][j+1];
		int v= val1.length-1;
		val1[v] = letters1[start1+i];
		val2[v] = letters2[start2+j];
		
	
		int numMM =0;
		int numInsDel = 0;
		LinkedList<CigarElement> cigarElems = new LinkedList<CigarElement>();
		boolean validBackPtr = backPtr!=AlignmentWriter.NONE;
		
		while(validBackPtr)
		{
			if(DEBUG)
			{
				log.info((backPtr));
			}
//			System.out.println(start1+i);
//			System.out.println(start2+j);
//			System.out.println();
//			
			if(backPtr==AlignmentWriter.DIAG)
			{
				
//				System.out.println("DIAG");
				
				val1[v] = letters1[start1+i];
				val2[v] = letters2[start2+j];
				i = i-1;
				j = j-1;
				
				final boolean iseq = val1[v]==val2[v];
				CigarElement elem;
				if(iseq)
				{
					elem = new CigarElement(1,CigarOperator.EQ);
				}
				else
				{
					numMM++;
					elem = new CigarElement(1,CigarOperator.X);
				}
				
				cigarElems.addFirst(elem);
			}
			else if(backPtr == AlignmentWriter.UP)
			{
//				System.out.println("UP");
				val1[v] = letters1[start1+i];
				val2[v] = negOneByte;
				i = i -1;
				
				numInsDel++;
				CigarElement elem = new CigarElement(1, CigarOperator.I);
				cigarElems.addFirst(elem);
				
			}
			else
			{
//				System.out.println("LEFT");
				val2[v] = letters2[start2+j];
				val1[v] = negOneByte;
				j = j - 1;
				
				numInsDel++;
				CigarElement elem = new CigarElement(1, CigarOperator.D);
				cigarElems.addFirst(elem);
			}
			v--;
			if(DEBUG)
			{
				log.info("i:"+i +", j:"+j);
			}
			int nextBackPtr = backptrMat[i+1][j+1];
			validBackPtr = nextBackPtr!=AlignmentWriter.NONE;
			if(!validBackPtr)
			{
				//Need to back of the last backward operation, as it will land us on the -1 row or column.
				if(backPtr==AlignmentWriter.DIAG)
				{
					i+=1;
					j+=1;
				}
				else if(backPtr==AlignmentWriter.LEFT)
				{
					j+=1;
				}
				else if(backPtr==AlignmentWriter.UP)
				{
					i+=1;
				}
				break;
			}
			else
			{
				backPtr = nextBackPtr;
			}
		}
		
		if(cigarElems.size()<len1)
		{
//			alignFullBand( letters1, 
//				start1,
//				len1,
//
//				letters2,
//				start2, 
//				len2,
//
//				mmPenalty,
//				gapPenalty,
//				bandwidth,
//
//				prevRow,
//				curRow,
//				 backptrMat,
//				val1, 
//				val2);
			throw new RuntimeException("bogus output-cigar too short");
		}
		v+=1;
		char[] cigarChar = CigarUtil.cigarArrayFromElements(cigarElems);
//		System.out.println(cigarChar);
		String cigarString=null;
		try
		{
		 cigarString = CigarUtil.cigarStringFromArray(cigarChar);
			cigarChar = CigarUtil.cigarArrayFromString(cigarString);
		}
		catch(Exception e)
		{
//			alignFullBand(letters1, start1, len1, letters2, start2, len2, mmPenalty, gapPenalty, bandwidth, prevRow, curRow, backptrMat, val1, val2);
			throw new RuntimeException(e);

		}
	
		Cigar cigar = new Cigar(cigarElems);
		
		if(DEBUG)
		{
			//TODO form cigar string here
			final byte[] s1 = new byte[val1.length -v];
			final byte[] s2 = new byte[val1.length -v];
			for(int x = 0; x<s1.length;x++)
			{
				s1[x] = (byte)val1[x+v];
				s2[x] = (byte)val2[x+v];
			}
			log.debug(cigarString);
			log.debug(new String(s1));
			log.debug(new String(s2));
			log.debug(start2+j);
			log.debug(start2+maxj);
			log.debug("delta: " +(len1-1 - maxVal));
		}
		final int subjectStart = start2+j;
//		final int matchLen     = maxj -j;
		intpair.start = subjectStart;
		intpair.end = maxj;
		intpair.maxval = maxVal;
		intpair.cigar  = cigar;
		intpair.numMM = numMM;
		intpair.numInsDel = numInsDel;
		intpair.cigarString = cigarString;
		return (intpair);
	}	
	
	public static Intpair alignFullBandLocal(final int[] letters1, 
			final int start1,
			final int len1,

			final int[] letters2,
			final int start2, 
			final int len2,

			final int mmPenalty,
			final int gapPenalty,
			final int bandwidth,

			int[] prevRow,
			int[] curRow,
			final int[][] backptrMat,
			final int[] val1, 
			final int[] val2)
	{
		if(DEBUG)
		{
			log.debug("aligning exact"+(toString(letters1, start1, len1)));
//			log.debug("aligning to "+(toString(letters2, start2, len2)));
		}
		
		Intpair intpair = new Intpair();
		int maxVal = 0;
		int[] tempRow=null;
		
//		Arrays.fill(prevRow,0,len2+1,0);
		Arrays.fill(prevRow,0);
//		Arrays.fill(curRow,0);
		int[] prevPtr = backptrMat[0];
		Arrays.fill(prevPtr,0);
		int maxi = 0;
		int maxj = 0;
//		Arrays.fill(val1, 0);
//		Arrays.fill(val2, 0);
		
		for(int i =0; i< len1; i++)
		{
			final int[] ptrRow = backptrMat[i+1];
			final int letter1 = letters1[i+start1];
			final int ulimit = i + bandwidth + 1;
			final int ulimitMax = ulimit<len2 ? ulimit : len2;
			final int llimit = i - bandwidth;
			final int llimitMin = llimit<0 ? 0 : llimit;

			for(int j = llimitMin; j<ulimitMax; j++)
			{
				final int match = letter1 == letters2[j+start2] ? 1 : mmPenalty;
				final int k = j + 1;
				final int diagVal = prevRow[k-1] + match;
				final int upVal   = prevRow[k]   + gapPenalty;
				final int leftVal = curRow[k-1]  + gapPenalty;
				if(diagVal>0)
				{
					if(diagVal>upVal)
					{
						if(diagVal>leftVal)
						{							
							ptrRow[k] = AlignmentWriter.DIAG;
							curRow[k]=diagVal;
							if(diagVal>maxVal)
							{
								maxi = i;
								maxj = j;
								maxVal=diagVal;
							}
							
						}
						else
						{
							ptrRow[k] = AlignmentWriter.LEFT;
							curRow[k]=leftVal;
							if(leftVal>maxVal)
							{
								maxi = i;
								maxj = j;
								
								maxVal = leftVal;
							}
						}
					}
					else
					{
						if(upVal>leftVal)
						{
							ptrRow[k] = AlignmentWriter.UP;
							curRow[k]=upVal;
							if(upVal>maxVal)
							{
								maxi = i;
								maxj = j;
								maxVal = upVal;
							}
						}
						else
						{
							ptrRow[k] = AlignmentWriter.LEFT;
							curRow[k] = leftVal;
							if(leftVal>maxVal)
							{
								maxi = i;
								maxj = j;
								maxVal = leftVal;
							}
						}
					}
				}
				else
				{
					if(upVal>0)
					{
						if(upVal>leftVal)
						{
							ptrRow[k] = AlignmentWriter.UP;
							curRow[k] = upVal;
							if(upVal>maxVal)
							{
								maxi = i;
								maxj = j;
								maxVal = upVal;
							}
						}
						else
						{
							ptrRow[k] = AlignmentWriter.LEFT;
							curRow[k] = leftVal;
							if(leftVal>maxVal)
							{
								maxi = i;
								maxj = j;
								maxVal = leftVal;
							}
						}
					}
					else
					{
						if(leftVal>0)
						{
							ptrRow[k] = AlignmentWriter.LEFT;
							curRow[k]=leftVal;
							if(leftVal>maxVal)
							{
								maxi = i;
								maxj = j;
								maxVal = leftVal;
							}
						}
						else
						{
							ptrRow[k] = AlignmentWriter.NONE;
							curRow[k] = 0;
						}
						//otherwise do nothing
					}
				}
			}
			//swap prev and current row.
			tempRow = prevRow;
			prevRow = curRow;
			curRow = tempRow;
		}
		
//		Cigar cigarString = new Cigar();
//		CigarUtil.cigarStringFromArray(arg0)
//		cigarString.add(new CigarElement(length, operator))
		
		int i = maxi; 
		int j = maxj;
		int backPtr = backptrMat[i+1][j+1];
		int v= val1.length-1;
		val1[v] = letters1[start1+i];
		val2[v] = letters2[start2+j];
		
	
		int numMM =0;
		int numInsDel = len1-1-i;
		LinkedList<CigarElement> cigarElems = new LinkedList<CigarElement>();
		cigarElems.addFirst(new CigarElement(numInsDel, CigarOperator.DELETION));
		
		boolean validBackPtr = backPtr!=AlignmentWriter.NONE;
		while(validBackPtr)
		{
//			System.out.println(start1+i);
//			System.out.println(start2+j);
//			System.out.println();
//			
			if(backPtr==AlignmentWriter.DIAG)
			{
				
//				System.out.println("DIAG");
				
				val1[v] = letters1[start1+i];
				val2[v] = letters2[start2+j];
				i = i-1;
				j = j-1;
				
				final boolean iseq = val1[v]==val2[v];
				CigarElement elem;
				if(iseq)
				{
					elem = new CigarElement(1,CigarOperator.EQ);
				}
				else
				{
					numMM++;
					elem = new CigarElement(1,CigarOperator.X);
				}
				
				cigarElems.addFirst(elem);
			}
			else if(backPtr == AlignmentWriter.UP)
			{
//				System.out.println("UP");
				val1[v] = letters1[start1+i];
				val2[v] = negOneByte;
				i = i -1;
				
				numInsDel++;
				CigarElement elem = new CigarElement(1, CigarOperator.I);
				cigarElems.addFirst(elem);
				
			}
			else
			{
//				System.out.println("LEFT");
				val2[v] = letters2[start2+j];
				val1[v] = negOneByte;
				j = j - 1;
				
				numInsDel++;
				CigarElement elem = new CigarElement(1, CigarOperator.D);
				cigarElems.addFirst(elem);
			}
			v--;
			int nextBackPtr = backptrMat[i+1][j+1];
			validBackPtr = nextBackPtr!=AlignmentWriter.NONE;
			if(!validBackPtr)
			{
				//Need to back of the last backward operation, as it will land us on the -1 row or column.
				if(backPtr==AlignmentWriter.DIAG)
				{
					i+=1;
					j+=1;
				}
				else if(backPtr==AlignmentWriter.LEFT)
				{
					j+=1;
				}
				else if(backPtr==AlignmentWriter.UP)
				{
					i+=1;
				}
				break;
			}
			else
			{
				backPtr = nextBackPtr;
			}
		}
		v+=1;
		
		int numLeadingInsDel = i;
		numInsDel+=numLeadingInsDel;
		cigarElems.addFirst(new CigarElement(numLeadingInsDel, CigarOperator.D));
		char[] cigarChar = CigarUtil.cigarArrayFromElements(cigarElems);
		String cigarString = CigarUtil.cigarStringFromArray(cigarChar);
		
		cigarChar = CigarUtil.cigarArrayFromString(cigarString);
		Cigar cigar = new Cigar(cigarElems);
		
		
		if(DEBUG)
		{
			//TODO form cigar string here
			final byte[] s1 = new byte[val1.length -v];
			final byte[] s2 = new byte[val1.length -v];
			for(int x = 0; x<s1.length;x++)
			{
				s1[x] = (byte)val1[x+v];
				s2[x] = (byte)val2[x+v];
			}
			log.debug(cigarString);
			log.debug(new String(s1));
			log.debug(new String(s2));
			log.debug(start2+j);
			log.debug(start2+maxj);
			log.debug("delta: " +(len1-1 - maxVal));
		}
		final int subjectStart = start2+j;
//		final int matchLen     = maxj -j;
		intpair.start = subjectStart;
		intpair.maxval = maxVal;
		intpair.cigar  = cigar;
		intpair.numMM = numMM;
		intpair.numInsDel = numInsDel;
		intpair.cigarString = cigarString;
		return (intpair);
	}	
}
