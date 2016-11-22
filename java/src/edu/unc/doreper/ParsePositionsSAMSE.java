package edu.unc.doreper;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ParsePositionsSAMSE
{
	private static int getRevStrandPos(final int endl, final int refGenomeLen) 
	{
		return refGenomeLen - endl - 1;
	}

	public final static int UNMAPPED = 4;
	public final static int PROPER_PAIR = 2;
	//	public final static int UNMAPPED_MATE = 8;
	public final static int NEGATIVE_STRAND = (int) Math.pow(2,(5-1));

	public final static boolean EVERY_OTHER_LINE = false;
	public static void main(String[] args) throws IOException
	{
		String samFile = args[0];
		String outFile = args[1];
		int numBases = Integer.parseInt(args[2]);

		BufferedReader reader = new BufferedReader(new FileReader(samFile));

		String line = null;
		int counter = 0; //read every other line

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		reader.readLine();
		reader.readLine();
		
		
		while((line=reader.readLine())!=null)
		{
			counter++;
			writer.write(Integer.toString(counter));
			writer.write(":");
			String[] tokens = line.split("\\s+");
			int flag = Integer.parseInt(tokens[1]);
			if((flag & UNMAPPED)!=UNMAPPED)
			{
				int position = Integer.parseInt(tokens[3]) - 1; //converting to 0 indexing
				int seqLen = tokens[9].length();

				//if the right hand read is on the negative strand
				writer.write("+");
				writer.write(Integer.toString(position));
				writer.write(",");

				position  = getRevStrandPos(position+seqLen-1, numBases);
				writer.write("-");
				writer.write(Integer.toString(position));
				//the right hand side is flipped with respect to the overall read. 


				writeAlternates(writer, tokens, numBases, seqLen);
			}
			else
			{
				writer.write("NA");
			}
			writer.write("\n");
		}
		writer.close();
	}

	private static void writeAlternates(BufferedWriter writer, String[] tokens, int numBases, int seqLen) throws IOException
	{
		String[] alternates = null;
		for(String token: tokens)
		{
			if(token.startsWith("XA:Z:"))
			{
				alternates = token.substring(5).split(";");
				break;
			}
		}
		if(alternates!=null)
		{
			for(String alternate: alternates)
			{
				writer.write(",");
				String[] altInfo = alternate.split(",");
				String altPosString = altInfo[1];
				int position = Integer.parseInt(altPosString.substring(1))-1;
				writer.write("+");
				writer.write(Integer.toString(position));

				writer.write(",");
				position  = getRevStrandPos(position+seqLen-1, numBases);
				writer.write("-");
				writer.write(Integer.toString(position));
			}
		}
		else
		{
			return;
		}
	}
}	
