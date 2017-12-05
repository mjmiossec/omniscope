import java.io.File;
import java.util.List;
//Samtools 
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Caller {
	
	private static String[] seqName;
	private static int[] seqSize;	

	public static void main(String[] args) {
		
		for(int i=0;i<args.length;i++){
			if(args[i].endsWith(".bam") || args[i].endsWith(".sam")){
				//load BAM file...
				SamReaderFactory samfactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
				SamReader samfile = samfactory.open(new File(args[i]));

				if(!compareSeq(samfile)){
					//Compare the sequences in the .bam files to make sure there is compatibility. Must replace break by a program error.
					System.out.println("bam files have incompatible sequences");
					//break; restore this later...important to take note of this...
				}

			}
			else if(args[i].endsWith(".fa") || args[i].contains("fasta")){
				//load reference FASTQ file
			}
		}

	}
	
	public static boolean compareSeq(SamReader samfile) {
		List<SAMSequenceRecord> seqdict = samfile.getFileHeader().getSequenceDictionary().getSequences();
		int size=seqdict.size();
		if(seqName == null){
			seqName = new String[size];
			seqSize = new int[size];
			for(int k=0;k<size;k++){
				seqName[k]=seqdict.get(k).getSequenceName();
				seqSize[k]=seqdict.get(k).getSequenceLength();
			}
			return true;
		}
		else{
			for(int k=0;k<size;k++){
				if(!seqName[k].equals(seqdict.get(k).getSequenceName()) || seqSize[k]!=seqdict.get(k).getSequenceLength()){
					return false;}
			}
			return true;
		}	
	}

}
