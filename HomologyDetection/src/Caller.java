import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
//Samtools library...
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Caller{

	private static SamReader samfile;
	private static SAMRecordIterator iterator;

	private static int last=0;
	private static int skip=0;
	private static Scanner refScan;
	private static String reference;
	private static String refSequence;
	private static String refRollover; //when dealing with indels, a reference can be on two separate lines, thus the next line needs to be recorded...
	//default values which the user can modify starting with posterior probabilities...
	private static int refLength;

	private static String[] seqName;
	private static int[] seqSize;
	private static ArrayList<Read> setofreads;

	private static Read tempRead;
	
	private static HashMap<String, BaseAtPosition> currbases;
	private static boolean IndelInPos; //Does the position contain an indel?
	private static int numberOfIndels=0; //Is there a non-negligible number of indels?

	public static void main(String[] args) throws FileNotFoundException {

		for(int i=0;i<args.length;i++){
			if(args[i].endsWith(".bam") || args[i].endsWith(".sam")){
				//load sam/bam file...
				samfile = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(args[i]));
				//does this serve any purpose when working with metagenomic data?
				List<SAMSequenceRecord> seqDict = samfile.getFileHeader().getSequenceDictionary().getSequences();
				int size=seqDict.size();
				seqName = new String[size];
				seqSize = new int[size];
				for(int k=0;k<size;k++){
					seqName[k]=seqDict.get(k).getSequenceName();
					seqSize[k]=seqDict.get(k).getSequenceLength();
				}
				//
				System.out.println("SAM/BAM file load ->"+args[i]);
			}		
			//load reference. Must be identical to the original reference.
			else if(args[i].endsWith(".fa") || args[i].endsWith(".fna") || args[i].endsWith(".fasta")){
				//load reference FASTQ file
				reference = args[i]; //23/07/13 changed to string instead of file. Used later.
				refScan = new Scanner(new File(reference));
				//FASTA files may vary in line length.
				String refHeader=refScan.nextLine();
				refSequence = refScan.nextLine();
				//length of each line...
				refLength = refSequence.length();
				refRollover = refScan.nextLine();
				//Don't close, will be used later...
				System.out.println("Reference Sequence used -"+refHeader);
			}			
		}	
		System.out.println("launch successful! Teasing out homologous sequences hitching a ride on this reference");

		for(int j=0;j<seqName.length;j++){
			int pos=seqSize[j];
			setupSequence(seqName[j], pos);
			//If there are large spaces between sequences in a reference file...
			if(j>0 && !refRollover.contains(">")){
				while(refScan.hasNextLine() && !refScan.nextLine().contains(">")){
					//A solution to spaces skipped between reference segments...
				}
				if(!refScan.hasNextLine()){break;}
				refSequence=refScan.nextLine();
				//At the end of file to avoid an error message!
				refRollover=refScan.nextLine();
				last=0;
			}

			String refBase=""; //the reference base at a given position
			for(int l=0;l<pos;l++){
				int totalReads=0; // the total number of reads, reset every iteration.
				if(l>0){refBase= getReferenceBase(l,0);} 
				getDeviationAtPosition(l, refBase);
			}
			//
			//
			//				totalreads=0; //the total of reads are reset every iteration...
			//				String checkRef="";
			//				if(l>0){checkRef= getReferenceBase(l,0);} //For each position, this is the reference, avoids being called over and over...
			//
			//				ArrayList <String> st= new ArrayList<String>(); //9/12/13 added here from 3rd loop over the same set of individuals, also redundant.
			//				for(int k=0;k<individual.size();k++){
			//					if(isBam){
			//						if(l==0){individual.get(k).clearRec();} // 28/10/13 removed redundant double if statement...
			//						individual.get(k).getDevForPositionBAM(l, checkRef);
			//					} //this replaces the commented lines...
			//					else{
			//						individual.get(k).getDevForPositionVec(l, checkRef);}
			//					totalreads=totalreads+individual.get(k).indivreads(); //9/12/13 originally in 2nd loop over the same set of individuals. Was redundant.
			//
			//					//stage 1, multiple sequence from multiple individuals...the following two segments can be condensed...this section has been moved from its own for loop to an already existing for loop.
			//					if(!checkRef.equals("N") && !checkRef.equals("S")){ //9/12/13 this if statement was originally inside the loop, more useful outside the loop.
			//						for(String keybase:individual.get(k).allbases()){
			//							if(!st.contains(keybase)&&!keybase.equals(checkRef)){ //12/09/2014...D signifies there's a deletion over that position, but not the start of one...
			//								st.add(keybase);	//8/1/13 changed every instance of keys.split()[1] with the String keybase...this bit needs working on...
			//							}
			//						}
			//					}
			//				}

			//Getting the total of reads...
			//for each separate deviation we do the following...
			//				for(int i=0;i<st.size();i++){
			//					//double observ=1.0; //set all coverage score to 1.0 by default...
			//
			//					boolean NoIndelIntersect=false; //18/09/2014 if false, no indel or deviations intersect with indels...
			//
			//					double score=0;
			//					int reads=0;
			//					String bdev=st.get(i);
			//					List<Double> basecalls = new ArrayList<Double>(); //07/11/14
			//					double currdev=0.0;
			//					for(int r=0;r<individual.size();r++){ //see if this loop can be gotten rid of!
			//						BaseAtPosition bp = individual.get(r).get(bdev);
			//						if(bp!=null){
			//							if((double)individual.get(r).indelInReads()/bp.getReads()>1 && individual.get(r).getCar()==1){NoIndelIntersect=!bp.indelIntersect();} //18/09/2014 if indel nearby and there is no intersection 24/09/2014 if deviation indel caused error, there should be more of the indel... //26/09/2014 added carrier status...
			//							//07/11/2014
			//							basecalls.addAll(bp.getBaseQ());
			//							//currdev=currdev+bp.getBaseQ();
			//							reads=reads+bp.getReads();
			//							//} 25/08/2014
			//						}
			//					}
			//
			//					Collections.sort(basecalls,Collections.reverseOrder()); //07/11/14
			//					for (int z=0;z<basecalls.size();z++) {
			//					    currdev= currdev+basecalls.get(z)*Math.pow(dependency,z); //dependency is factored in here!
			//					}
			//					
			//					score=calculator(currdev,Oprior); //This is our final calculation for pooled.
			//
			//					//if(score>=0.9){ //10/02/14 required to avoid cases where variants are suspected based on accumulated lack of evidence... 25/08/2014
			//
			//					//STAGE 2 begins, multiple sequences from each individual...
			//
			//					String DP4="DP4="; //store reads and their direction...
			//					String DP="DP="; //12/05/14
			//					double poolQS=-10*Math.log10(1-score); //09/09/2014
			//					if(Double.isNaN(poolQS)){ //Sometimes, instead of score being 1, it is 1+some insignificant number, therefore use NaN instead of Infinite catch.
			//						poolQS=999;
			//					}
			//					String AF="AF="; //12/05/14
			//					String QS="QS="; //26/08/14
			//					double jointscore=1.0; //10/12/13 the final score!
			//					double indivscore; //26/08/2014 individual scores to be reported 
			//					//filters...
			//					boolean carrierfreq=false;
			//					boolean noCarrier=true;
			//
			//					for(int y=0;y<individual.size();y++){
			//						double newscore=0;
			//						BaseAtPosition basepos = individual.get(y).get(bdev);
			//
			//						String strands=""; //28/10/2014 moved over from below the following segment
			//						if(basepos!=null){ //10/12/13
			//							//double newdev=basepos.getBaseQ();
			//							//double newdev=0.0;
			//							//07/11/14
			//							List<Double> singlebasecalls = basepos.getBaseQ();
			//							Collections.sort(singlebasecalls,Collections.reverseOrder()); //07/11/14
			//							//25/11/14 Need bases that are not bdev
			//							List<Double> singlenondevcalls = new ArrayList<Double>();
			//							//for(String bases:individual.get(y).allbases()){
			//								//if(bases!=bdev){ //ref et al.
			//								if(individual.get(y).haskey(checkRef)){
			//									singlenondevcalls.addAll(individual.get(y).get(checkRef).getBaseQ());		
			//							}
			//								Collections.sort(singlenondevcalls,Collections.reverseOrder()); //2/12/14
			//							//}
			//							//
			//							
			//							
			////							for (int z=0;z<singlebasecalls.size();z++) {
			////							    newdev= newdev+singlebasecalls.get(z)*Math.pow(dependency,z); //dependency is factored in here!
			////							}
			//							//25/11/2014
			//							newscore=genos(singlebasecalls,singlenondevcalls,individual.get(y).indivreads(),score,l);
			//							//25/11/2014
			//							
			//							//obtaining a final probability...
			//							//newscore=calculator(newdev,Sprior,score); //individual calculation
			//							
			//							//beggining of moved segment
			//							if(individual.get(y).getCar()==1){noCarrier=false;}
			//							strands=individual.get(y).refreads(checkRef)+","+basepos.getFor()+","+basepos.getRev()+";";
			//
			//							double freq=(double)basepos.getReads()/individual.get(y).indivreads();
			//							if(freq<=0.05 && individual.get(y).getCar()==1){carrierfreq=true;}
			//							AF=AF+dec.format(freq)+";"; //28/10/2014 IndiLowFreq tells us if the current non-carrier individual has a lowfreq.
			//							//end of moved segment
			//						}
			//						else{
			//
			//						
			//						//beginning of moved segment...
			//						strands=individual.get(y).refreads(checkRef)+",0,0;";
			//						AF=AF+0+";";
			//						//end of moved segment...
			//						} //1/10/2014 excesslowQ can only be applied if the deviation is seen in all individuals even if not much...
			//
			//						//10/12/13 important update! //10/1/14 PPV considered to be 1.
			//						double observ=individual.get(y).coverage();
			//						if(individual.get(y).getCar()==1){ //variant of interest must be present...
			//							double scoreToAdd=newscore+(1-newscore)*observ; //14/08/14
			//							indivscore=scoreToAdd;
			//							jointscore=jointscore*scoreToAdd;
			//							//jointscore=jointscore*(newscore+(1-newscore)*observ);
			//						}
			//						//else if(!indiLowFreqNC){ //28/10/2014 skip if non-carrier low frequency...
			//						else{
			//							if(observ==1){
			//								jointscore=0; //09/04/2014
			//								indivscore=0; //26/08/14
			//							}
			//							else{
			//								jointscore=jointscore*(1-newscore); //09/04/2014
			//								indivscore=1-newscore; //26/08/14
			//							}
			//							//jointscore=jointscore*(1-newscore)*(1-observ);
			//						}
			//						//else{
			//						//	indivscore=1; //4/11/2014 corrected, should be 1 not 999.
			//						//}//carrier status 0.
			//						//10/12/13 important update!
			////						String strands="";
			////						if(basepos!=null){
			////							if(individual.get(y).getCar()==1){noCarrier=false;}
			////							strands=individual.get(y).refreads(checkRef)+","+basepos.getFor()+","+basepos.getRev()+";";
			////							if(individual.get(y).getCar()==1){
			////								if(!excessNegative){//23/09/2014 count only if carrier, for non-carrier, lowfreqnoncarrier flag might be more appropriate.
			////									//if(individual.get(y).indivreads()>10){ only necessary without devSeen?
			////									excessNegative=basepos.getExcess();
			////									//} //1/10/2014
			////								}
			////							} //15/09/2014 too many low quality reads or quality read 0...if excessnegative false check for true. Once true, always true...
			////
			////							double freq=(double)basepos.getReads()/individual.get(y).indivreads();
			////							if(freq<=0.08){if(individual.get(y).getCar()==1){carrierfreq=true;}else if(basepos.getReads()<=5){if(noncarrierfreq){multiplenoncarrierfreq=true;}else{noncarrierfreq=true;}}}
			////							AF=AF+dec.format(freq)+";";
			////						}
			////						else{
			////							strands=individual.get(y).refreads(checkRef)+",0,0;";
			////							AF=AF+0+";";
			////						}
			//						if(!outMode){ //Is this mode quite necessary?
			//							out[y].write(toVCF(j,l,bdev,jointscore,checkRef,DP4+strands, ""+reads, ""+totalreads/reads, QS, poolQS, "")+"\n"); //I added observ here but it doesn't belong here!
			//							//I'm unsure if this section is even still needed...it can be reworked at any rate...	
			//						}
			//						else if(outMode){
			//							DP4=DP4+strands;
			//							DP=DP+individual.get(y).indivreads()+";"; //12/05/14
			//							indivscore=-10*Math.log10(1-indivscore);
			//							if(Double.isInfinite(indivscore) ){indivscore=999;}
			//							QS=QS+dec2.format(indivscore)+";";
			//
			//
			//						}
			//					}
			//					if(outMode){
			//						if(jointscore>=minQual && !noCarrier){
			//							String filter;
			//							if(noCarrier){filter="noCarrier";}else if(NoIndelIntersect){filter="IndelArtefact";}else if(carrierfreq){filter="lowFreqCarrier";}else{filter=".";}
			//							//if(applyFilters && poolQS>=40 && !noCarrier){ //implementing only pool QS for now! 17/09/2014 if deviation within deletion, disqualify
			//								out[0].write(toVCF(j,l,bdev,jointscore,checkRef,DP4, DP, AF, QS, poolQS, filter)+"\n");
			//							//}
			//						} //26/08/14
			//					}
			//					//} 25/08/2014
			//				}
			//			}

		}
	}

	public static void setupSequence(String seqName, int pos){
		//get all the reads that overlap with a sequence in the reference.
		if(iterator!=null){
			iterator.close();}
		setofreads = new ArrayList<Read>();
		iterator=samfile.queryOverlapping(seqName, 0, pos);
		skip=0;
	}

	public static String getReferenceBase(int pos, int extra) throws FileNotFoundException {
		int adj=pos-(last*refLength);
		int iterations=adj/refLength; //doesn't round up, that's the trick...

		if(adj%refLength==0){iterations=iterations-1;}
		for(int i=1;i<=iterations;i++){
			refSequence=refRollover;
			if(refScan.hasNextLine()){
				refRollover=refScan.nextLine();}}
		last=last+iterations; //move this to the equation below??
		String ref=""+refSequence.charAt(pos-(last*refLength)-1); //-1 because index starts with 0. Char to String...	

		for(int i=1;i<=extra;i++){
			//for deletions...check if at the edge of a line so that you can rollover
			if(pos-(last*refLength)-1+i>=refLength){
				ref=ref+refRollover.charAt(pos-((last+1)*refLength)-1+i);
			}
			else{
				ref=ref+refSequence.charAt(pos-(last*refLength)-1+i);}
		}
		return ref.toUpperCase();
	}

	public static void getDeviationAtPosition(int l, String ref){ 
		//all the deviations from that particular position
		clear(); //Need to establish if this is useful in this version.
		//what does tempRead do again? I need to write more comments!

		while((iterator.hasNext() || tempRead!=null) && l==skip){			
			if(tempRead!=null){
				setofreads.add(tempRead);
				tempRead=null;
				if(!iterator.hasNext()){
					break;}
			}

			SAMRecord sr=iterator.next();
			try{//reads that fail basic requirements.
				if(sr.getDuplicateReadFlag()||sr.getReadFailsVendorQualityCheckFlag()||sr.getMappingQuality()==255){
					continue;}
			}catch(java.lang.NullPointerException e){}

			Read r = getRead(sr);
			if(r==null){continue;}
			int start=sr.getAlignmentStart();
			if(start==l){
				setofreads.add(r);}
			else{ 
				tempRead = r;
				skip=start;
				break;
			}		
		}
		for(int x=0;x<setofreads.size();x++){
			Read r = setofreads.get(x);
			//solution to removal redundancy problem...
			if(r.done(l)<2){
				setofreads.remove(x); //if it's on the last position, it won't be needed next iteration
				x--;
				if(r.done(l)==1){
					continue;} //if it's way past the last position, skip!
			} 
			Base base = r.getBase(l,ref);
			boolean isIndel=false; 
			if(base!=null){
				if(base.getBase().length()>1){isIndel=true;} 
				update(base.getBase(),base.getQual(),isIndel, r.indel());} //Is the change an indel? Is there an indel in this read generally? if first true, second should be true.
		}
	}

	public static void clear(){
		currbases.clear(); //clears when a new position is being studied...
		IndelInPos=false;
		numberOfIndels=0;

	}
	
	public static void update(String base, double bqphred, boolean hasIndel, boolean indelInRead){ //Added indel awareness for various functions...18/09/2014 added a value for the presence of an indel somewhere in the read...
		
		if(indelInRead){numberOfIndels++;} //Important not to disqualify deviations because of error indel...
		
		if(currbases.containsKey(base)){
			currbases.get(base).update(bqphred, indelInRead);}
		else{
			if(hasIndel){IndelInPos=hasIndel;} //only needs to be true once to be triggered permanently...
			currbases.put(base, new BaseAtPosition(bqphred, indelInRead));} //Signal whether individual bases share a read with an indel...
	}
	
	public static Read getRead(SAMRecord sr){

		//catching problematic CIGARs... 
		if(sr.getCigarString().equals("*")){return null;}
		//...and determining soft clips length
		int startSoft = 0;
		if(sr.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.S)){
			//temporary solution. D(eletion) is not supposed to be a starting cigar, but I have a doubt about I(nsertion).
			if(sr.getCigar().getCigarElement(1).getOperator().equals(CigarOperator.I) || sr.getCigar().getCigarElement(1).getOperator().equals(CigarOperator.D)){return null;}
			else{startSoft=sr.getCigar().getCigarElement(0).getLength();}
		}
		//storing CIGAR
		ArrayList<Integer> cigar = new ArrayList<Integer>();
		boolean indel=false;
		for (CigarElement cig: sr.getCigar().getCigarElements()){
			CigarOperator cigop = cig.getOperator();
			//signal presence of one or more indels (only needs to be flipped once)...MAKE SURE THIS NEW FORMAT WORKS!
			if(!indel){indel= cigop.isIndel();}
			//get length of each CIGAR, add sign if deletion...
			if(cigop.equals(CigarOperator.M)||cigop.equals(CigarOperator.I)){
				cigar.add(cig.getLength());}
			else if(cigop.equals(CigarOperator.D)){
				cigar.add(-cig.getLength());
			}	
		}
		return new Read(cigar ,sr.getReadString(), sr.getBaseQualityString(), sr.getMappingQuality(), sr.getAlignmentStart(), sr.getAlignmentEnd(), startSoft, indel);
	}
}