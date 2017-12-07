import java.util.ArrayList;

public class Read {

	private String read;
	private ArrayList<Integer> qread;
	private double mapqual;
	private int start;
	private int finish;
	private ArrayList<Integer> cigar;
	private int adjust=0; //Readjustment for insertion, to match ref to read correctly...
	private int totalCigs=0; //Total length covered by cig so far...
	private int currcigsize;
	private String currcigop="M"; //always going to start with a match...
	private boolean indel=false; //Indel present in read... 


	public Read(ArrayList<Integer> cigar, String read, String ASCIIqread, int mapread, int start, int finish, int startsoft, boolean indel){
		this.cigar=cigar;
		this.read=read;
		this.start=start;
		this.finish=finish;
		mapqual=Math.pow(10, -mapread/10);
		adjust=startsoft;
		this.currcigsize= cigar.get(0);
		this.indel=indel;
		//convert all the base qualities from ASCII to Phred just once.
		for (char quals: ASCIIqread.toCharArray()){
			qread.add((int)quals-33);
		}
	}

	public Base getBase(int l, String ref){
		int basesCovered=l-start; //bases of the reference covered on this read...
		int pos=basesCovered+adjust; //basesCovered gives the position relative to the read, adding adjust gives the position relative to the current cigar.
		if(basesCovered<currcigsize-1+totalCigs){ //-1: currcigsize has n elements but these elements are from 0 to n-1.
			if(currcigop=="M"){
				String base= read.substring(pos, pos+1);
				if(base==ref){return new Base(base,1);} //this segment was meant for speedup in other program, consider removing...
				int qbase=qread.get(pos);
				return new Base(base,combineBaseAndMap(qbase));
			} //ignore insertion and deletions, added when pos==currcigsize...
			else if(currcigop=="D"){
				return new Base("D",0); //16/09/2014...meant for the frequency to be correct. Can disqualify poorly covered deviations!
			}
		}
		else if(basesCovered==currcigsize-1+totalCigs){
			totalCigs=totalCigs+currcigsize;
			cigar.remove(0); //done with this cigar...onto either insertion or deletion cigar...
			if(currcigop.equals("M")){
				String base= read.substring(pos, pos+1);
				int qbase=qread.get(pos);
				if(!cigar.isEmpty()){
					if(cigar.get(0)<0){ //deletion...
						currcigsize=Math.abs(cigar.get(0));
						currcigop="D";
						String stars = new String(new char[currcigsize]).replace("\0", "*"); //adds '*' per deleted base
						base=base+stars;
						return new Base(base,combineBaseAndMap(qbase));
					}
					else{ //insertion...once done, move on to next cig and don't put in counter as actual insertion doesn't map to reference. No need to update currcigsize and op either.
						adjust=adjust+cigar.get(0);
						base=base+read.substring(pos+1, pos+cigar.get(0)+1);
						for(int i=0;i<cigar.get(0);i++){qbase=qbase+qread.get(pos);} //only for loop in this part!
						cigar.remove(0);
						if(!cigar.isEmpty()){
							currcigsize= cigar.get(0);
							this.currcigop="M";}
						return new Base(base,Math.round(combineBaseAndMap(qbase)));}
				}
				else{
					return new Base(base,Math.round(combineBaseAndMap(qbase)));}
			}
			else if(currcigop.equals("D")){ //currcigop should never be insertion...It's a deletion area so add nothing...but switch for next round...
				adjust=adjust-currcigsize;
				if(!cigar.isEmpty()){
					currcigsize=cigar.get(0);
					currcigop="M";
				}
				return new Base("D",0); //over deletion area 12/09/2014 previously returning nothing, but that gives the wrong frequency so returning D to signify deletion area...
			}
		}
		else{System.out.println("Houston, we've had a problem");
		return null;}

		return null;
	}

	public int done(int l){		
		if(l==finish){return 0;} //signals last base in the read string. 
		else if(l>finish){return 1;}
		else{return 2;}	
	}

	public double combineBaseAndMap(int qbase){
		//Mapping and base quality are given as one value which requires conversion...
		double basequal=Math.pow(10, -qbase/10);
		double finalQual=basequal+mapqual-basequal*mapqual;
		//and deconversion...
		return Math.round(-10*Math.log10(finalQual));} //2/12/14 discretizing values...(e.g 29.995 -> 30)

	public boolean indel(){
		return indel; //Indel is present in the read...
	}
}