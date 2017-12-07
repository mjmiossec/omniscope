import java.util.ArrayList;
import java.util.List;

public class BaseAtPosition {

	private List<Double> baseQ= new ArrayList<Double>();
	private int reads;
	private boolean indelIntersect;
	
	//private boolean unmapped=false; //15/09/2014 Base present on unmapped reads therefore extremely likely false

	public BaseAtPosition(double baseQ, boolean indelIntersect){
		this.reads=1;
		this.baseQ.add(baseQ);
		this.indelIntersect=indelIntersect; //Does it coincide with an indel?
	}	

	public BaseAtPosition(double baseQ, int reads){
		this.baseQ.add(baseQ);
		this.reads=reads;
	}	

	public List<Double> getBaseQ(){
		return baseQ;
	}

	public int getReads(){
		return reads;
	}

	public void update(double nbaseQ, boolean indelUpdate){
		if(!indelIntersect){indelIntersect=indelUpdate;}
		// + because log...
		baseQ.add(nbaseQ);
		reads++;
	}
	
	public boolean indelIntersect(){
		return indelIntersect;
	}
}