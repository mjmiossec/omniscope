public class Base {
	public String base;
	public double qual;

	public Base(String base, double qual){
		this.base=base;
		if(qual<1){this.qual=1;}
		else{this.qual=qual;}
	}

	public String getBase(){
		return base;
	}

	public double getQual(){
		return qual;
	}
}
