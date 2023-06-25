class Center{
	//Representation of a cluster center
	public double[] cLoc, pLoc;
	public int id;

	public Center(double[] center, double[] point, int i){
		id = i;
		cLoc = new double[center.length];
		for(int j = 0; j < center.length; j++){
			cLoc = center;
		}
		pLoc = new double[point.length];
		for(int j = 0; j < center.length; j++){
			pLoc = point;
		}
	}
}