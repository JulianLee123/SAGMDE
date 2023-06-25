import java.util.*;
import java.io.FileNotFoundException;
import java.io.File;
import java.io.PrintWriter;


public class PointGen {
	//generates synthetic datasets as described in the paper
	public static double max(double a, double b){
		if(a > b)
			return a;
		else
			return b;
	}
	public static double min(double a, double b){
		if(a < b)
			return a;
		else
			return b;
	}

	public static void main(String[] args){
		Random rand = new Random();
		int numCenters = 30, numPoints = 600;
		double[][] points = new double[numPoints][2];
		double[][] centers = new double[numCenters][2];
		for(int i = 0; i < numCenters; i++){
			centers[i][0] = 0.1+rand.nextDouble()*0.8;
			centers[i][1] = 0.1+rand.nextDouble()*0.8;
			double spread = 0.02*rand.nextDouble() + 0.02;//   15-->0.04
			for(int j = 0; j < numPoints/numCenters; j++){
				for(int d = 0; d < points[0].length; d++){
					points[numPoints/numCenters*i+j][d] = min(max(rand.nextGaussian()*spread + centers[i][d],0),1);
				}
			}
		}
		try (PrintWriter writer = new PrintWriter(new File("synthetic.csv"))) {
			StringBuilder sb = new StringBuilder();
			for(int i = 0; i < numPoints; i++){
				//System.out.println(i);
				sb.append(i);
				sb.append(',');
				sb.append(points[i][0]);
				sb.append(',');
				sb.append(points[i][1]);
				sb.append(",,,");
				if(i < centers.length){
					sb.append(centers[i][0]);
					sb.append(",");
					sb.append(centers[i][1]);
				}
				sb.append('\n');
			}
			writer.write(sb.toString());
		} catch (FileNotFoundException e) {
			System.out.println(e.getMessage());
		}
	}
}
