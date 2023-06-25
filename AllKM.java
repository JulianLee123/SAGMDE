import java.util.*;
import java.lang.Math;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class AllKM {

	public static double log2(double x)
	{
		return Math.log(x + 1e-10) / (Math.log(2) + 1e-10);
	}
	
	//finds the distance between two points with the same number of dimensions
	public static double dist(double[] p1, double[] p2){
		float distance = 0;
		for(int i = 0; i < p1.length; i++){
			distance += Math.pow(Math.abs(p1[i]-p2[i]),2);
		}
		return Math.sqrt(distance);
	}

	//forgy initialization for centers of clusters
	public static double[][] forgy(double[][] points, int K){
		Random rand = new Random();
		double[][] centers = new double[K][points[0].length];
		//initialize centers by the forgy method
		boolean[] isCenter = new boolean[points.length];
		for(int i = 0; i < points.length; i++){
			isCenter[i] = false;
		}
		for(int c = 0; c < K; c++){
			int idx = rand.nextInt(points.length);
			while(isCenter[idx])
				idx = rand.nextInt(points.length);
			isCenter[idx] = true;
			//assign center to point
			for(int d = 0; d < points[0].length; d++){
				centers[c][d] = points[idx][d];
			}
		}
		return centers;
	}

	//random partition initialization for centers of clusters
	public static double[][] randomPartition(double[][] points, int K){//assigns points to clusters randomly and then computes the centroid of each cluster as the centers
		//initialization
		Random rand = new Random();
		double[][] centers = new double[K][points[0].length];
		int[] numPts = new int[K];//number of points in the given cluster
		for(int c = 0; c < K; c++){
			for(int d = 0; d < points[0].length; d++){
				centers[c][d] = 0;
			}
			numPts[c] = 0;
		}
		//assign points to clusters; centers will be equal to the sum of all points in the given cluster
		for(int i = 0; i < points.length; i++){
			int cluster = rand.nextInt(K);
			numPts[cluster]++;
			for(int d = 0; d < points[0].length; d++){
				centers[cluster][d] += points[i][d];
				if(Double.isNaN(points[i][d]))
					System.out.println(points[i][d] + " " + i + " " + d);
			}
		}
		for(int c = 0; c < K; c++){
			for(int d = 0; d < points[0].length; d++){
				centers[c][d] /= numPts[c];
			}
		}
		return centers;
	}

	//generate cluster centers at random locations
	public static double[][] randomGeneration(double[][] points, int K){
		//initialization
		Random rand = new Random();
		double[][] centers = new double[K][points[0].length];
		for(int c = 0; c < K; c++){
			for(int d = 0; d < points[0].length; d++){
				centers[c][d] = rand.nextDouble();
			}
		}
		return centers;
	}

	public static double[][] CPCL(double[][] points, int K){//K >= actual K
		double[][] centers = forgy(points,K);//new double[K][points[0].length];
		double[][] oldCenters = new double[centers.length][centers[0].length];
		double[] wins = new double[K];
		int T = 0, totalWins = 0;
		double e = 1, epsilon = 0.00001, lr = 0.001;//learning rate
		for(int c = 0; c < K; c++){
			wins[c] = 1;
			totalWins++;
		}
		while(e > epsilon){
			for(int c = 0; c < centers.length; c++){
				for(int d = 0; d < centers[0].length; d++){
					oldCenters[c][d] = centers[c][d];
				}
			}
			for(int i = 0; i < points.length; i++){
				double minI = 100000;
				int c = 0;
				//find winner
				for(int j = 0; j < K; j++){
					if(minI > wins[j]/totalWins*Math.pow(dist(points[i],centers[j]),2)){
						minI = wins[j]/totalWins*Math.pow(dist(points[i],centers[j]),2);
						c = j;
					}
				}
				Vector<Center> Sc = new Vector<Center>();
				for(int j = 0; j < K; j++){
					if(dist(centers[j],centers[c]) <= dist(centers[c],points[i]) && c != j)
						Sc.add(new Center(centers[j],centers[c],j));
				}
				Collections.sort(Sc, new CenterComparator());
				int Qu = (int)(Sc.size()*Math.min(1.0,lr*wins[c]));
				for(int m = 0; m < Qu; m++){//cooperate
					int j = Sc.get(m).id;//get the center location
					double pu = dist(centers[c],points[i])/Math.max(dist(centers[c],points[i]),dist(centers[j],points[i]));
					for(int d = 0; d < points[0].length; d++)
						centers[j][d] += lr*pu*(points[i][d] - centers[j][d]);
				}
				for(int m = Qu; m < Sc.size(); m++){//penalize
					int j = Sc.get(m).id;//get the center location
					double pu = dist(centers[c],points[i])/dist(centers[j],points[i]);
					for(int d = 0; d < points[0].length; d++)
						centers[j][d] -= lr*pu*(points[i][d] - centers[j][d]);
				}
				for(int d = 0; d < points[0].length; d++)
					centers[c][d] += lr*(points[i][d] - centers[c][d]);
				wins[c]++;
				totalWins++;
				T++;
				e = 0;
				for(int k = 0; k < centers.length; k++)
					e += dist(oldCenters[k],centers[k]);
			}
		}
		return centers;
	}

	public static double[][] KM(double[][] points, int K){
		int maxIt = 1000000;
		double[][] centers = randomPartition(points,K);
		int it = -1;
		while(++it < maxIt){
			int[] pointToCluster = pointsToClusters(points,centers); //index of the center to which point i is closest to
			int[] pointsPerCluster = new int[K];
			for(int i = 0; i < points.length; i++)
			pointsPerCluster[pointToCluster[i]]++;
			//calculate new centers
			double[][] newCenters = new double[K][points[0].length];
			//each center is equal to the sum of all points
			for(int i = 0; i < points.length; i++){
				for(int d = 0; d < points[0].length; d++){
						newCenters[pointToCluster[i]][d] += points[i][d];
				}
			}
			//divide each center by the number of points
			boolean done = true;
			for(int c = 0; c < K; c++){
				for(int d = 0; d < points[0].length; d++){
					newCenters[c][d] /= pointsPerCluster[c];
					if(Math.abs(newCenters[c][d] - centers[c][d]) > 0.0001)
						done = false;
				}
			}
			if(done)
				break;
			centers = newCenters;
		}
		return centers;
	}

	public static double[][] BH(double[][] points, int K){
		//randomPartition initialization
		Random rand = new Random();
		double[][] centers = new double[K][points[0].length];
		int[] pointToCluster = new int[points.length];
		int[] numPts = new int[K];//number of points in the given cluster
		for(int c = 0; c < K; c++){
			for(int d = 0; d < points[0].length; d++){
				centers[c][d] = 0;
			}
			numPts[c] = 0;
		}
		//assign points to clusters; centers will be equal to the sum of all points in the given cluster
		for(int i = 0; i < points.length; i++){
			int cluster = rand.nextInt(K);
			pointToCluster[i] = cluster;
			numPts[cluster]++;
			for(int d = 0; d < points[0].length; d++){
				centers[cluster][d] += points[i][d];
			}
		}
		for(int c = 0; c < K; c++){
			for(int d = 0; d < points[0].length; d++){
				centers[c][d] /= numPts[c];
			}
		}
		//variables for the algorithm
		int[] pointsPerCluster = new int[K];
		for(int i = 0; i < points.length; i++){
			pointsPerCluster[pointToCluster[i]]++;
		}
		int emptyCluster = -1;//if there is an empty cluster, this will store the index of the corresponding center
		double currPerf = KMPerf(points,centers), newPerf = KMPerf(points,centers);
		//recommended values from brown and huntley
		int MAXITER = 4*points.length, numTemp = 200;
		double x = 0.75, b = 0.125, EPSILON = 0.00000000001, m = 0, u = 0, minU = 1000000;
		boolean trial = true;
		//will be recomputed during the trial run
		double T = 100, minT = 0.000001, t = 0.95;

		//for display
		double[][] bestCenters = new double[K][points[0].length];
		double bestPerf = 100000000;
		while(T > minT){
			for(int it = 0; it < MAXITER; it++){
				if(!trial && currPerf < bestPerf){
					bestPerf = currPerf;
					for(int k = 0; k < K; k++){
						for(int d = 0; d < points[0].length; d++){
							bestCenters[k][d] = centers[k][d];
						}
					}
				}
				//shift one point at random to another cluster at random
				int movedPoint = rand.nextInt(points.length), oldCluster = pointToCluster[movedPoint], oldEmptyCluster = emptyCluster;
				pointsPerCluster[oldCluster]--;
				int newCluster = rand.nextInt(K);
				if(emptyCluster != -1){//shift point to empty cluster
					pointToCluster[movedPoint] = emptyCluster;
					newCluster = emptyCluster;
					emptyCluster = -1;
				}
				else{
					pointToCluster[movedPoint] = newCluster;
				}
				pointsPerCluster[newCluster]++;
				if(pointsPerCluster[oldCluster] == 0)
					emptyCluster = oldCluster;
				//calculate new centers
				double[][] newCenters = new double[K][points[0].length];
				for(int i = 0; i < points.length; i++){
					for(int d = 0; d < points[0].length; d++){
						newCenters[pointToCluster[i]][d] += points[i][d];
					}
				}
				for(int c = 0; c < K; c++){
					for(int d = 0; d < points[0].length; d++){
						newCenters[c][d] /= pointsPerCluster[c];
					}
				}
				newPerf = KMPerf(points,newCenters);
				double randDouble = rand.nextDouble();
				if(newPerf - currPerf <= 0 || (!trial && Math.pow(Math.E,(currPerf - newPerf)/T) > randDouble)){//accept new solution
					centers = newCenters;
					currPerf = newPerf;
				}
				else{
					if(trial){
						m++;//cost increase
						u += newPerf - currPerf;
						if(newPerf - currPerf < minU)
						minU = newPerf - currPerf;
					}
					pointToCluster[movedPoint] = oldCluster;
					pointsPerCluster[newCluster]--;
					pointsPerCluster[oldCluster]++;
					emptyCluster = oldEmptyCluster;
				}
			}
			T *= t;
			if(trial){//compute temperatures
				u /= m;
				T = u/Math.log(m/(x*m-(1-x)*(MAXITER-m)));
				minT = -1 * b * u / Math.log(EPSILON);
				t = Math.pow(minT/T,1.0/(double)numTemp);
				trial = false;
			}
		}
		return centers;
	}

	public static double[][] SAKM(double[][] points, int K){
		//randomPartition initialization
		Random rand = new Random();
		double[][] centers = new double[K][points[0].length];
		int[] pointToCluster = new int[points.length];
		int[] numPts = new int[K];//number of points in the given cluster
		for(int c = 0; c < K; c++){
			for(int d = 0; d < points[0].length; d++){
				centers[c][d] = 0;
			}
			numPts[c] = 0;
		}
		//assign points to clusters; centers will be equal to the sum of all points in the given cluster
		for(int i = 0; i < points.length; i++){
			int cluster = rand.nextInt(K);
			pointToCluster[i] = cluster;
			numPts[cluster]++;
			for(int d = 0; d < points[0].length; d++){
				centers[cluster][d] += points[i][d];
			}
		}
		for(int c = 0; c < K; c++){
			for(int d = 0; d < points[0].length; d++){
				centers[c][d] /= numPts[c];
			}
		}
		//variables for the algorithm
		double currPerf = KMPerf(points,centers), newPerf = KMPerf(points,centers);
		int MAXITER = points.length, numTemp = 200;
		double T = 100, t = 0.95;
		//for display
		double[][] bestCenters = new double[K][points[0].length];
		double bestPerf = 100000000;
		boolean done = false;
		while(true){
			if(numTemp-- == 0)
				break;
			for(int it = 0; it < MAXITER; it++){
				if(currPerf < bestPerf){
					bestPerf = currPerf;
					for(int k = 0; k < K; k++){
						for(int d = 0; d < points[0].length; d++){
							bestCenters[k][d] = centers[k][d];
						}
					}
				}
				//pertubation
				int[] newPointToCluster = new int[points.length];
				for(int i = 0; i < points.length; i++){
					newPointToCluster[i] = pointToCluster[i];
				}
				done = true;
				for(int i = 0; i < points.length; i++){
					int Ck = rand.nextInt(K);
					if(newPointToCluster[i] == Ck)
					continue;
					double acceptProb = Math.pow(Math.E,-Math.max(dist(points[i],centers[Ck])-dist(points[i],centers[newPointToCluster[i]]),0)/T);
					if(acceptProb > 0.00001)
					done = false;
					if(acceptProb > rand.nextDouble())
					newPointToCluster[i] = Ck;
				}
				if(done)
					break;
				//calculate new centers
				double[][] newCenters = new double[K][points[0].length];
				for(int i = 0; i < points.length; i++){
					for(int d = 0; d < points[0].length; d++){
						newCenters[newPointToCluster[i]][d] += points[i][d];
					}
				}
				int[] pointsPerCluster = new int[K];
				for(int i = 0; i < points.length; i++){
					pointsPerCluster[newPointToCluster[i]]++;
				}
				for(int c = 0; c < K; c++){
					for(int d = 0; d < points[0].length; d++){
						newCenters[c][d] /= pointsPerCluster[c];
					}
				}
				newPerf = KMPerf(points,newCenters);
				double randDouble = rand.nextDouble();
				if(newPerf - currPerf <= 0 || Math.pow(Math.E,(currPerf - newPerf)/T) > randDouble){//accept new solution
					centers = newCenters;
					pointToCluster = newPointToCluster;
				}
			}
			T *= t;
			if(done)
			break;
		}
		return centers;
	}

	public static double[][] SAGM(double[][] points, int K){
		double[][] centers = randomPartition(points,K);
		//variables for the algorithm
		Random rand = new Random();
		double currPerf = KMPerf(points,centers), newPerf = KMPerf(points,centers);
		//recommended values
		int MAXITER = 2*points.length;
		double T = 0.002, minT = 0.0001, t = 0.95;
		//for display
		while(T > minT){
			for(int it = 0; it < MAXITER; it++){
				//calculate new centers
				double[][] newCenters = new double[K][points[0].length];
				for(int c = 0; c < K; c++){
					for(int d = 0; d < points[0].length; d++){
						newCenters[c][d] = centers[c][d] + rand.nextGaussian()*0.001;
					}
				}
				newPerf = KMPerf(points,newCenters);
				double randDouble = rand.nextDouble();
				if(newPerf - currPerf <= 0 || Math.pow(Math.E,(currPerf - newPerf)/T) > randDouble){//accept new solution
					centers = newCenters;
					currPerf = newPerf;
				}
			}
			T *= t;
		}
		return centers;
	}

	//calculates KMPerf of one run of the distortion equalization perturbation
	public static double DEPerf(double[][] points, double[][] centers, double[] globalMaxVals, double[] globalMinVals){//min and max vals for entire dataset
		double[] minVals = new double[globalMinVals.length];
		double[] maxVals = new double[globalMaxVals.length];
		int K = centers.length;
		Random rand = new Random();
		int replaced = rand.nextInt(K);
		for(int d = 0; d < points[0].length; d++){
			centers[replaced][d] = rand.nextDouble() * (globalMaxVals[d] - globalMinVals[d]) + globalMinVals[d];
		}
		int[] pointToCluster = pointsToClusters(points,centers); //index of the center to which point i is closest to
		double[] distortion = new double[K];
		double[] UI = new double[K];//utility
		int[] lowUtilIdxs = new int[K];//contains index of cells with low utility
		int[] highUtilIdxs = new int[K];
		int numLowUtil = 0, numHighUtil = 0;
		double meanDistortion = 0, sumUtil = 0;//sum of utilities that are > 1
		for(int i = 0; i < points.length; i++){
			distortion[pointToCluster[i]] += dist(points[i],centers[pointToCluster[i]]);
		}
		for(int k = 0; k < K; k++){
			meanDistortion += distortion[k];
		}
		meanDistortion /= K;
		for(int k = 0; k < K; k++){
			UI[k] = distortion[k]/meanDistortion;
			if(UI[k] < 1){
				lowUtilIdxs[numLowUtil++] = k;
			}
			else{
				sumUtil += UI[k];
				highUtilIdxs[numHighUtil++] = k;
			}
		}
		for(int i = 0; i < numLowUtil; i++){
			pointToCluster = pointsToClusters(points,centers);
			int[] testingCounts = new int[K];
			for(int g = 0; g < K; g++){
				testingCounts[g] = 0;
			}
			for(int g = 0; g < points.length; g++){
				testingCounts[pointToCluster[g]]++;
			}
			//pick highUtilIdx
			double random = rand.nextDouble() * sumUtil;
			int highUtilIdxIdx = 0, highUtilIdx, lowUtilIdx = lowUtilIdxs[i];
			while(highUtilIdxIdx < numHighUtil){
				random -= UI[highUtilIdxs[highUtilIdxIdx]];
				if(random < 0)
					break;
				highUtilIdxIdx++;
			}
			highUtilIdx = highUtilIdxs[highUtilIdxIdx];
			//find closest center from the old location of the moved center
			int closestIdx = 0;
			double closestDist = 10000;
			for(int k = 0; k < K; k++){
				if(k != lowUtilIdx && dist(centers[lowUtilIdx],centers[k]) < closestDist){
					closestIdx = k;
					closestDist = dist(centers[lowUtilIdx],centers[k]);
				}
			}
			//Codeword shift
			for(int d = 0; d < points[0].length; d++){
				maxVals[d] = -10000;
				minVals[d] = 10000;
			}
			double[] origHighUtilCenter = new double[points[0].length];
			double[] origLowUtilCenter = new double[points[0].length];
			int numSubPts = 0;
			for(int j = 0; j < points.length; j++){
				if(pointToCluster[j] == highUtilIdx){
					numSubPts++;
					for(int d = 0; d < points[0].length; d++){
						if(maxVals[d] < points[j][d])
							maxVals[d] = points[j][d];
						else if(minVals[d] > points[j][d])
							minVals[d] = points[j][d];
					}
				}
			}
			for(int d = 0; d < points[0].length; d++){
				origHighUtilCenter[d] = centers[highUtilIdx][d];
				origLowUtilCenter[d] = centers[lowUtilIdx][d];
				centers[highUtilIdx][d] = 2/3.0*(maxVals[d] - minVals[d])+minVals[d];
				centers[lowUtilIdx][d] = 1/3.0*(maxVals[d] - minVals[d])+minVals[d];
			}
			//center correction
			double[][] changeCenters = new double[2][points[0].length];
			double[][] subPoints = new double[numSubPts][points[0].length];//points that belonged to the high utility center
			for(int d = 0; d < points[0].length; d++){
				changeCenters[0][d] = centers[highUtilIdx][d];
				changeCenters[1][d] = centers[lowUtilIdx][d];
			}
			int subPtIdx = 0;
			for(int j = 0; j < points.length; j++){
				if(pointToCluster[j] == highUtilIdx){
					for(int d = 0; d < points[0].length; d++){
						subPoints[subPtIdx][d] = points[j][d];
					}
					subPtIdx++;
				}
			}
			if(subPoints.length == 0){
				continue;
			}
			changeCenters = KMHelper(subPoints,changeCenters,0.001);
			for(int d = 0; d < points[0].length; d++){
				centers[highUtilIdx][d] = changeCenters[0][d];
				centers[lowUtilIdx][d] = changeCenters[1][d];
			}
			//reassign points
			int[] origpointToCluster = new int[points.length];
			for(int j = 0; j < points.length; j++){
				origpointToCluster[j] = pointToCluster[j];
				if(pointToCluster[j] == lowUtilIdx)
					pointToCluster[j] = closestIdx;
				if(pointToCluster[j] == highUtilIdx)
					if(dist(centers[highUtilIdx],points[j]) < dist(centers[lowUtilIdx],points[j]))
						pointToCluster[j] = highUtilIdx;
					else
						pointToCluster[j] = lowUtilIdx;
			}
			//recalculate other center
			for(int d = 0; d < points[0].length; d++){
				centers[closestIdx][d] = 0;
			}
			int numPts = 0;
			for(int j = 0; j < points.length; j++){
				origpointToCluster[j] = pointToCluster[j];
				if(pointToCluster[j] == closestIdx){
					numPts++;
					for(int d = 0; d < points[0].length; d++){
						centers[closestIdx][d] += points[j][d];
					}
				}
			}
			for(int d = 0; d < points[0].length; d++){
				if(numPts != 0)
					centers[closestIdx][d] /= numPts;
			}
			//validation
			double lowUtilDistort = 0, closestDistort = 0, highUtilDistort = 0, countCheck = 0;
			for(int j = 0; j < points.length; j++){
				if(pointToCluster[j] == closestIdx)
					closestDistort += dist(points[j],centers[pointToCluster[j]]);
				if(pointToCluster[j] == lowUtilIdx){
					lowUtilDistort += dist(points[j],centers[pointToCluster[j]]);
					countCheck += 1;
				}
				if(pointToCluster[j] == highUtilIdx){
					highUtilDistort += dist(points[j],centers[pointToCluster[j]]);
					countCheck += 1;
				}
			}
			if(closestDistort + lowUtilDistort + highUtilDistort > distortion[closestIdx] + distortion[lowUtilIdx] + distortion[highUtilIdx]){//get old solution back
				for(int j = 0; j < points.length; j++){
					pointToCluster[j] = origpointToCluster[j];
				}
				for(int d = 0; d < points[0].length; d++){
					centers[highUtilIdx][d] = origHighUtilCenter[d];
					centers[lowUtilIdx][d] = origLowUtilCenter[d];
				}
			}
			else{//update distortions: odds of reselecting the previous high util distortion should be reduced
				//remove highUtilIdx
				highUtilIdxs[highUtilIdxIdx] = highUtilIdxs[--numHighUtil];
				sumUtil -= UI[highUtilIdx];
				if(numHighUtil == 0)
					break;
			}
		}
		//calculate new performance value and see which solution to keep
		return KMPerf(points,centers);
	}

	public static double[][] SAGMDE(double[][] points, int K, boolean display){
		//for distortions
		double[] maxVals = new double[points[0].length];//max value for each attribute
		double[] minVals = new double[points[0].length];//min value for each attribute
		for(int d = 0; d < points[0].length; d++){
			maxVals[d] = -10000;
			minVals[d] = 10000;
			for(int i = 0; i < points.length; i++){
				if(maxVals[d] < points[i][d])
					maxVals[d] = points[i][d];
				if(minVals[d] > points[i][d])
					minVals[d] = points[i][d];
			}
		}
		double[][] centers = randomPartition(points,K);
		//variables for the algorithm
		Random rand = new Random();
		double currPerf = KMPerf(points,centers), newPerf = KMPerf(points,centers);
		//recommended values
		int MAXITER = (int)(2*points.length);
		double T = 0.001*1.5/*from 1.65*/, minT = 0.000001, alpha = /*0.98*/0.99, t2 = 3*2/*from 2.13*/, alpha2 = 0.9925, saveT = minT;///play around with t2 should be somewhere in the 1-3 range and the if statement

		//for display
		double[][] bestCenters = new double[K][points[0].length];
		double bestPerf = 100000000;
		while(T > minT){
			if(display)
				System.out.println(KMPerf(points,centers));
			for(int it = 0; it < MAXITER; it++){
				if(currPerf < bestPerf){
					bestPerf = currPerf;
					for(int k = 0; k < K; k++){
						for(int d = 0; d < points[0].length; d++){
							bestCenters[k][d] = centers[k][d];
						}
					}
					saveT = T;
				}
				//calculate new centers
				double[][] newCenters = new double[K][points[0].length];
				for(int c = 0; c < K; c++){
					for(int d = 0; d < points[0].length; d++){
						newCenters[c][d] = centers[c][d] + rand.nextGaussian()*Math.sqrt(T)/10;
					}
				}
				newPerf = KMPerf(points,newCenters);
				double randDouble = rand.nextDouble();
				if(newPerf - currPerf <= 0 || Math.pow(Math.E,(currPerf - newPerf)/T) > randDouble){//accept new solution
					centers = newCenters;
					currPerf = newPerf;
				}
				if(it % 20 == 0){
					double[][] oldCenters = new double[K][points[0].length];
					for(int k = 0; k < K; k++){
						for(int d = 0; d < points[0].length; d++){
							oldCenters[k][d] = centers[k][d];
						}
					}
					int replaced = rand.nextInt(K);
					for(int d = 0; d < points[0].length; d++){
						centers[replaced][d] = rand.nextDouble() * (maxVals[d] - minVals[d]) + minVals[d];
					}
					int[] pointToCluster = pointsToClusters(points,centers); //index of the center to which point i is closest to
					double[] distortion = new double[K];
					double[] UI = new double[K];//utility
					int[] lowUtilIdxs = new int[K];//contains index of cells with low utility
					int[] highUtilIdxs = new int[K];
					int numLowUtil = 0, numHighUtil = 0;
					double meanDistortion = 0, sumUtil = 0;//sum of utilities that are > 1
					for(int i = 0; i < points.length; i++){
						distortion[pointToCluster[i]] += dist(points[i],centers[pointToCluster[i]]);
					}
					for(int k = 0; k < K; k++){
						meanDistortion += distortion[k];
					}
					meanDistortion /= K;
					for(int k = 0; k < K; k++){
						UI[k] = distortion[k]/meanDistortion;
						if(UI[k] < 1)
						lowUtilIdxs[numLowUtil++] = k;
						else{
						sumUtil += UI[k];
						highUtilIdxs[numHighUtil++] = k;
						}
					}
					for(int i = 0; i < numLowUtil; i++){
						pointToCluster = pointsToClusters(points,centers);
						int[] testingCounts = new int[K];
						for(int g = 0; g < K; g++){
							testingCounts[g] = 0;
						}
						for(int g = 0; g < points.length; g++){
							testingCounts[pointToCluster[g]]++;
						}
						//pick highUtilIdx
						double random = rand.nextDouble() * sumUtil;
						int highUtilIdxIdx = 0, highUtilIdx, lowUtilIdx = lowUtilIdxs[i];
						while(highUtilIdxIdx < numHighUtil){
							random -= UI[highUtilIdxs[highUtilIdxIdx]];
							if(random < 0)
								break;
							highUtilIdxIdx++;
						}
						highUtilIdx = highUtilIdxs[highUtilIdxIdx];
						//find closest center from the old location of the moved center
						int closestIdx = 0;
						double closestDist = 10000;
						for(int k = 0; k < K; k++){
							if(k != lowUtilIdx && dist(centers[lowUtilIdx],centers[k]) < closestDist){
								closestIdx = k;
								closestDist = dist(centers[lowUtilIdx],centers[k]);
							}
						}
						//Codeword shift
						for(int d = 0; d < points[0].length; d++){
							maxVals[d] = -10000;
							minVals[d] = 10000;
						}
						double[] origHighUtilCenter = new double[points[0].length];
						double[] origLowUtilCenter = new double[points[0].length];
						int numSubPts = 0;
						for(int j = 0; j < points.length; j++){
							if(pointToCluster[j] == highUtilIdx){
								numSubPts++;
								for(int d = 0; d < points[0].length; d++){
									if(maxVals[d] < points[j][d])
										maxVals[d] = points[j][d];
									else if(minVals[d] > points[j][d])
										minVals[d] = points[j][d];
								}
							}
						}
						for(int d = 0; d < points[0].length; d++){
							origHighUtilCenter[d] = centers[highUtilIdx][d];
							origLowUtilCenter[d] = centers[lowUtilIdx][d];
							centers[highUtilIdx][d] = 2/3.0*(maxVals[d] - minVals[d])+minVals[d];
							centers[lowUtilIdx][d] = 1/3.0*(maxVals[d] - minVals[d])+minVals[d];
						}
						//center correction
						double[][] changeCenters = new double[2][points[0].length];
						double[][] subPoints = new double[numSubPts][points[0].length];//points that belonged to the high utility center
						for(int d = 0; d < points[0].length; d++){
							changeCenters[0][d] = centers[highUtilIdx][d];
							changeCenters[1][d] = centers[lowUtilIdx][d];
						}
						int subPtIdx = 0;
						for(int j = 0; j < points.length; j++){
							if(pointToCluster[j] == highUtilIdx){
								for(int d = 0; d < points[0].length; d++){
									subPoints[subPtIdx][d] = points[j][d];
								}
								subPtIdx++;
							}
						}
						if(subPoints.length == 0){
							continue;
						}
						changeCenters = KMHelper(subPoints,changeCenters,0.001);
						for(int d = 0; d < points[0].length; d++){
							centers[highUtilIdx][d] = changeCenters[0][d];
							centers[lowUtilIdx][d] = changeCenters[1][d];
						}
						//reassign points
						int[] origpointToCluster = new int[points.length];
						for(int j = 0; j < points.length; j++){
							origpointToCluster[j] = pointToCluster[j];
							if(pointToCluster[j] == lowUtilIdx)
								pointToCluster[j] = closestIdx;
							if(pointToCluster[j] == highUtilIdx){
								if(dist(centers[highUtilIdx],points[j]) < dist(centers[lowUtilIdx],points[j]))
									pointToCluster[j] = highUtilIdx;
								else
									pointToCluster[j] = lowUtilIdx;
							}
						}
						//recalculate other center
						for(int d = 0; d < points[0].length; d++){
							centers[closestIdx][d] = 0;
						}
						int numPts = 0;
						for(int j = 0; j < points.length; j++){
							origpointToCluster[j] = pointToCluster[j];
							if(pointToCluster[j] == closestIdx){
								numPts++;
								for(int d = 0; d < points[0].length; d++)
									centers[closestIdx][d] += points[j][d];
							}
						}
						for(int d = 0; d < points[0].length; d++){
							if(numPts != 0)
								centers[closestIdx][d] /= numPts;
						}
						//validation
						double lowUtilDistort = 0, closestDistort = 0, highUtilDistort = 0, countCheck = 0;
						for(int j = 0; j < points.length; j++){
							if(pointToCluster[j] == closestIdx)
								closestDistort += dist(points[j],centers[pointToCluster[j]]);
							if(pointToCluster[j] == lowUtilIdx){
								lowUtilDistort += dist(points[j],centers[pointToCluster[j]]);
								countCheck += 1;
							}
							if(pointToCluster[j] == highUtilIdx){
								highUtilDistort += dist(points[j],centers[pointToCluster[j]]);
								countCheck += 1;
							}
						}
						double pd1 = lowUtilDistort + highUtilDistort, pd2 = distortion[lowUtilIdx] + distortion[highUtilIdx];
						//System.out.println(pd1 + " " + pd2);
						if(closestDistort + lowUtilDistort + highUtilDistort > distortion[closestIdx] + distortion[lowUtilIdx] + distortion[highUtilIdx]){//get old solution back
						for(int j = 0; j < points.length; j++){
							pointToCluster[j] = origpointToCluster[j];
						}
						for(int d = 0; d < points[0].length; d++){
							centers[highUtilIdx][d] = origHighUtilCenter[d];
							centers[lowUtilIdx][d] = origLowUtilCenter[d];
						}
						}
						else{//update distortions: odds of reselecting the previous high util distortion should be reduced
						highUtilIdxs[highUtilIdxIdx] = highUtilIdxs[--numHighUtil];
						sumUtil -= UI[highUtilIdx];
						if(numHighUtil == 0)
							break;
						}
					}
					//calculate new performance value and see which solution to keep
					newPerf = KMPerf(points,centers);
					if(newPerf - currPerf <= 0 || Math.pow(Math.E,(currPerf - newPerf)/t2) > rand.nextDouble()){//accept solution
						currPerf = newPerf;
					}
					else{
						for(int k = 0; k < K; k++){
							for(int d = 0; d < points[0].length; d++){
								centers[k][d] = oldCenters[k][d];
							}
						}
					}
				}
			}
			t2 *= alpha2;
			T *= alpha;
		}

		centers = bestCenters;
		currPerf = bestPerf;
		T = saveT;
		while(T > minT){
			for(int it = 0; it < MAXITER; it++){
				//calculate new centers
				double[][] newCenters = new double[K][points[0].length];
				for(int c = 0; c < K; c++){
					for(int d = 0; d < points[0].length; d++){
						newCenters[c][d] = centers[c][d] + rand.nextGaussian()*Math.sqrt(T)/10;
					}
				}
				newPerf = KMPerf(points,newCenters);
				double randDouble = rand.nextDouble();
				if(newPerf - currPerf <= 0 || Math.pow(Math.E,(currPerf - newPerf)/T) > randDouble){//accept new solution
					centers = newCenters;
					currPerf = newPerf;
				}
			}
			T *= alpha;
		}
		return centers;
	}

	public static double[][] KMHelper(double[][] points, double[][] centers, double tune){
		int it = -1, maxIt = 300;
		while(++it < maxIt){
			int[] pointToCluster = pointsToClusters(points,centers); //index of the center to which point i is closest to
			int[] pointsPerCluster = new int[centers.length];
			for(int i = 0; i < points.length; i++){
				pointsPerCluster[pointToCluster[i]]++;
			}
			double[][] newCenters = new double[centers.length][points[0].length];
			//each center is equal to the sum of all points
			for(int i = 0; i < points.length; i++){
				for(int d = 0; d < points[0].length; d++){
						newCenters[pointToCluster[i]][d] += points[i][d];
				}
			}
			//divide each center by the number of points
			boolean done = true;
			for(int c = 0; c < centers.length; c++){
				for(int d = 0; d < points[0].length; d++){
					newCenters[c][d] /= pointsPerCluster[c];
					if(Math.abs(newCenters[c][d] - centers[c][d]) > tune)
						done = false;
				}
			}
			if(done)
				break;
			centers = newCenters;
		}
		return centers;
	}

	/*Assuming there are more points than clusters
		Returns a K by D array where K is the number of clusters and D is the number of dimensions
		Each row in the array represents a point in the D-dimensional space equivalent to the mean
		of the ith cluster.
		Clusters are formed such that each point belongs to the closest cluster.*/
	public static double[][] KHM(double[][] points, int K, double p){
		int maxIt = 10000000;
		double[][] centers = randomPartition(points,K);
		//Initialize variables for main algorithm
		double[][] dist = new double[points.length][K];//dist[a][b] contains the distance between point a and center b
		int[] minIdx = new int[points.length];//dist[a][minIdx[a]] contains the minimum distance between point a and a center
		//intermediaries for harmonic averages calculation
		double[] A = new double[points.length];
		double[][] Q = new double[points.length][K];
		double[] QQ = new double[K];
		double[][] P = new double[points.length][K];
		double[][] oldCenters = new double[K][points[0].length];
		for(int c = 0; c < K; c++){
			for(int d = 0; d < points[0].length; d++){
				oldCenters[c][d] = 0;
			}
		}
		//main loop
		int it = 0;
		for(; it < maxIt; it++){
			//check for convergence
			boolean done = false;
			for(int c = 0; c < K; c++){
				for(int d = 0; d < points[0].length; d++){
					if(Math.abs(centers[c][d] - oldCenters[c][d]) > 0.0001){
						c = K;
						break;
					}
				}
				if(c == K - 1)
					done = true;
			}
			if(done){
				break;
			}
			for(int c = 0; c < K; c++){
				for(int d = 0; d < points[0].length; d++){
						oldCenters[c][d] = centers[c][d];
				}
			}
			//reinitialize temporary arrays
			for(int i = 0; i < dist.length; i++){
				minIdx[i] = 0;
				A[i] = 0;
				if(i < K)
					QQ[i] = 0;
				for(int c = 0; c < K; c++){
					Q[i][c] = 0;
					P[i][c] = 0;
					dist[i][c] = 0;
				}
			}
			//update distances
			for(int i = 0; i < dist.length; i++){
				for(int c = 0; c < K; c++){
					dist[i][c] = dist(points[i],centers[c]);
					if(dist[i][minIdx[i]] > dist[i][c])
						minIdx[i] = c;
				}
			}
			//compute harmonic averages using formulas provided by B. Zhang. Generalized k-harmonic means – boosting in unsupervised learning. Technical Report HPL-2000-137, Hewlett-Packard Labs, 2000.
			for(int i = 0; i < dist.length; i++){
				for(int c = 0; c < K; c++){
					if(c != minIdx[i])//prevent division by 0
						A[i] += Math.pow(dist[i][minIdx[i]]/dist[i][c],p);
				}
			}
			for(int i = 0; i < dist.length; i++){
				for(int c = 0; c < K; c++){
					Q[i][c] = Math.pow(dist[i][minIdx[i]],p-2)/Math.pow(A[i]+1,2);
					if(c != minIdx[i])//prevent division by 0
						Q[i][c] *= Math.pow(dist[i][minIdx[i]]/dist[i][c],p+2);
				}
			}
			for(int c = 0; c < K; c++){
				for(int i = 0; i < dist.length; i++){
					QQ[c] += Q[i][c];
				}
			}
			for(int i = 0; i < dist.length; i++){
				for(int c = 0; c < K; c++){
					P[i][c] = Q[i][c] / QQ[c];
				}
			}
			for(int c = 0; c < K; c++){
				for(int d = 0; d < points[0].length; d++){
					centers[c][d] = 0;
				}
			}
			for(int c = 0; c < K; c++){
				for(int i = 0; i < dist.length; i++){
					for(int d = 0; d < points[0].length; d++){
						centers[c][d] += P[i][c] * points[i][d];
					}
				}
			}
		}
		return centers;
	}


	public static double[][] KHMHelper(double[][] points, double[][] centers, double p){
		int K = centers.length;
		//Initialize variables for main algorithm
		double[][] dist = new double[points.length][K];//dist[a][b] contains the distance between point a and center b
		int[] minIdx = new int[points.length];//dist[a][minIdx[a]] contains the minimum distance between point a and a center
		//intermediaries for harmonic averages calculation
		double[] A = new double[points.length];
		double[][] Q = new double[points.length][K];
		double[] QQ = new double[K];
		double[][] P = new double[points.length][K];
		//update distances
		for(int i = 0; i < dist.length; i++){
			for(int c = 0; c < K; c++){
				dist[i][c] = dist(points[i],centers[c]);
				if(dist[i][minIdx[i]] > dist[i][c])
					minIdx[i] = c;
			}
		}
		//compute harmonic averages using formulas provided by B. Zhang. Generalized k-harmonic means – boosting in unsupervised learning. Technical Report HPL-2000-137, Hewlett-Packard Labs, 2000.
		for(int i = 0; i < dist.length; i++){
			for(int c = 0; c < K; c++){
				if(c != minIdx[i])//prevent division by 0
					A[i] += Math.pow(dist[i][minIdx[i]]/dist[i][c],p);
			}
		}
		for(int i = 0; i < dist.length; i++){
			for(int c = 0; c < K; c++){
				Q[i][c] = Math.pow(dist[i][minIdx[i]],p-2)/Math.pow(A[i]+1,2);
				if(c != minIdx[i])//prevent division by 0
					Q[i][c] *= Math.pow(dist[i][minIdx[i]]/dist[i][c],p+2);
			}
		}
		for(int c = 0; c < K; c++){
			for(int i = 0; i < dist.length; i++){
				QQ[c] += Q[i][c];
			}
		}
		for(int i = 0; i < dist.length; i++){
			for(int c = 0; c < K; c++){
				P[i][c] = Q[i][c] / QQ[c];
			}
		}
		for(int c = 0; c < K; c++){
			for(int d = 0; d < points[0].length; d++){
				centers[c][d] = 0;
			}
		}
		for(int c = 0; c < K; c++){
			for(int i = 0; i < dist.length; i++){
				for(int d = 0; d < points[0].length; d++){
					centers[c][d] += P[i][c] * points[i][d];
				}
			}
		}
		return centers;
	}


	public static double[][] SAKHMC(double[][] points, int K, double p){
		Random rand = new Random();
		double[][] centers = randomPartition(points,K);
		double[][] bestCenters = new double[K][points[0].length];
		double T = 100, minT = 0.000001, t = 0.95;
		double perfVal = 100000000, oldPerfVal;
		int MAXITER = 100;
		double bestPerf = 100000000;
		double[] maxVals = new double[points[0].length];//max value for each attribute
		double[] minVals = new double[points[0].length];//min value for each attribute
		for(int d = 0; d < points[0].length; d++){
			maxVals[d] = -10000;
			minVals[d] = 10000;
			for(int i = 0; i < points.length; i++){
				if(maxVals[d] < points[i][d])
					maxVals[d] = points[i][d];
				if(minVals[d] > points[i][d])
					minVals[d] = points[i][d];
			}
		}
		while(T > minT){
			for(int it = 0; it < MAXITER; it++){
				if(perfVal < bestPerf){
					bestPerf = perfVal;
					for(int k = 0; k < K; k++){
						for(int d = 0; d < points[0].length; d++){
							bestCenters[k][d] = centers[k][d];
						}
					}
					System.out.println();
					System.out.println(KMPerf(points,centers));
					System.out.println(T);
				}
				oldPerfVal = KHMPerf(points,centers,3.5);
				double[][] oldCenters = new double[K][points[0].length];
				for(int k = 0; k < K; k++){
					for(int d = 0; d < points[0].length; d++){
						oldCenters[k][d] = centers[k][d];
						//centers[k][d] += rand.nextDouble()*0.08-0.04;
					}
				}
				int replaced = rand.nextInt(K);
				for(int d = 0; d < points[0].length; d++){
					centers[replaced][d] = rand.nextDouble() * (maxVals[d] - minVals[d]) + minVals[d];
				}
				int[] pointToCluster = pointsToClusters(points,centers); //index of the center to which point i is closest to
				double[] distortion = new double[K];
				double[] UI = new double[K];//utility
				int[] lowUtilIdxs = new int[K];//contains index of cells with low utility
				int[] highUtilIdxs = new int[K];
				int numLowUtil = 0, numHighUtil = 0;
				double meanDistortion = 0, sumUtil = 0;//sum of utilities that are > 1
				for(int i = 0; i < points.length; i++){
						distortion[pointToCluster[i]] += dist(points[i],centers[pointToCluster[i]]);
				}
				for(int k = 0; k < K; k++){
					meanDistortion += distortion[k];
				}
				meanDistortion /= K;
				for(int k = 0; k < K; k++){
					UI[k] = distortion[k]/meanDistortion;
					if(UI[k] < 1)
						lowUtilIdxs[numLowUtil++] = k;
					else{
						sumUtil += UI[k];
						highUtilIdxs[numHighUtil++] = k;
					}
				}
				for(int i = 0; i < numLowUtil; i++){
					double random = rand.nextDouble() * sumUtil;
					int highUtilIdx = 0, lowUtilIdx = lowUtilIdxs[i];
					while(highUtilIdx < numHighUtil){
						random -= UI[highUtilIdxs[highUtilIdx]];
						if(random < 0)
							break;
						highUtilIdx++;
					}
					highUtilIdx = highUtilIdxs[highUtilIdx];
					if(highUtilIdx < 1)
						continue;
					//Codeword shift
					for(int d = 0; d < points[0].length; d++){
						maxVals[d] = -10000;
						minVals[d] = 10000;
					}
					double[] origHighUtilCenter = new double[points[0].length];
					double[] origLowUtilCenter = new double[points[0].length];
					int numSubPts = 0;
					for(int j = 0; j < points.length; j++){
						if(pointToCluster[j] == highUtilIdx){
							numSubPts++;
							for(int d = 0; d < points[0].length; d++){
								if(maxVals[d] < points[j][d])
									maxVals[d] = points[j][d];
								else if(minVals[d] > points[j][d])
									minVals[d] = points[j][d];
							}
						}
					}
					for(int d = 0; d < points[0].length; d++){
						origHighUtilCenter[d] = centers[highUtilIdx][d];
						origLowUtilCenter[d] = centers[lowUtilIdx][d];
						centers[highUtilIdx][d] = 2/3.0*(maxVals[d] - minVals[d])+minVals[d];
						centers[lowUtilIdx][d] = 1/3.0*(maxVals[d] - minVals[d])+minVals[d];
					}
					//center correction
					double[][] changeCenters = new double[2][points[0].length];
					double[][] subPoints = new double[numSubPts][points[0].length];//points that belonged to the high utility center
					for(int d = 0; d < points[0].length; d++){
						changeCenters[0][d] = centers[highUtilIdx][d];
						changeCenters[1][d] = centers[lowUtilIdx][d];
					}
					int subPtIdx = 0;
					for(int j = 0; j < points.length; j++){
						if(pointToCluster[j] == highUtilIdx){
							for(int d = 0; d < points[0].length; d++){
								subPoints[subPtIdx][d] = points[j][d];
							}
							subPtIdx++;
						}
					}
					changeCenters = KMHelper(subPoints,changeCenters,0.001);
					for(int d = 0; d < points[0].length; d++){
						centers[highUtilIdx][d] = changeCenters[0][d];
						centers[lowUtilIdx][d] = changeCenters[1][d];
					}
					//find closest center from the old location of the moved center
					int closestIdx = 0;
					double closestDist = 10000;
					for(int k = 0; k < K; k++){
						if(k != lowUtilIdx && dist(centers[lowUtilIdx],centers[k]) < closestDist){
							closestIdx = k;
							closestDist = dist(centers[lowUtilIdx],centers[k]);
						}
					}
					//reassign points
					int[] origpointToCluster = new int[points.length];
					for(int j = 0; j < points.length; j++){
						origpointToCluster[j] = pointToCluster[j];
						if(pointToCluster[j] == lowUtilIdx)
							pointToCluster[j] = closestIdx;
						if(pointToCluster[j] == highUtilIdx){
							if(dist(centers[highUtilIdx],points[j]) < dist(centers[lowUtilIdx],points[j]))
								pointToCluster[j] = highUtilIdx;
							else
								pointToCluster[j] = highUtilIdx;
						}
					}
					//validation
					double lowUtilDistort = 0, closestDistort = 0, highUtilDistort = 0;
					for(int j = 0; j < points.length; j++){
						if(pointToCluster[j] == closestIdx)
							closestDistort += dist(points[j],centers[pointToCluster[j]]);
						if(pointToCluster[j] == lowUtilIdx)
							lowUtilDistort += dist(points[j],centers[pointToCluster[j]]);
						if(pointToCluster[j] == highUtilIdx)
							highUtilDistort += dist(points[j],centers[pointToCluster[j]]);
					}
					if(closestDistort + lowUtilDistort + highUtilDistort > distortion[closestIdx] + distortion[lowUtilIdx] + distortion[highUtilIdx]){//get old solution back
						for(int j = 0; j < points.length; j++){
							pointToCluster[j] = origpointToCluster[j];
						}
						for(int d = 0; d < points[0].length; d++){
							centers[highUtilIdx][d] = origHighUtilCenter[d];
							centers[lowUtilIdx][d] = origLowUtilCenter[d];
						}
					}
					else{//update distortions: odds of reselecting the previous high util distortion should be reduced
						meanDistortion *= K;
						double meanChange = closestDistort + lowUtilDistort + highUtilDistort - distortion[closestIdx] - distortion[lowUtilIdx] - distortion[highUtilIdx];
						meanDistortion += meanChange;
						meanDistortion /= K;
						sumUtil += distortion[highUtilIdx]/meanDistortion - UI[highUtilIdx];
						if(UI[closestIdx] > 1)
							sumUtil += distortion[closestIdx]/meanDistortion - UI[closestIdx];
						UI[closestIdx] = distortion[closestIdx]/meanDistortion;
						UI[highUtilIdx] = distortion[highUtilIdx]/meanDistortion;
						UI[lowUtilIdx] = distortion[lowUtilIdx]/meanDistortion;
					}
				}
				//calculate new performance value and see which solution to keep
				perfVal = KHMPerf(points,centers,3.5);
				if(perfVal > oldPerfVal && Math.pow(Math.E,(oldPerfVal-perfVal)/T) < rand.nextDouble()){//get old solution back
					for(int k = 0; k < K; k++){
						for(int d = 0; d < points[0].length; d++){
							centers[k][d] = oldCenters[k][d];
						}
					}
					perfVal = oldPerfVal;
				}
			}
			T *= t;
			//centers = KHMHelper(points,centers,3.5);
		}
		return bestCenters;
	}



	////////GENERATE DATA STRICTLY FOR TESTING
	//generates cluster based on the gaussian distribution given a number of points for the cluster, number of dimensions, and location of the center
	//pass in array to add points to, and specify index where to start adding points
	public static void genCluster(int numPts, int dim, double[] center, double[][] pts, int startIdx){
		Random rand = new Random();
		for(int i = 0; i < numPts; i++){
			for(int d = 0; d < dim; d++){
				pts[startIdx+i][d] = rand.nextGaussian() + center[d];
			}
		}
	}

	public static double[][] pointGen(int r, int c, int ptsPerCluster){//2d
		double[] center = new double[2];
		double[][] pts = new double[r*c*ptsPerCluster][2];
		center[0] = center[1] = 4*Math.sqrt(2);
		for(int i = 0; i < r; i++){
			for(int j = 0; j < c; j++){
				genCluster(ptsPerCluster,2,center,pts,i*ptsPerCluster*c+j*ptsPerCluster);
				center[1] += 4*Math.sqrt(2);
			}
			center[0] += 4*Math.sqrt(2);
			center[1] = 4*Math.sqrt(2);
		}
		return pts;
	}

	public static int[] pointsToClusters(double[][] points, double[][] centers){
		int[] ptsToCenter = new int[points.length];
		for(int i = 0; i < points.length; i++){
			double minDist = 100000;
			for(int c = 0; c < centers.length; c++){
				double currDist = dist(points[i],centers[c]);
				if(minDist > currDist){
					minDist = currDist;
					ptsToCenter[i] = c;
				}
			}
		}
		return ptsToCenter;
	}

	public static double KMPerf(double[][] points, double[][] centers){//measures how well algorithm performed based on KM performance function
		double perfVal = 0;
		int[] ptsToCenter = pointsToClusters(points,centers);
		for(int i = 0; i < points.length; i++){
			perfVal += Math.pow(dist(points[i],centers[ptsToCenter[i]]),2);
		}
		return perfVal;
	}

	public static double KMPerf(double[][] points, double[][] centers, int[] ptsToCenter){//measures how well algorithm performed based on KM performance function
		double perfVal = 0;
		for(int i = 0; i < points.length; i++)
			perfVal += Math.pow(dist(points[i],centers[ptsToCenter[i]]),2);
		return perfVal;
	}

	public static double KHMPerf(double[][] points, double[][] centers,double p){//measures how well algorithm performed based on KHM performance function
		double perfVal = 0;
		for(int i = 0; i < points.length; i++){
			double denominator = 0;
			for(int k = 0; k < centers.length; k++){
				denominator += 1/Math.pow(dist(points[i],centers[k]),p);
			}
			perfVal += centers.length/denominator;
		}
		return perfVal;
	}

	public static double NMI(String[] cls, double[][] points, double[][] centers){//returns an array of H
		double NMI = 0;
		int[] ptsToCenter = pointsToClusters(points,centers);
		Map<String,Integer> classInt=new HashMap<String,Integer>();//maps class names to integers
		int foundCategories = 0;
		for(int i = 0; i < cls.length; i++){
			if(!classInt.containsKey(cls[i])){
				classInt.put(cls[i],foundCategories++);
			}
		}
		int[][] countClassCluster = new int[foundCategories][centers.length];
		int[] numPerCluster = new int[centers.length];//number of points in each  cluster
		int[] numPerClass = new int[foundCategories];
		double[][] PYC = new double[foundCategories][centers.length];//for conditional entropy: [w][c] is probability of a data point being classified as c in cluster w
		double[] PY = new double[foundCategories];//probability of a class
		double[] PC = new double[centers.length];//probability of a cluster
		double HY = 0;//Entropy for class labels
		double HC = 0;//Entropy for class labels
		double[] HYC = new double[centers.length];//for each cluster
		double HYCTot = 0, IYC = 0;
		for(int i = 0; i < points.length; i++){
			countClassCluster[classInt.get(cls[i])][ptsToCenter[i]]++;
			numPerCluster[ptsToCenter[i]]++;
			numPerClass[classInt.get(cls[i])]++;
		}
		for(int w = 0; w < foundCategories; w++){
			PY[w] = numPerClass[w]/(double)points.length;
			HY -= PY[w]*log2(PY[w]);
		}
		for(int c = 0; c < centers.length; c++){
			PC[c] = numPerCluster[c]/(double)points.length;
			HC -= PC[c]*log2(PC[c]);
		}
		for(int w = 0; w < foundCategories; w++){
			for(int c = 0; c < centers.length; c++){
				if(numPerCluster[c] != 0)
					PYC[w][c] = countClassCluster[w][c]/(double)numPerCluster[c];
				else
					PYC[w][c] = 0;//occasionally may occur for conventional KM
			}
		}
		for(int c = 0; c < centers.length; c++){
			for(int w = 0; w < foundCategories; w++){
				HYC[c] -= PYC[w][c]*log2(PYC[w][c]);
			}
			HYC[c] *= PC[c];
		}
		for(int c = 0; c < centers.length; c++){
			HYCTot += HYC[c];
		}
		IYC = HY - HYCTot;
		NMI = (2*IYC)/(HY+HC);
		return NMI;
	}

	public static long nChooseTwo(int n){
		if(n < 2)
			return 0;
		return n*(n+1)/2;
	}

	public static int[] centersToLabels(String[] cls, Map<String,Integer> classInt, double[][] points, double[][] centers){//associates center with closest label; multiple centers may be associated with the same label; not used
		//find locations of truecenters
		int K = centers.length;
		double[][] trueCenters = new double[K][points[0].length];//centroid of points belonging to each class
		int[] centerToLabel = new int[K];
		int[] numPts = new int[K];//number of points in the given cluster
		for(int l = 0; l < K; l++){
			for(int d = 0; d < points[0].length; d++){
				trueCenters[l][d] = 0;
			}
			numPts[l] = 0;
			centerToLabel[l] = 0;
		}
		//assign points to clusters; centers will be equal to the sum of all points in the given cluster
		for(int i = 0; i < points.length; i++){
			numPts[classInt.get(cls[i])]++;
			for(int d = 0; d < points[0].length; d++){
				trueCenters[classInt.get(cls[i])][d] += points[i][d];
			}
		}
		for(int l = 0; l < K; l++){
			for(int d = 0; d < points[0].length; d++){
				trueCenters[l][d] /= numPts[l];
			}
		}
		//find nearest true center to each center found by the algorithm
		for(int c = 1; c < K; c++){
			for(int l = 0; l < K; l++){
				if(dist(centers[c],trueCenters[l]) < dist(centers[c],trueCenters[centerToLabel[c]]))
					centerToLabel[c] = l;
			}
		}
		return centerToLabel;
	}

	public static double RI(String[] cls, double[][] points, double[][] centers, boolean normalized){//http://faculty.washington.edu/kayee/pca/supp.pdf
		long TP = 0, TN = 0, FP = 0, FN = 0, sumIJ = 0, sumI = 0, sumJ = 0;
		int[] ptsToCenter = pointsToClusters(points,centers);
		Map<String,Integer> classInt=new HashMap<String,Integer>();//maps class names to integers
		int foundCategories = 0;
		for(int i = 0; i < cls.length; i++){
			if(!classInt.containsKey(cls[i])){
				classInt.put(cls[i],foundCategories++);
			}
		}
		int[][] countClassCluster = new int[foundCategories][centers.length];
		int[] numPerCluster = new int[centers.length];//number of points in each  cluster
		int[] numPerClass = new int[foundCategories];
		for(int i = 0; i < points.length; i++){
			countClassCluster[classInt.get(cls[i])][ptsToCenter[i]]++;
			numPerCluster[ptsToCenter[i]]++;
			numPerClass[classInt.get(cls[i])]++;
		}

		//(TP + TN)/(TP+TN+FP+FN)
		for(int l = 0; l < foundCategories; l++){
			for(int c = 0; c < centers.length; c++){
				System.out.print(countClassCluster[l][c] + " ");
			}
			System.out.println();
		}
		if(!normalized){
			double total = 0;
			for(int c = 0; c < centers.length; c++){
				int max = 0;
				for(int l = 0; l < foundCategories; l++){
					if(countClassCluster[l][c] > max)
						max = countClassCluster[l][c];
				}
				total += max;
			}
			return total/points.length;
		}

		//calculate intermediates
		for(int l = 0; l < foundCategories; l++){
			for(int c = 0; c < centers.length; c++){
			sumIJ += nChooseTwo(countClassCluster[l][c]);
			}
		}
		for(int c = 0; c < centers.length; c++){
			sumJ += nChooseTwo(numPerCluster[c]);
		}
		for(int l = 0; l < foundCategories; l++){
			sumI += nChooseTwo(numPerClass[l]);
		}

		TP = sumIJ;
		FN = sumI - sumIJ;
		FP = sumJ - sumIJ;
		TN = nChooseTwo(points.length) - TP - FN - FP;
		if(!normalized)
			return (TP + TN) / (double) (TP + FN + FP + TN);
		double temp = sumI/(double)nChooseTwo(points.length)*sumJ;
		return (sumIJ - temp)/(0.5*(sumI+sumJ)-temp);
	}

	public static void minMaxNormalize(double[][] points){
		double[] maxVals = new double[points[0].length];
		double[] minVals = new double[points[0].length];
		for(int d = 0; d < points[0].length; d++){
			minVals[d] = 100000;
		}
		for(int i = 0; i < points.length; i++){
			for(int d = 0; d < points[0].length; d++){
				if(maxVals[d] < points[i][d])
					maxVals[d] = points[i][d];
				if(minVals[d] > points[i][d])
					minVals[d] = points[i][d];
			}
		}
		for(int i = 0; i < points.length; i++){
			for(int d = 0; d < points[0].length; d++){
				points[i][d] = (points[i][d] - minVals[d])/(maxVals[d]-minVals[d]);
			}
		}
	}

	public static double[][] extractCenters(double[][] centers){
		Vector<Integer> condensedIdx = new Vector<Integer>();
		for(int c = 0; c < centers.length; c++){
			boolean duplicate = false;
			for(int i = 0; i < condensedIdx.size(); i++){
				if(dist(centers[c],centers[condensedIdx.get(i)]) < 0.01)
					duplicate = true;
			}
			if(duplicate)
				break;
			if(centers[c][0] > 1.9 && centers[c][0] < 2.1)
				break;
			condensedIdx.add(c);
		}
		double[][] retCenters = new double[condensedIdx.size()][centers[0].length];
		for(int i = 0; i < condensedIdx.size(); i++){
			for(int d = 0; d < centers[0].length; d++){
				retCenters[i][d] = centers[condensedIdx.get(i)][d];
			}
		}
		return retCenters;
	}

	public static void runner(int dataset, int algorithm, int t){
		//t is the number of iterations
		System.out.println();
		switch(algorithm){
			case 0:
				System.out.print("KM ");
				break;
			case 1:
				System.out.print("BH ");
				break;
			case 2:
				System.out.print("SAKM ");
				break;
			case 3:
				System.out.print("SAGM ");
				break;
			case 6:
				System.out.print("CPCL ");
				break;
			default:
				System.out.print("SAGMDE ");
		}
		//get extract and true cluster information from csv
		BufferedReader br = null;
		String line = "";
		String csvFile = "datasets/";
		int dim;
		int numPts;
		int numClusters;
		switch(dataset){
			case 0:
				csvFile += "glass.csv";
				System.out.print("glass");
				dim = 9;
				numPts = 214;
				numClusters = 6;
				break;
			case 1:
				csvFile += "ecoli.csv";
				System.out.print("ecoli");
				dim = 7;
				numPts = 336;
				numClusters = 8;
				break;
			case 2:
				csvFile += "iris.csv";
				System.out.print("iris");
				dim = 4;
				numPts = 150;
				numClusters = 3;
				break;
			case 3:
				csvFile += "wine.csv";
				System.out.print("wine");
				dim = 13;
				numPts = 178;
				numClusters = 3;
				break;
			case 4:
				csvFile += "breastCancer.csv";
				System.out.print("brestCancer");
				dim = 9;
				numPts = 683;
				numClusters = 2;
				break;
			case 5:
				csvFile += "yeast.csv";
				System.out.print("yeast");
				dim = 8;
				numPts = 1484;
				numClusters = 10;
				break;
			case 6:
				csvFile += "sat.csv";
				System.out.print("sat");
				dim = 36;
				numPts = 6435;
				numClusters = 6;
				break;
			case 7:
				csvFile += "synthetic(K=30).csv";
				System.out.print("syn(K=30)");
				dim = 2;
				numPts = 600;
				numClusters = 30;
				break;
			case 8:
				csvFile += "synthetic(K=15).csv";
				System.out.print("syn(K=15)");
				dim = 2;
				numPts = 600;
				numClusters = 15;
				break;
			case 9:
				csvFile += "synthetic(K=20).csv";
				System.out.print("syn(K=20)");
				dim = 2;
				numPts = 600;
				numClusters = 20;
				break;
			case 10:
				csvFile += "synthetic_control.csv";
				System.out.print("synCon");
				dim = 60;
				numPts = 600;
				numClusters = 6;
				break;
			case 11:
				csvFile += "mfeat.csv";
				System.out.print("mfeat");
				dim = 649;
				numPts = 2000;
				numClusters = 10;
				break;
			case 12:
				csvFile += "miceProtein.csv";
				System.out.print("mice");
				dim = 68;
				numPts = 1077;
				numClusters = 8;
				break;
			case 14:
				csvFile += "sonar.csv";
				System.out.print("sonar");
				dim = 60;
				numPts = 208;
				numClusters = 2;
				break;
			case 15:
				csvFile += "landCover.csv";
				System.out.print("land");
				dim = 147;
				numPts = 675;
				numClusters = 9;
				break;
			case 16:
				csvFile += "LSVT_voice_rehabilitation.csv";
				System.out.print("voice");
				dim = 309;
				numPts = 126;
				numClusters = 2;
				break;
			case 17:
				csvFile += "seeds_dataset.csv";
				System.out.print("seeds");
				dim = 7;
				numPts = 210;
				numClusters = 3;
				break;
			default:
				return;
		}
		System.out.println("");
		double[][] points = new double[numPts][dim];
		String[] trueCluster = new String[numPts];//the actual cluster the point should belong to based on the dataset
		int i = 0;
		try {
			br = new BufferedReader(new FileReader(csvFile));
			while ((line = br.readLine()) != null) {
			boolean incomplete = false;
			String[] l = line.split(",");
			for(int d = 0; d < dim; d++){
				if(l[d+1].equals("?"))
					incomplete = true;
			}
			if(!incomplete){
				if(i == 0 && dataset == 100)
					trueCluster[i] = l[0].substring(1);
				else
					trueCluster[i] = l[0];
				for(int d = 0; d < dim; d++){
					points[i][d] = Double.parseDouble(l[d+1]);
				}
				i++;
			}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
					try {
						br.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
			}
		}
		if(algorithm != 6)
			minMaxNormalize(points);
		/*for(i = 0; i < points.length; i++)
			System.out.println(points[i][0] + "," + points[i][1] + "," + points[i][2]);*/
		double[][][] allCenters = new double[t][numClusters*10][points[0].length];
		long[] times = new long[t];
		for(i = 0; i < t; i++){
			System.out.println("IT:" + i);
			long millis=System.currentTimeMillis();
			double[][] centers;
			switch(algorithm){
			case 0:
				centers = KM(points,numClusters);
				break;
			case 1:
				centers = BH(points,numClusters);
				break;
			case 2:
				centers = SAKM(points,numClusters);
				break;
			case 3:
				centers = SAGM(points,numClusters);
				break;
			case 4:
				centers = SAGMDE(points,numClusters,false);
				break;
			case 10:
				centers = KHM(points,numClusters,3.5);
				break;
			default:
				centers = SAGMDE(points,numClusters,false);
			}
			times[i]=System.currentTimeMillis() - millis;
			for(int j = 0; j < centers.length; j++){
				for(int d = 0; d < centers[0].length; d++){
					allCenters[i][j][d] = centers[j][d];
				}
			}
			allCenters[i][centers.length][0] = 2;
		}
		if(algorithm != 5){
			System.out.println("Results:");

			for(i = 0; i < t; i++){
				System.out.println(times[i]);
			}
			System.out.println();

			for(i = 0; i < t; i++){
				System.out.println(NMI(trueCluster,points,extractCenters(allCenters[i])));//System.out.println(KMPerf(points,allCenters[i]));
			}
			System.out.println();

			for(i = 0; i < t; i++){
				System.out.println(KMPerf(points,extractCenters(allCenters[i])));
			}
			System.out.println();
		}
		System.out.println("\n\n");
	}

	public static void main(String[] args){
		for(int dataset = 0; dataset < 16; dataset++){
			if(dataset != 2)
				continue;//just test 1 dataset for demo
			for(int alg = 0; alg <= 4; alg++){
				System.out.println("------------------------------------------------------------------------------");
				runner(dataset,alg,20);
			}
			System.out.println("********************************************************************************");
		}
	}
}