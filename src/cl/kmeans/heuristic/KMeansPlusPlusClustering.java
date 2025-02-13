package cl.kmeans.heuristic;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import cl.data.Cluster;
import cl.data.GlobalParam;
import cl.util.Quadruple;

/**
 * KMeansPlusPlus, obtained from:
 * https://github.com/JasonAltschuler/KMeansPlusPlus/blob/master/src/KMeans.java
 * 
 *
 */
public class KMeansPlusPlusClustering {

	/***********************************************************************
	 * Data structures
	 **********************************************************************/
	// user-defined parameters
	private int k;                // number of centroids
	private double[][] points;    // n-dimensional data points. 

	// optional parameters
	private int iterations;       // number of times to repeat the clustering. Choose run with lowest WCSS
	private boolean pp;           // true --> KMeans++. false --> basic random sampling
	private double epsilon;       // stops running when improvement in error < epsilon
	private boolean useEpsilon;   // true  --> stop running when marginal improvement in WCSS < epsilon
	// false --> stop running when 0 improvement
	private boolean L1norm;       // true --> L1 norm to calculate distance; false --> L2 norm

	// calculated from dimension of points[][]
	private int m;                // number of data points   (# of pixels for PhenoRipper)  
	private int n;                // number of dimensions    (# of channels for PhenoRipper)

	// output
	private double[][] centroids; // position vectors of centroids                      dim(2): (k) by (number of channels)
	private int[] assignment;     // assigns each point to nearest centroid [0, k-1]    dim(1): (number of pixels)
	private double WCSS;          // within-cluster sum-of-squares. Cost function to minimize

	// timing information
	private long start;
	private long end;

	/***********************************************************************
	 * Constructors
	 **********************************************************************/

	/**
	 * Empty constructor is private to ensure that clients have to use the 
	 * Builder inner class to create a KMeans object.
	 */
	private KMeansPlusPlusClustering() {} 

	/**
	 * The proper way to construct a KMeans object: from an inner class object.
	 * @param builder See inner class named Builder
	 */
	private KMeansPlusPlusClustering(Builder builder) {
		// start timing
		start = System.currentTimeMillis();

		// use information from builder
		k = builder.k;
		points = builder.points;
		iterations = builder.iterations;
		pp = builder.pp;
		epsilon = builder.epsilon;
		useEpsilon = builder.useEpsilon;
		L1norm = builder.L1norm;

		// get dimensions to set last 2 fields
		m = points.length;
		n = points[0].length;

		// run KMeans++ clustering algorithm
		run();

		end = System.currentTimeMillis();
	}


	/**
	 * Builder class for constructing KMeans objects.
	 * 
	 * For descriptions of the fields in this (inner) class, see outer class.
	 */
	public static class Builder {
		// required
		private final int k;
		private final double[][] points;

		// optional (default values given)
		private int iterations     = 100;
		private boolean pp         = true;
		private double epsilon     = .001;
		private boolean useEpsilon = true;
		private boolean L1norm = false;

		/**
		 * Sets required parameters and checks that are a sufficient # of distinct
		 * points to run KMeans.
		 */
		public Builder(int k, double[][] originalPoints) {
			double[][] points = new double[originalPoints.length][originalPoints[0].length];
			for(int i = 0; i < originalPoints.length; i++) {
				for(int j = 0; j < originalPoints[0].length; j++) {
					points[i][j] = originalPoints[i][j];
				}
			}
			
			// check dimensions are valid
			if (k > points.length)
				throw new IllegalArgumentException("Required: # of points >= # of clusters");

			// check that there is a sufficient # of distinct points to run KMeans
			HashSet<double[]> hashSet = new HashSet<double[]>(k);
			int distinct = 0;

			for (int i = 0; i < points.length; i++) {
				if (!hashSet.contains(points[i])) {
					distinct++;
					if (distinct >= k)
						break;
					hashSet.add(points[i]);
				}
			}

			if (distinct < k)
				throw new IllegalArgumentException("Required: # of distinct points >= # of clusters");

			this.k = k;
			this.points = points;
		}


		/**
		 * Sets optional parameter. Default value is 100. 
		 */
		public Builder iterations(int iterations) {
			if (iterations < 1) 
				throw new IllegalArgumentException("Required: non-negative number of iterations. Ex: 50");
			this.iterations = iterations;
			return this;
		}

		/**
		 * Sets optional parameter. Default value is true.
		 */
		public Builder pp(boolean pp) {
			this.pp = pp;
			return this;
		}

		/**
		 * Sets optional parameter. Default value is .001.
		 */
		public Builder epsilon(double epsilon) {
			if (epsilon < 0.0)
				throw new IllegalArgumentException("Required: non-negative value of epsilon. Ex: .001"); 

			this.epsilon = epsilon;
			return this;
		}

		/**
		 * Sets optional parameter. Default value is true.
		 */
		public Builder useEpsilon(boolean useEpsilon) {
			this.useEpsilon = useEpsilon;
			return this;
		}

		/**
		 * Sets optional parameter. Default value is true
		 */
		public Builder useL1norm(boolean L1norm) {
			this.L1norm = L1norm;
			return this;
		}

		/**
		 * Build a KMeans object
		 */
		public KMeansPlusPlusClustering build() {
			if(points.length==2) {
				iterations = 1;
				pp = false;
				useEpsilon = false;
			}
			return new KMeansPlusPlusClustering(this);
		}
	}


	/***********************************************************************
	 * KMeans clustering algorithm
	 **********************************************************************/

	/** 
	 * Run KMeans algorithm
	 */
	private void run() {
		// for choosing the best run
		double bestWCSS = Double.POSITIVE_INFINITY;
		double[][] bestCentroids = new double[0][0];
		int[] bestAssignment = new int[0];

		// run multiple times and then choose the best run
		for (int n = 0; n < iterations; n++) {
			cluster();

			// store info if it was the best run so far
			if (WCSS < bestWCSS) {
				bestWCSS = WCSS;
				bestCentroids = centroids;
				bestAssignment = assignment;
			}
		}

		// keep info from best run
		WCSS = bestWCSS;
		centroids = bestCentroids;
		assignment = bestAssignment;
	}


	/**
	 * Perform KMeans clustering algorithm once.
	 */
	private void cluster() {
		int	stopIter = 100;
		// continue to re-cluster until marginal gains are small enough
		chooseInitialCentroids();
		assignmentStep();
		calcWCSS();
		int count = 0;
		double prevWCSS;
		do {  
			if(count > stopIter) {
				break;
			}
			count++;
			
			assignmentStep();   // assign points to the closest centroids

			updateStep();       // update centroids

			prevWCSS = WCSS;    // check if cost function meets stopping criteria
			calcWCSS();
		} while (!stop(prevWCSS));
	}


	/** 
	 * Assigns to each data point the nearest centroid.
	 */
	private void assignmentStep() {
		assignment = new int[m];

		double tempDist;
		double minValue;
		int minLocation;

		for (int i = 0; i < m; i++) {
			minLocation = 0;
			minValue = Double.POSITIVE_INFINITY;
			for (int j = 0; j < k; j++) {
				tempDist = distance(points[i], centroids[j]);
				if (tempDist < minValue) {
					minValue = tempDist;
					minLocation = j;
				}
			}

			assignment[i] = minLocation;
		}

	}


	/** 
	 * Updates the centroids.
	 */
	private void updateStep() {
		// reuse memory is faster than re-allocation
		for (int i = 0; i < k; i++)
			for (int j = 0; j < n; j++)
				centroids[i][j] = 0;

		int[] clustSize = new int[k];

		// sum points assigned to each cluster
		for (int i = 0; i < m; i++) {
			clustSize[assignment[i]]++;
			for (int j = 0; j < n; j++)
				centroids[assignment[i]][j] += points[i][j];
		}

		// store indices of empty clusters
		HashSet<Integer> emptyCentroids = new HashSet<Integer>();

		// divide to get averages -> centroids
		for (int i = 0; i < k; i++) {
			if (clustSize[i] == 0)
				emptyCentroids.add(i);

			else
				for (int j = 0; j < n; j++)
					centroids[i][j] /= clustSize[i];
		}

		// gracefully handle empty clusters by assigning to that centroid an unused data point
		if (emptyCentroids.size() != 0) {
			HashSet<double[]> nonemptyCentroids = new HashSet<double[]>(k - emptyCentroids.size());
			for (int i = 0; i < k; i++)
				if (!emptyCentroids.contains(i))
					nonemptyCentroids.add(centroids[i]);

			for (int i : emptyCentroids)
				while (true) {
					int rand = GlobalParam.RANDOM.nextInt(points.length);
					if (!nonemptyCentroids.contains(points[rand])) {
						nonemptyCentroids.add(points[rand]);
						centroids[i] = points[rand];
						break;
					}
				}

		}

	}


	/***********************************************************************
	 * Choose initial centroids
	 **********************************************************************/
	/**
	 * Uses either plusplus (KMeans++) or a basic randoms sample to choose initial centroids
	 */
	private void chooseInitialCentroids() {
		if (pp)
			plusplus();
		else
			basicSample();
	}

	/** 
	 * Choose two centroids
	 */
	private void basicSample() {
		if(m!=2) {
			throw new IllegalArgumentException("Only use this option when there are only 2 data points");
		}
		
		centroids = new double[k][n];
		double[][] copy = points;

		for (int i = 0; i < k; i++) {
			centroids[i] = copy[i];
		}
	}

	/** 
	 * Randomly chooses (without replacement) k data points as initial centroids using a
	 * weighted probability distribution (proportional to D(x)^2 where D(x) is the 
	 * distance from a data point to the nearest, already chosen centroid). 
	 */
	private void plusplus() {
		centroids = new double[k][n];       
		double[] distToClosestCentroid = new double[m];
		double[] weightedDistribution  = new double[m];  // cumulative sum of squared distances

		int choose = 0;

		for (int c = 0; c < k; c++) {

			// first centroid: choose any data point
			if (c == 0)
				choose = GlobalParam.RANDOM.nextInt(m);

			// after first centroid, use a weighted distribution
			else {

				// check if the most recently added centroid is closer to any of the points than previously added ones
				for (int p = 0; p < m; p++) {
					// gives chosen points 0 probability of being chosen again -> sampling without replacement
					double tempDistance = Distance.L2(points[p], centroids[c - 1]); // need L2 norm here, not L1

					// base case: if we have only chosen one centroid so far, nothing to compare to
					if (c == 1)
						distToClosestCentroid[p] = tempDistance;

					else { // c != 1 
						if (tempDistance < distToClosestCentroid[p])
							distToClosestCentroid[p] = tempDistance;
					}

					// no need to square because the distance is the square of the euclidean dist
					if (p == 0)
						weightedDistribution[0] = distToClosestCentroid[0];
					else weightedDistribution[p] = weightedDistribution[p-1] + distToClosestCentroid[p];

				}

				// choose the next centroid
				double rand = GlobalParam.RANDOM.nextDouble();
				for (int j = m - 1; j > 0; j--) {
					// starts at the largest bin. EDIT: not actually the largest
					if (rand > weightedDistribution[j - 1] / weightedDistribution[m - 1]) { 
						choose = j; // one bigger than the one above
						break;
					}
					else // Because of invalid dimension errors, we can't make the forloop go to j2 > -1 when we have (j2-1) in the loop.
						choose = 0;
				}
			}  

			// store the chosen centroid
			for (int i = 0; i < n; i++)
				centroids[c][i] = points[choose][i];
		}   
	}


	/***********************************************************************
	 * Cutoff to stop clustering
	 **********************************************************************/    

	/**
	 * Calculates whether to stop the run
	 * @param prevWCSS error from previous step in the run
	 * @return
	 */
	private boolean stop(double prevWCSS) {
		if (useEpsilon)
			return epsilonTest(prevWCSS);
		else
			return prevWCSS == WCSS; 
		// could this take infinite amount of time? double compare...
		// I think not because WCSS is calc in same way as prevWCSS (if data structs don't change)
	}

	/**
	 * Signals to stop running KMeans when the marginal improvement in WCSS
	 * from the last step is small.
	 * @param prevWCSS error from previous step in the run
	 * @return
	 */
	private boolean epsilonTest(double prevWCSS) {
		return epsilon > 1 - (WCSS / prevWCSS);
	}

	/***********************************************************************
	 * Utility functions
	 **********************************************************************/
	/**
	 * Calculates distance between two n-dimensional points.
	 * @param x
	 * @param y
	 * @return
	 */
	private double distance(double[] x, double[] y) {
		return L1norm ? Distance.L1(x, y) : Distance.L2(x, y);
	}

	private static class Distance {

		/**
		 * L1 norm: distance(X,Y) = sum_i=1:n[|x_i - y_i|]
		 * <P> Minkowski distance of order 1.
		 * @param x
		 * @param y
		 * @return
		 */
		public static double L1(double[] x, double[] y) {
			if (x.length != y.length) throw new IllegalArgumentException("dimension error");
			double dist = 0;
			for (int i = 0; i < x.length; i++) 
				dist += Math.abs(x[i] - y[i]);
			return dist;
		}

		/**
		 * L2 norm: distance(X,Y) = sqrt(sum_i=1:n[(x_i-y_i)^2])
		 * <P> Euclidean distance, or Minkowski distance of order 2.
		 * @param x
		 * @param y
		 * @return
		 */
		public static double L2(double[] x, double[] y) {
			if (x.length != y.length) throw new IllegalArgumentException("dimension error");
			double dist = 0;
			for (int i = 0; i < x.length; i++)
				dist += Math.abs((x[i] - y[i]) * (x[i] - y[i]));
			return dist;
		}
	}

	/** 
	 * Calculates WCSS (Within-Cluster-Sum-of-Squares), a measure of the clustering's error.
	 */
	private void calcWCSS() {
		double WCSS = 0;
		int assignedClust;

		for (int i = 0; i < m; i++) {
			assignedClust = assignment[i];
			WCSS += distance(points[i], centroids[assignedClust]);
		}     

		this.WCSS = WCSS;
	}

	/***********************************************************************
	 * Accessors
	 ***********************************************************************/
	public int[] getAssignment() {
		return assignment;
	}
	
	/**
	 * Get sum-of-squares for each centroid
	 * @return
	 */
	public double[] getClusterWCSS() {
		double[] result = new double[k];
		for (int i = 0; i < m; i++) {
			int assignedClust = assignment[i];
			result[assignedClust] += distance(points[i], centroids[assignedClust]);
		}   
		return result;
	}

	public double[][] getCentroids() {
		return centroids;
	}

	public double getWCSS() {
		return WCSS;
	}

	public String getTiming() {
		return "KMeans++ took: " + (double) (end - start) / 1000.0 + " seconds";
	}

	public List<Cluster> getClusters() {
		List<Cluster> result = new ArrayList<>();

		for(int i = 0; i < k; i++) {
			boolean[] cluster = new boolean[m];
			for(int j = 0; j < m; j++) {
				if(assignment[j]==i) {
					cluster[j] = true;
				}
			}
			result.add(new Cluster(cluster));
		}
		return result;
	}

	/**
	 * Returns UK , Par, and lb, dwp
	 * @return
	 */
	public Quadruple<RealMatrix, RealMatrix, RealMatrix, RealMatrix> getUltraMetricMatrix() {
		Map<Integer, List<Integer>> mapCluster = new LinkedHashMap<>();
		for(int i = 0; i < k; i++) {
			mapCluster.put(i, new ArrayList<>());
		}
		
		for(int i = 0; i < m; i++) {
			mapCluster.get(assignment[i]).add(i);
		}
		
		// Find out which cluster should go first
		Map<Integer, Integer> lowest = new LinkedHashMap<>();
		for(int i = 0; i < k; i++) {
			int lowestValue = Integer.MAX_VALUE;
			for(Integer j: mapCluster.get(i)) {
				if(j < lowestValue) {
					lowestValue = j;
				}
			}
			lowest.put(i, lowestValue);
		}
		
		// Sort based on values
		Map<Integer,Integer> sorted =
				lowest.entrySet().stream()
			       .sorted(Map.Entry.comparingByValue())
			       .collect(Collectors.toMap(
			          Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new));
		
		double[][] Uk = new double[m][k];
		double[][] Par = new double[m][1];
		double[][] lb = new double[k][1];
		double[] originalCWSS = getClusterWCSS();
		double[][] dwp = new double[k][1];
		int count = 0;
		for(Entry<Integer, Integer> entry: sorted.entrySet()) {
			int lowestValue = entry.getValue();
			lb[count][0] = lowestValue;
			dwp[count][0] = originalCWSS[entry.getKey()];
			for(Integer j : mapCluster.get(entry.getKey())) {
				Par[j][0] = lowestValue;
				Uk[j][count] = 1;
			}
			count++;
		}
		return new Quadruple<>(MatrixUtils.createRealMatrix(Uk), MatrixUtils.createRealMatrix(lb), MatrixUtils.createRealMatrix(Par), MatrixUtils.createRealMatrix(dwp));
	}
}