package cl.hierarchical.heuristic;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import cl.data.Cluster;
import cl.data.Instance;
import cl.data.Solution;
import cl.data.type.SolverType;
import cl.kmeans.heuristic.KMeansPlusPlusClustering;
import cl.util.MatrixOperations;
import cl.util.Quadruple;

/**
 * HMC algorithm based on Hierarchical Means Clustering by Maurizio Vichi, Carlo Cavicchia, and Patrick J. F. Groenen (2022)
 * 
 * The data structures have been adjusted to easily connect with the rest of the code
 */
public class HMC {
	private Instance instance;
	private int n; // Number of points
	private int K_start; // Number of clusters
	private int K_end;

	// HMC parameters
	private int Parsi = 2; // 1 if only parsimoneous hierarchical tree is required; 2 if a hierarchical complete (till K_end) tree is required 
	private int mlst = 100; // number of multistart for K-means 
	
	// Solution
	private Solution solution;
	private SolverType solverType;

	public HMC(Instance instance, int K_start, int K_end, SolverType solverType) {
		if(K_end<instance.getNumClusters() || K_end>instance.getNumPoints()) {
			throw new IllegalArgumentException("K_end has to be in the interval [K, numPoints]");
		}
		
		this.instance = instance;
		this.n = instance.getNumPoints();
		this.K_start = K_start;
		this.K_end = K_end;
		this.solverType = solverType;
	}

	public void solve() {
		// Keep track of different solutions
		Map<Integer, List<Cluster>> hierarchicalClusters = new LinkedHashMap<>();
		
		// Initialise variables
		RealMatrix X = MatrixUtils.createRealMatrix(instance.getCoordinates());		
		List<Integer> indKI = IntStream.range(0, K_start).boxed().collect(Collectors.toList());
		RealMatrix Par = MatrixUtils.createRealMatrix(n, n);
		RealMatrix incr = MatrixUtils.createRealMatrix(n-1, 1);
		double ofhmc1 = 0;
		double AAA = 1;

		// first apply k-means with K clusters (using multistart)
		KMeansPlusPlusClustering clustering = new KMeansPlusPlusClustering.Builder(K_start, instance.getCoordinates()).iterations(mlst).build();
		List<Cluster> clusters = clustering.getClusters();
		hierarchicalClusters.put(K_start-1, clusters);
		
		// Calculate Ultrametric UK
		double dwK = clustering.getWCSS();
		Quadruple<RealMatrix, RealMatrix, RealMatrix, RealMatrix> quadruple = clustering.getUltraMetricMatrix();
		RealMatrix UK = quadruple.first; // compute the membership matrix for the partition in K clusters  
		RealMatrix lb = quadruple.second;
		Par.setColumnMatrix(K_start-1, quadruple.third);
		RealMatrix dwp = quadruple.fourth;
		RealMatrix sUK = MatrixOperations.sumElementsVertically(UK);

		// Define matrices
		RealMatrix XmK = MatrixOperations.diagonal(MatrixOperations.toThePower(sUK, -1)).multiply(UK.transpose().multiply(X)); // compute the centroid matrix

		// compute dissimilarity matrix D associated to the K-cluster partition, necessary for the agglomerative part of the algorithm
		RealMatrix C = X.subtract(UK.multiply(XmK));
		double[][] d = new double[K_start][K_start];

		for(int k = 0; k < K_start; k++) {
			for(int i = k+1; i < K_start; i++) {
				Set<Integer> indK = new LinkedHashSet<>(indKI);
				indK.remove(k);
				indK.remove(i);

				RealMatrix uac = UK.getColumnMatrix(k).add(UK.getColumnMatrix(i)); // Merge clusters k and i
				RealMatrix UKm1 = MatrixUtils.createRealMatrix(n, K_start-1); // clusters of the partition in k-1 clusters in common with the partition in k clusters
				int count = 0;
				for(Integer j: indK) {
					UKm1.setColumnMatrix(count, UK.getColumnMatrix(j));
					count++;
				}

				// partition in k-1 clusters
				UKm1.setColumnMatrix(count, uac);
				RealMatrix sUKm1 = MatrixOperations.sumElementsVertically(UKm1);

				// compute centroid matrix of the partition in k clusters
				RealMatrix XmKm1 = MatrixOperations.diagonal(MatrixOperations.toThePower(sUKm1, -1)).multiply(UKm1.transpose().multiply(X));
				RealMatrix CC = X.subtract(UKm1.multiply(XmKm1));

				// compute the increase in the deviance merging cluster k and i
				d[k][i] = CC.transpose().multiply(MatrixOperations.diagonal(uac)).multiply(CC).getTrace() - C.transpose().multiply(MatrixOperations.diagonal(UK.getColumnMatrix(k))).multiply(C).getTrace() - C.transpose().multiply(MatrixOperations.diagonal(UK.getColumnMatrix(i))).multiply(C).getTrace();
				d[i][k] = d[k][i];
			}
		}
		RealMatrix D = MatrixUtils.createRealMatrix(d);

		// compute agglomerative part of the analysis on matrix D
		RealMatrix Aa  = MatrixUtils.createRealMatrix(K_start - 1, 1);
		RealMatrix Ba  = MatrixUtils.createRealMatrix(K_start - 1, 1);
		RealMatrix AA = MatrixUtils.createRealMatrix(K_start - 1, 1);
		RealMatrix BB = MatrixUtils.createRealMatrix(K_start - 1, 1);
		RealMatrix Ad  = MatrixUtils.createRealMatrix(n - K_start, 1);
		RealMatrix Bd  = MatrixUtils.createRealMatrix(n - K_start, 1);
		RealMatrix P  = MatrixUtils.createRealMatrix(K_start, 1);
		RealMatrix PP = MatrixUtils.createRealMatrix(K_start, 1);
		RealMatrix QQ = MatrixOperations.sumElementsVertically(UK);
		RealMatrix levfusa = MatrixUtils.createRealMatrix(K_start - 1, 1);

		// 'leader' of class to which object i belongs
		for(int i = 0; i < K_start; i++) {
			PP.setEntry(i, 0, i);
		}

		// start the agglomerative part of the algorithm
		// K = number of steps + 1 to which Ward's model has to be stopped
		for(int istep = 0; istep < K_start-1; istep++) {
			// compute the minimum distance
			int ic = -1;
			int jc = -1;
			double dmin = Double.POSITIVE_INFINITY;
			for(int i = 0; i < K_start-1; i++) {
				if(PP.getEntry(i, 0)==i) {
					for(int j = i+1; j < K_start; j++) {
						if(PP.getEntry(j, 0)==j) {
							if(D.getEntry(i, j) < dmin) {
								ic = i;
								jc = j;
								dmin = D.getEntry(i, j);
							}
						}
					}
				}
			}

			for(int j = 0; j < K_start; j++) {
				if(PP.getEntry(j, 0)==jc) {
					PP.setEntry(j, 0, ic);
				}
			}

			if(lb.getEntry(ic, 0)> lb.getEntry(jc, 0)) {
				double swp = lb.getEntry(ic, 0);
				lb.setEntry(ic, 0, lb.getEntry(jc, 0));
				lb.setEntry(jc, 0, swp);
			}
			Aa.setEntry(istep, 0, lb.getEntry(ic, 0));
			Ba.setEntry(istep, 0, lb.getEntry(jc, 0));
			AA.setEntry(istep, 0, ic);
			BB.setEntry(istep, 0, jc);
			levfusa.setEntry(istep, 0, 2*dmin);
			ofhmc1 = (K_start - istep - 1) * levfusa.getEntry(istep, 0)/2d + ofhmc1;
			incr.setEntry((K_start - istep - 1), 0, (K_start - istep - 1)*levfusa.getEntry(istep, 0)/2d);

			List<Integer> isw1 = findIndices(Par.getColumnMatrix(K_start-istep-1), lb.getEntry(ic, 0), false);
			List<Integer> isw2 = findIndices(Par.getColumnMatrix(K_start-istep-1), lb.getEntry(jc, 0), false);
			Par.setColumn(K_start-istep-2, Par.getColumn(K_start-istep-1));
			for(Integer i: isw1) {
				Par.setEntry(i, K_start-istep-2, lb.getEntry(ic, 0));
			}
			for(Integer i: isw2) {
				Par.setEntry(i, K_start-istep-2, lb.getEntry(ic, 0));
			}

			// Update hierarchical clustering
			List<Cluster> updatedClusters = new ArrayList<>();
			Map<Integer, boolean[]> temp = new LinkedHashMap<>();
			for(int m = 0; m < n; m++) {
				int val = (int) Par.getEntry(m, K_start-istep-2);
				if(!temp.containsKey(val)) {
					temp.put(val, new boolean[n]);
				}
				temp.get(val)[m] = true;
			}
			for(boolean[] e : temp.values()) {
				updatedClusters.add(new Cluster(e));	
			}
			hierarchicalClusters.put(K_start - istep - 1 -1, updatedClusters);
			
			// update distances 
			for(int i = 0; i < K_start; i++) {
				if(i!=ic && PP.getEntry(i, 0)==i) {
					double ds = (1./(QQ.getEntry(ic, 0) + QQ.getEntry(jc, 0) + QQ.getEntry(i, 0))) * ((QQ.getEntry(ic, 0) + QQ.getEntry(i, 0))*D.getEntry(ic,i) + (QQ.getEntry(i, 0) + QQ.getEntry(jc, 0))*D.getEntry(jc,i) - QQ.getEntry(i, 0)*dmin);
					D.setEntry(ic, i, ds);
					D.setEntry(i, ic, ds);
				}
			}
			QQ.setEntry(ic, 0, QQ.getEntry(ic, 0) + QQ.getEntry(jc, 0));
			List<Integer> ilbjc = findIndices(lb, lb.getEntry(jc, 0), false);
			lb.setEntry(ilbjc.get(0), 0, lb.getEntry(ic, 0));
		}

		// Write results
		RealMatrix Ult = MatrixUtils.createRealMatrix(K_start, K_start);
		for(int i=0; i < K_start ; i++) {
			P.setEntry(i, 0, i);
		}
		for(int k = 0; k < K_start-1; k++) {
			for(int i = (int) AA.getEntry(k, 0); i < K_start; i++) {
				if(P.getEntry(i, 0)==AA.getEntry(k, 0)) {
					for(int j = (int) BB.getEntry(k, 0); j < K_start; j++) {
						if(P.getEntry(j, 0)==BB.getEntry(k, 0)) {
							Ult.setEntry(i, j, levfusa.getEntry(k, 0));
							Ult.setEntry(j, i, levfusa.getEntry(k, 0));
						}
					}
				}
			}
			for(int i = (int) AA.getEntry(k, 0); i < K_start; i++) {
				if(P.getEntry(i, 0)==BB.getEntry(k, 0)) {
					P.setEntry(i, 0, AA.getEntry(k, 0));
				}
			}
		}
		// ultrametric associted to the tree
		Ult = (UK.multiply(Ult).multiply(UK.transpose())).add(UK.multiply(MatrixOperations.diagonal(dwp)).multiply(UK.transpose())).subtract(MatrixOperations.diagonal(MatrixOperations.selectDiagonalElements(UK.multiply(MatrixOperations.diagonal(dwp).multiply(UK.transpose())))));

		RealMatrix A = AA.transpose();
		RealMatrix B = BB.transpose();
		RealMatrix levfus = levfusa;

		if(Parsi!=1) { // if Parsi=1 only the parsimoneous dendrogram 
			// Extend dwp
			RealMatrix dwp_temp = MatrixUtils.createRealMatrix(K_end, 1);
			for(int j = 0; j < dwp.getRowDimension(); j++) {
				dwp_temp.setEntry(j, 0, dwp.getEntry(j, 0));
			}
			dwp = dwp_temp;
			
			// compute divisive part
			RealMatrix levfusd = MatrixUtils.createRealMatrix(n - K_start,1);
			int ist = 0;
			int Kn = K_start;
			RealMatrix ddscl = MatrixUtils.createRealMatrix(K_end, 1);
			RealMatrix Par2 = MatrixUtils.createRealMatrix(n, K_end);
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < K_end; j++) {
					Par2.setEntry(i, j, -1);
				}
			}

			RealMatrix dw2 = MatrixUtils.createRealMatrix(2, K_end);
			for(int k = 0; k < Kn; k++) {
				List<Integer> clustk = findIndices(UK.getColumnMatrix(k), 1, false); // identify units of cluster k 
				if(clustk.size()>1) { // if # cluster is 1 cannot be split
					double[][] data = new double[clustk.size()][instance.getNumDimensions()];
					int count = 0;
					for(Integer j: clustk) {
						data[count] = instance.getCoordinates()[j];
						count++;
					}
					KMeansPlusPlusClustering clustering2 = new KMeansPlusPlusClustering.Builder(2, data).iterations(mlst).build();
					Quadruple<RealMatrix, RealMatrix, RealMatrix, RealMatrix> quadruple2 = clustering2.getUltraMetricMatrix();
					RealMatrix P2 = quadruple2.first;
					dw2.setColumnMatrix(k, quadruple2.fourth);
					ddscl.setEntry(k, 0, (dwp.getEntry(k, 0) - clustering2.getWCSS()));
					int count2 = 0;
					for(Integer j : clustk) {
						if(P2.getEntry(count2, 0)==1) {
							Par2.setEntry(j, k, 0);
						}
						else {
							Par2.setEntry(j, k, 1);
						}
						count2++;
					}
				}
				else {
					Par2.setEntry(clustk.get(0), k, 0);
					ddscl.setEntry(k, 0, 0);
				}
			}
			
			int countK = K_start+1;
			for(int j = K_end; j > K_start; j--) {	
				double dwmax = Double.NEGATIVE_INFINITY;
				int pmax = -1;
				for(int m = 0; m < ddscl.getRowDimension(); m++) {
					if(ddscl.getEntry(m, 0)>dwmax) {
						dwmax = ddscl.getEntry(m, 0);
						pmax = m;
					}
				}

				if(dwmax>0) {
					List<Integer> isw3 = findIndices(Par2.getColumnMatrix(pmax), -1, true);
					double usw = Par2.getEntry(isw3.get(0), pmax);
					double sw2 = dw2.getEntry(1, pmax);
					if(usw>0) { // reorder cluster
						List<Integer> isw1 = findIndices(Par2.getColumnMatrix(pmax), 0, false);
						List<Integer> isw2 = findIndices(Par2.getColumnMatrix(pmax), 1, false);
						for(Integer m: isw1) {
							Par2.setEntry(m, pmax, 1);
						}
						for(Integer m: isw2) {
							Par2.setEntry(m, pmax, 0);
						}
						dw2.setEntry(1, pmax, dw2.getEntry(0, pmax));
						dw2.setEntry(0, pmax, sw2);
					}
					List<Integer> isw1 = findIndices(Par2.getColumnMatrix(pmax), 0, false);
					List<Integer> isw2 = findIndices(Par2.getColumnMatrix(pmax), 1, false);
					double mpcpc1 = Collections.min(isw1);
					double mpcpc2 = Collections.min(isw2);
					levfusd.setEntry(j-K_start-1, 0, 2*dwmax);
					ofhmc1 = ofhmc1 + (K_start + ist)*levfusd.getEntry(j - K_start-1, 0)/2d; // of as weighted sum of the incrementals
					incr.setEntry(K_start + ist -1, 0, (K_start + ist)*levfusd.getEntry(j - K_start-1, 0)/2d);
					ist++;
					Ad.setEntry(j - K_start-1, 0, mpcpc1);
	                Bd.setEntry(j - K_start-1, 0, mpcpc2);
	                Par.setColumn(K_start + ist-1, Par.getColumn(K_start + ist - 1));
	                for(int m: isw1) {
	                	Par.setEntry(m, K_start + ist-1, mpcpc1);
	                }
	                for(int m: isw2) {
	                	Par.setEntry(m, K_start + ist-1, mpcpc2);  
	                }
	                
	                // Update hierarchical information
					// Find out which cluster has been split
					Cluster splitCluster = null;
					List<Cluster> updateClusters = new ArrayList<>();
					for(Cluster cluster: clusters) {
						if(cluster.getPointId()[0]==isw1.get(0) || cluster.getPointId()[0]==isw2.get(0)) {
							splitCluster = cluster;
						}
						else {
							updateClusters.add(cluster);
						}
					}
					// Create two new clusters
					boolean[] cluster1 = new boolean[n];
					boolean[] cluster2 = new boolean[n];
					for(Integer m : isw1) {
						cluster1[m] = true;
					}
					for(Integer m : isw2) {
						cluster2[m] = true;
					}
					updateClusters.add(new Cluster(cluster1));
					updateClusters.add(new Cluster(cluster2));
					// Add clusters to hierarchy
					hierarchicalClusters.put(countK-1, updateClusters);
					clusters = updateClusters;
					countK++;
	                
	                dwK = dwK - dwmax; // with total variance 
	                dwp.setEntry(pmax, 0, dw2.getEntry(0, pmax)); // and individual variances
	                dwp.setEntry(Kn, 0, dw2.getEntry(1, pmax));
	              
	                int kk = 0;
	                RealMatrix P2o = Par2.getColumnMatrix(pmax);
	                for(int k: new int[]{pmax, Kn}) {
	                	List<Integer> icl1 = findIndices(P2o, kk, false);
	                	kk++;
	                	if(icl1.size()>1) {  // if # cluster is 1 cannot be split
	                		// split the cluster into two parts
	                		double[][] data = new double[icl1.size()][instance.getNumDimensions()];
	    					int count = 0;
	    					for(Integer m: icl1) {
	    						data[count] = instance.getCoordinates()[m];
	    						count++;
	    					}
	    	        		KMeansPlusPlusClustering clustering2 = new KMeansPlusPlusClustering.Builder(2, data).iterations(mlst).build();
	    	        		Quadruple<RealMatrix, RealMatrix, RealMatrix, RealMatrix> quadruple2 = clustering2.getUltraMetricMatrix();
	                		RealMatrix P2 = quadruple2.first;
	                		dw2.setColumnMatrix(k, quadruple2.fourth);
	                        for(int m = 0; m < n; m++) {
                            	Par2.setEntry(m, k, -1);
                            }
                            int count2 = 0;
        					for(Integer m : icl1) {
        						if(P2.getEntry(count2, 0)==1) {
        							Par2.setEntry(m, k, 0);
        						}
        						else {
        							Par2.setEntry(m, k, 1);
        						}
        						count2++;
        					}
                            ddscl.setEntry(k, 0, (dwp.getEntry(k, 0) - clustering2.getWCSS()));
	                	}
	                	else {
	                		for(int m = 0; m < n; m++) {
                            	Par2.setEntry(m, k, -1);
                            }
	                        Par2.setEntry(icl1.get(0), k, 0);
	                        ddscl.setEntry(k, 0, dwp.getEntry(k, 0));
	                	}
	                }
	                Kn++; // the partition is now in K+1 cluster 
				} 
				else {
					AAA = 2;
					throw new IllegalArgumentException("Warning! There are ties among dissimilarities");
				}
			}
				
			// Get hierarchical information
			hierarchicalClusters = hierarchicalClusters.entrySet().stream()
                    .sorted(Map.Entry.comparingByKey())
                    .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new));
			solution = new Solution(instance, hierarchicalClusters, instance.computeHierarchicalClusterCost(hierarchicalClusters));
			solution.setSolverType(solverType);
		}
	}

	public Solution getSolution() {
		return solution;
	}
	
	/**
	 * Input is a column vector
	 * @param mat
	 * @param value
	 * @return
	 */
	private List<Integer> findIndices(RealMatrix mat, double value, boolean largerThan) {
		if(mat.getColumnDimension()!=1) {
			throw new IllegalArgumentException("Input should be a column");
		}

		List<Integer> result = new ArrayList<>();
		for(int i = 0; i < mat.getRowDimension(); i++) {
			if(largerThan) {
				if(mat.getEntry(i, 0)>value) {
					result.add(i);
				}
			}
			else {
				if(mat.getEntry(i, 0)==value) {
					result.add(i);
				}
			}
		}
		return result;
	}
}
