package cl.hierarchical.model.extended;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.jgrapht.Graph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cl.data.Cluster;
import cl.data.GlobalParam;
import cl.data.Instance;
import cl.data.ObjectiveWeightFunction;
import cl.hierarchical.model.hp.CliqueFinder_Knapsack;
import cl.util.Pair;
import cl.util.Triple;
import ilog.concert.IloException;


public class PricingProblemColgenMedoid implements PricingProblem {
	private static Logger log = LoggerFactory.getLogger(PricingProblemColgenMedoid.class);
	
	private Instance instance;
	private DualInformation duals;
	private BranchInformation bi;
	private Set<Integer> exclusionSet = new LinkedHashSet<>();

	private Map<Cluster, Double> resultPotential; //Cluster and reducedCost
	private Triple<List<Integer>, Double, Double> bestClique;
	private Map<Integer, Double> bestCliquePerLevel; // Only store the best objective
	private double bestRC;

	// Use for enumeration
	private int countEnumerate;
	
	// Used for branching
	private BranchInformationExtra bix;
	
	private int columnLimit = GlobalParam.PP_NUM_COLUMNS;

	public PricingProblemColgenMedoid(Instance instance) {
		this.instance = instance;
	}
	
	public void setColumnLimit(int columnLimit) {
		this.columnLimit = columnLimit;
	}
	
	public double getBestRC() {
		return bestRC;
	}

	public void setDuals(DualInformation duals) {
		this.duals = duals;
	}

	public void setBranchInformation(BranchInformation bi) {
		this.bi = bi;
	}

	public void setExclusionSet(Set<Integer> exclusionSet) {
		this.exclusionSet = exclusionSet;
	}

	public void setBix(BranchInformationExtra bix) {
		this.bix = bix;
	}

	public Map<Integer, Set<Cluster>> generateColumns() {
		return generateColumns(false, null);
	}

	/**
	 * Reduced cost (y_t) = c_t - \sum_k a_kt lambda_k - z
	 * Find minimum RC
	 * @return
	 * @throws IloException 
	 */
	public Map<Integer, Set<Cluster>> generateColumns(boolean enumerate, Map<Integer, Double> bounds) {
		if(enumerate) {
			throw new IllegalArgumentException("Enumerate has not been implemented yet for Medoids");
		}
		
		Map<Integer, Set<Cluster>> result = new LinkedHashMap<>();
		bestCliquePerLevel = new LinkedHashMap<>();
		bestRC = Double.POSITIVE_INFINITY;

		for(int i = 1; i < instance.getStartEndLevelList().size(); i++) { // We skip the first cluster, since it contains all elements
			if(instance.getStartEndLevelList().get(i).second==0) { // We skip clusters that end at level 0
				continue;
			}
			
			if(enumerate && bounds!=null) {
				if(!bounds.containsKey(i)) {
					continue;
				}
				if(bounds.get(i) > GlobalParam.NEGATIVE_RC_THRESHOLD) {
					continue;
				}
				GlobalParam.ENUMERATE_GRAPH_THRESHOLD = - bounds.get(i);	
				GlobalParam.ENUMERATE_RC_COST_THRESHOLD = GlobalParam.NEGATIVE_RC_THRESHOLD;
				System.out.println("Start enumeration "+ i + " "+ GlobalParam.ENUMERATE_GRAPH_THRESHOLD + " "+ GlobalParam.ENUMERATE_RC_COST_THRESHOLD);
			}
			
			// Reset BIX
			bix = null;
			try {
				generatePotentialColumns(i, enumerate);
			} catch (IloException e) {
				e.printStackTrace();
			}

			if(enumerate) {
				System.out.println(i + " clusters found " + resultPotential.size());
				int limit = resultPotential.size();
				if(resultPotential.size()<GlobalParam.ENUMERATE_PP_NUM_COLUMNS_ADD_LIMIT || GlobalParam.ENUMERATE_ADD_ALL_COLUMNS) {
					// Add everything
				}
				else if(resultPotential.size()>GlobalParam.ENUMERATE_PP_NUM_COLUMNS) {
					limit = Math.min(limit, GlobalParam.ENUMERATE_PP_NUM_COLUMNS);
				}
				
				Map<Cluster,Double> sorted =
						resultPotential.entrySet().stream() // TODO: Only keep negative RC smaller than epsilon
						.sorted(Map.Entry.comparingByValue())
						.limit(limit)
						.collect(Collectors.toMap(
								Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new)); //Ascending order
				result.put(i, sorted.keySet());
			}
			else {
				Map<Cluster,Double> sorted =
						resultPotential.entrySet().stream() // TODO: Only keep negative RC smaller than epsilon
						.sorted(Map.Entry.comparingByValue())
						.limit(columnLimit)
						.collect(Collectors.toMap(
								Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new)); //Ascending order
				result.put(i, sorted.keySet());
				
				for(Double obj: sorted.values()) {
					if(obj < bestRC) {
						bestRC = obj;
					}
				}
			}
		}
		log.info("Lowest RC {}", bestRC);
		return result;
	}

	public Pair<Triple<List<Integer>, Double, Double>, Map<Cluster, Double>> generateColumns(int level, boolean enumerate) {
		try {
			generatePotentialColumns(level, enumerate);
		} catch (IloException e) {
			e.printStackTrace();
		}

		Map<Cluster,Double> sorted =
				resultPotential.entrySet().stream() // TODO: Only keep negative RC smaller than epsilon
				.sorted(Map.Entry.comparingByValue())
				.limit(columnLimit)
				.collect(Collectors.toMap(
						Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new)); //Ascending order

		return new Pair<>(bestClique, sorted);
	}

	public List<List<Integer>> createBranchMapping(int level) {
		// Check if there is branching
		// Find out which points are enforced together
		// See branching in Aloise (2012)
		List<List<Integer>> branchMapping = new ArrayList<>(); // Key is the new ID
		List<Integer> redundant = new ArrayList<>();
		for(int i = 0; i < instance.getNumPoints(); i++) {
			for(int j = i+1 ; j < instance.getNumPoints(); j++) {
				if(bi.getEnforced(level).contains(new Pair<>(i, j))) {
					boolean found = false;
					for(List<Integer> oldIDs: branchMapping) {
						if(oldIDs.contains(i) || oldIDs.contains(j)) {
							if(!oldIDs.contains(j)) {
								oldIDs.add(j);
								redundant.add(j);
							}
							else if(!oldIDs.contains(i)) {
								oldIDs.add(i);
								redundant.add(i);
							}
							found = true;
							break;
						}
					}
					if(!found) {
						List<Integer> oldIDs = new ArrayList<>();
						oldIDs.add(i);
						oldIDs.add(j);
						redundant.add(i);
						redundant.add(j);
						branchMapping.add(oldIDs);
					}
				}
			}
		}

		// Check if we should merge any of the oldIDs
		while(true) {
			boolean overlapFound = false;
			outerloop: for(int i = 0; i < branchMapping.size(); i++) {
				for(int j = i+1; j < branchMapping.size(); j++) {
					Set<Integer> oldIDs1 = new LinkedHashSet<>(branchMapping.get(i));
					Set<Integer> oldIDs2 = new LinkedHashSet<>(branchMapping.get(j));
					
					oldIDs1.retainAll(oldIDs2); // calculate intersection
					if(!oldIDs1.isEmpty()) {
						// intersection is non empty
						oldIDs1 = new LinkedHashSet<>(branchMapping.get(i));
						oldIDs2 = new LinkedHashSet<>(branchMapping.get(j));
						
						oldIDs1.addAll(oldIDs2);
						branchMapping.remove(j);
						branchMapping.remove(i);
						branchMapping.add(new ArrayList<>(oldIDs1));
						overlapFound = true;
						break outerloop;
					}
				}
			}
			
			if(!overlapFound) {
				break;
			}
		}
		
		// Sort oldIDs
		for(List<Integer> oldIDs: branchMapping) {
			Collections.sort(oldIDs);
		}

		// New points: first the enforced points, then the old points
		for(int i = 0; i < instance.getNumPoints(); i++) {
			if(redundant.contains(i)) {
				continue;
			}

			List<Integer> oldIDs = new ArrayList<>();
			oldIDs.add(i);
			branchMapping.add(oldIDs);
		}
		return branchMapping;
	}
	
	/**
	 * Return false if there are inconsistencies
	 * @param level
	 * @return
	 */
	private boolean createReducedInstance(int level) {
		bix = new BranchInformationExtra();
		if(bi==null) { // Only in the rootnode
			bix.setBranch(false);
			bix.setReducedInstance(instance);
		}
		else {
			bix.setBranch(true);

			List<List<Integer>> branchMapping = createBranchMapping(level);
			// Check if there are no inconsistencies
			for(List<Integer> oldIDs: branchMapping) {
				for(int m = 0; m < oldIDs.size(); m++) {
					for(int n = m+1; n < oldIDs.size(); n++) {
						Pair<Integer, Integer> pair = new Pair<>(oldIDs.get(m), oldIDs.get(n));
						if(bi.getForbidden(level).contains(pair)) {
							return false;
						}
					}
				}
			}
			int numPoints = branchMapping.size();

			// Precalculate weights
			List<Integer> branchMappingWeights = new ArrayList<>();
			for(List<Integer> oldIDs: branchMapping) {
				branchMappingWeights.add(oldIDs.size());
			}
			
			double dist = 0;
			Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
			dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
			
			// Fill the distance matrix
			double[][] regularDistanceMatrix = new double[numPoints][numPoints];
			double[][] squaredEuclideanDistanceMatrix = new double[numPoints][numPoints];
			for(int i = 0; i < numPoints; i++) {
				for(int j = i; j < numPoints; j++) {
					if(i==j) {
						regularDistanceMatrix[i][i] = 0;
						squaredEuclideanDistanceMatrix[i][i] = 0;	
						List<Integer> oldIDs = branchMapping.get(i);
						for(int m = 0; m < oldIDs.size(); m++) {
							for(int n = m+1; n < oldIDs.size(); n++) {
								squaredEuclideanDistanceMatrix[i][i] += (dist * instance.getSquaredDistanceMatrix()[oldIDs.get(m)][oldIDs.get(n)]);
							}
						}
						regularDistanceMatrix[i][i] = Math.sqrt(squaredEuclideanDistanceMatrix[i][i]);
					}
					else {
						regularDistanceMatrix[i][j] = 0;
						squaredEuclideanDistanceMatrix[i][j] = 0;	
						List<Integer> oldIDs1 = branchMapping.get(i);
						List<Integer> oldIDs2 = branchMapping.get(j);
						for(int m = 0; m < oldIDs1.size(); m++) {
							for(int n = 0; n < oldIDs2.size(); n++) {
								squaredEuclideanDistanceMatrix[i][j] += (dist * instance.getSquaredDistanceMatrix()[oldIDs1.get(m)][oldIDs2.get(n)]);
								Pair<Integer, Integer> pair2;
								if(oldIDs1.get(m)<oldIDs2.get(n)) {
									pair2 = new Pair<>(oldIDs1.get(m), oldIDs2.get(n));
								}
								else {
									pair2 = new Pair<>(oldIDs2.get(n), oldIDs1.get(m));
								}
								if(bi.getForbidden(level).contains(pair2)) {
									regularDistanceMatrix[i][j] = GlobalParam.BIGM8;
									squaredEuclideanDistanceMatrix[i][j] = GlobalParam.BIGM8;
									regularDistanceMatrix[j][i] = GlobalParam.BIGM8;
									squaredEuclideanDistanceMatrix[j][i] = GlobalParam.BIGM8;
								}
							}
						}
						regularDistanceMatrix[i][j] = Math.sqrt(squaredEuclideanDistanceMatrix[i][j]);
					}
				}
			}
			Instance reducedInstance = new Instance(instance.getInstanceName()+"_branch", numPoints, instance.getNumDimensions(), instance.getNumClusters(), squaredEuclideanDistanceMatrix, regularDistanceMatrix);

			// Construct reducedDuals
			double[] dualsPoint = new double[reducedInstance.getNumPoints()];
			double[] dualsSqrtPoint = new double[reducedInstance.getNumPoints()];
			for(int j = 0; j < reducedInstance.getNumPoints(); j++) {
				dualsPoint[j] = 0;
				for(Integer i: branchMapping.get(j)) {
					dualsPoint[j] += duals.getDualsPoint()[level][i];
				}
				dualsSqrtPoint[j] = Math.sqrt(dualsPoint[j]);
			}

			bix.setReducedInstance(reducedInstance);
			bix.setReducedDuals(new DualInformation(dualsPoint, dualsSqrtPoint));
			bix.setBranchMapping(branchMapping);
			bix.setBranchMappingWeights(branchMappingWeights);
		}
		return true;
	}


	/**
	 * Algorithm 2 from Aloise (2012)
	 * @throws IloException 
	 */
	private void generatePotentialColumns(int level, boolean enumerate) throws IloException {
		resultPotential = new LinkedHashMap<>();
		bestClique = null;
		
		if(bix==null) { //Do it once for each level
			if(!createReducedInstance(level)) {
				return;
			}
		}
		
		double dist = 0;
		if(!bix.isBranch()) {
			Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
			dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
			dist = Math.sqrt(dist);
		}

		if(enumerate) {
			findAllCliques(level);
		}
		else {
			findOptimalClique(level);
		}
	}

	/**
	 * Use cluster enumeration to find all cliques that have an objective lower than the optimality gap
	 * @param level
	 * @param network
	 */
	private void findAllCliques(int level) {
		countEnumerate = 0;
		
		boolean cont = true;
		while(cont) {
			if(countEnumerate>=GlobalParam.ENUMERATE_PP_NUM_COLUMNS_MAXIMUM_IN_GRAPH) {
				break;
			}
		}
	}
	
	/**
	 * Find all cliques containing best vertex
	 * Then remove this vertex from network
	 * 
	 * @param level
	 * @param network
	 * @param bestVertex
	 */
	private void findAllCliquesContainingVertex(int level, Graph<Integer, Edge> network, int bestVertex) {
		// Find neighbors (excluding current vertex)
		Set<Integer> Gi = new LinkedHashSet<>();
		for(Edge e: network.edgesOf(bestVertex)) {
			Gi.add(e.getOrigin());
			Gi.add(e.getDestination());
		}
		Gi.remove(bestVertex);

		// Sort the neighbors
		List<Integer> Gi_list = new ArrayList<>(Gi);
		Collections.sort(Gi_list); //Sort

		// Initialise clique
		boolean[] currentClique = new boolean[instance.getNumPoints()];
		currentClique[bestVertex] = true;
		List<Integer> remainingVertices = new ArrayList<>(Gi_list);

		// Find all cliques that contain the vertex 
		double dualCost = - duals.getDualsMaxCluster()[level];
		recursiveSearch(network, level, currentClique, remainingVertices, bestVertex, 0, 0, dualCost);			

		// Remove the vertex
		network.removeVertex(bestVertex);
	}

	private void recursiveSearch(Graph<Integer, Edge> network, int level, boolean[] currentClique, List<Integer> remainingVertices, int currentVertex, int currentSize, double currentCost, double currentDualCost) {
		if(countEnumerate>=GlobalParam.ENUMERATE_PP_NUM_COLUMNS_MAXIMUM_IN_GRAPH) {
			return;
		}
		
		// Check if this is indeed a clique, only the last vertex has to be checked
		double nextCost = currentCost * currentSize;
		int nextSize = currentSize + 1;
		List<Integer> curPoints = new ArrayList<>();
		for(int i = 0; i < instance.getNumPoints(); i++) {
			if(currentClique[i] && i!=currentVertex) {
				if(!network.containsEdge(i, currentVertex)) {
					// Not a valid clique, stop
					return;
				}
				double dist = 0;
				if(!bix.isBranch()) {
					Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
					dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
				}
				nextCost += dist * instance.getSquaredDistanceMatrix()[i][currentVertex];
				curPoints.add(i);
			}
		}
		nextCost = nextCost / ((double) nextSize);

		// Actual current cost
		double nextDualCost = currentDualCost - duals.getDualsPoint()[level][currentVertex]; 
	
		// Check if we should add this cluster
		double obj = nextCost + nextDualCost;
		double actualObj = obj;
		Cluster cluster = new Cluster(currentClique);
		if(actualObj < GlobalParam.ENUMERATE_RC_COST_THRESHOLD) {
			countEnumerate++;
			resultPotential.put(cluster, obj); 
		}

		if(remainingVertices.isEmpty()) {
			return;
		}

		// Simulate what would happen if we always add the closest point (multiple times) and the largest dual (sorted)
		// If at least one of the situations is lower than optimality gap we continue
		// Only when there are 2 or more objects remaining
		if(remainingVertices.size()>1) {
			boolean boundSuccess = false;

			// Find the closest node 
			double shortestDistance = Double.POSITIVE_INFINITY;
			for(Integer j: remainingVertices) {
				double tempDistance = 0;
				for(Integer i: curPoints) {
					double dist = 0;
					if(!bix.isBranch()) {
						Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
						dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
					}
					tempDistance += dist * instance.getSquaredDistanceMatrix()[i][j];
				}
				if(tempDistance<shortestDistance) {
					shortestDistance = tempDistance;
				}
			}

			List<Double> remainingDuals = new ArrayList<>();
			for(Integer i: remainingVertices) {
				remainingDuals.add(duals.getDualsPoint()[level][i]); 
			}
			Collections.sort(remainingDuals, Collections.reverseOrder()); // Largest dual first

			// Start simulation
			double tempNextCost = nextCost;
			int tempNextSize = nextSize;
			double tempNextDualCost = nextDualCost;

			for(Double tempDual: remainingDuals) {
				tempNextCost = tempNextCost * tempNextSize;
				tempNextSize++;
				tempNextCost += shortestDistance;
				tempNextCost = tempNextCost / ((double) tempNextSize);			
				tempNextDualCost = tempNextDualCost - tempDual; 

				// Total cost
				if(tempNextCost + tempNextDualCost < GlobalParam.ENUMERATE_RC_COST_THRESHOLD) {
					boundSuccess = true;
				}
			}

			if(!boundSuccess) {
				return;
			}
		}

		List<Integer> nextRemainingVertices = new ArrayList<>(remainingVertices);

		// Add new vertex
		while(!nextRemainingVertices.isEmpty()) {
			int nextVertex = nextRemainingVertices.remove(0);
			boolean[] nextClique = currentClique.clone();
			nextClique[nextVertex] = true;
			recursiveSearch(network, level, nextClique, nextRemainingVertices, nextVertex, nextSize, nextCost, nextDualCost);
		}
	}

	private void findOptimalClique(int level) {
		// Step 1
		double bestObj = Double.POSITIVE_INFINITY;
		Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
		
		if(!bix.isBranch()) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				// Assume i is the medoid
				double obj = 0;
				List<Integer> pointId = new ArrayList<>();
				pointId.add(i);
				obj -= duals.getDualsPoint()[level][i];
				for(int j = 0; j < instance.getNumPoints(); j++) {
					if(i==j) {
						continue;
					}
					double dist = GlobalParam.DISTANCE_FUNCTION.getSquaredDistanceBetweenPoints(instance, i, j) * ObjectiveWeightFunction.getWeight(pair.first, pair.second);
					if(dist <= duals.getDualsPoint()[level][j]) {
						// We should add this point
						pointId.add(j);
						obj += dist - duals.getDualsPoint()[level][j];
					}
				}
				
				// Add all columns with actual negative RC 
				obj -= duals.getDualsMaxCluster()[level];
				if(obj < GlobalParam.NEGATIVE_RC_THRESHOLD) {
					// Create cluster
					Collections.sort(pointId);
					int[] pointId_array = pointId.stream().mapToInt(Integer::intValue).toArray();
					Cluster cluster = new Cluster(pointId_array, instance.getNumPoints());
					resultPotential.put(cluster, obj);
				}
				if(obj < bestObj) {
					bestObj = obj;
					bestClique = new Triple<>(pointId, obj, obj);
				}
			}
		}
		else {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				// Assume i is the medoid
				double obj = 0;
				List<Integer> pointId = new ArrayList<>();
				pointId.add(i);
				obj -= duals.getDualsPoint()[level][i];
				
				// Add all the nodes enforced by i
				for(List<Integer> list: bix.getBranchMapping()) {
					if(list.contains(i)) {
						for(Integer j: list) {
							if(i==j) {
								continue;
							}
							double dist = GlobalParam.DISTANCE_FUNCTION.getSquaredDistanceBetweenPoints(instance, i, j) * ObjectiveWeightFunction.getWeight(pair.first, pair.second);
							obj += dist - duals.getDualsPoint()[level][j];
							pointId.add(j);
						}
					}
				}
				
				for(List<Integer> list: bix.getBranchMapping()) {
					if(list.contains(i)) {
						// We already added these point
						continue;
					}
					
					// We have to decide whether to add all points in list
					double dist = 0;
					double dualValue = 0;
					for(Integer j: list) {
						dist += GlobalParam.DISTANCE_FUNCTION.getSquaredDistanceBetweenPoints(instance, i, j) * ObjectiveWeightFunction.getWeight(pair.first, pair.second);
						dualValue += duals.getDualsPoint()[level][j];
					}
					
					if(dist <= dualValue) {
						// We should add this point
						pointId.addAll(list);
						obj += dist - dualValue;
					}
				}
				
				// Check if any forbidden variables are violated
				Collections.sort(pointId);
				boolean violationDetected = false;
				outerloop: for(Integer m: pointId) {
					for(Integer n: pointId) {
						if(m>=n) {
							continue;
						}
						
						if(bi.getForbidden(level).contains(new Pair<>(m, n))) {
							violationDetected = true;
							break outerloop;
						}
					}
				}
				if(violationDetected) {
					GlobalParam.COUNTER_KNAPSACK++;
					// Solve a knapsack problem
					CliqueFinder_Knapsack cliqueFinder = null;
					try {
						cliqueFinder = new CliqueFinder_Knapsack(instance, level, i, duals, bi, bix);
						cliqueFinder.solve();
						Triple<List<Integer>, Double, Double> clique = cliqueFinder.getClique();
						obj = clique.second;
						pointId = clique.first;
						Collections.sort(pointId);
						cliqueFinder.cleanUp();
					} catch (IloException e) {
						e.printStackTrace();
					}
				}
				
				// Add all columns with actual negative RC 
				obj -= duals.getDualsMaxCluster()[level];
				if(obj < GlobalParam.NEGATIVE_RC_THRESHOLD) {
					// Create cluster
					int[] pointId_array = pointId.stream().mapToInt(Integer::intValue).toArray();
					Cluster cluster = new Cluster(pointId_array, instance.getNumPoints());
					resultPotential.put(cluster, obj);
				}
				if(obj < bestObj) {
					bestObj = obj;
					bestClique = new Triple<>(pointId, obj, obj);
				}
			}
		}
	}

	/**
	 * Officially what we calculate here is beta_i.
	 * However, alpha_i = beta_i^2
	 * And the bound is defined as sqrt(alpha_i)
	 * So we simply return beta_i
	 * 
	 * @param level
	 * @param i
	 * @return
	 */
	private double calcClusterEnumerationIndividualBound(int level, int i) {
		return (-2*duals.getDualsSqrtPoint()[level][i]+ Math.sqrt(4*duals.getDualsPoint()[level][i]+4*GlobalParam.ENUMERATE_GRAPH_THRESHOLD)) / 2d;
	}

	public Map<Integer, Double> getBestCliquePerLevel() {
		return bestCliquePerLevel;
	}
	
	@Override
	public Map<Integer, Map<Cluster, Double>> getRC() {
		throw new IllegalArgumentException("Not implemented yet");
	}
}