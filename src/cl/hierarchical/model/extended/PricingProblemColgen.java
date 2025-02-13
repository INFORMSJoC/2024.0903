package cl.hierarchical.model.extended;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cl.data.Cluster;
import cl.data.GlobalParam;
import cl.data.Instance;
import cl.data.ObjectiveWeightFunction;
import cl.data.type.PricingProblemStrategy;
import cl.hierarchical.model.hp.CliqueFinder_Dinkelbach;
import cl.util.Pair;
import cl.util.Triple;
import ilog.concert.IloException;


public class PricingProblemColgen implements PricingProblem {
	private static Logger log = LoggerFactory.getLogger(PricingProblemColgen.class);
	
	private Instance instance;
	private PricingProblemStrategy pricingProblemStrategy;
	private DualInformation duals;
	private BranchInformation bi;
	private Set<Integer> exclusionSet = new LinkedHashSet<>();
	private Set<Integer> restrictedGraph = new LinkedHashSet<>(); // Use for PP strategy with restricted graph

	private Map<Cluster, Double> resultPotential; //Cluster and reducedCost
	private Triple<List<Integer>, Double, Double> bestClique;
	private Map<Integer, Double> bestCliquePerLevel; // Only store the best objective
	private double bestRC;
	private Map<Integer, Map<Cluster, Double>> allRC;
	
	// Use for enumeration
	private int countEnumerate;

	// Used for branching
	private BranchInformationExtra bix;

	// For pricing problem strategy
	private List<List<Integer>> pricingProblemStrategyOrder;
	
	private int columnLimit = GlobalParam.PP_NUM_COLUMNS;

	public PricingProblemColgen(Instance instance, PricingProblemStrategy pricingProblemStrategy) {
		this.instance = instance;
		this.pricingProblemStrategy = pricingProblemStrategy;

		if(pricingProblemStrategy!=PricingProblemStrategy.None) {
			initPricingProblemStrategyOrder();
		}
	}

	public void setColumnLimit(int columnLimit) {
		this.columnLimit = columnLimit;
	}

	public double getBestRC() {
		return bestRC;
	}

	private void initPricingProblemStrategyOrder() {
		pricingProblemStrategyOrder = new ArrayList<>();

		Map<Pair<Integer, Integer>, Integer> startEndLevelMap = new LinkedHashMap<>();
		for(int i = 1; i < instance.getStartEndLevelList().size(); i++) { // We skip the first cluster, since it contains all elements
			if(instance.getStartEndLevelList().get(i).second==0) { // We skip clusters that end at level 0
				continue;
			}
			startEndLevelMap.put(instance.getStartEndLevelList().get(i), i);			
		}

		if(pricingProblemStrategy==PricingProblemStrategy.EndFirstStartSecond) {
			List<Integer> singleRestart = new ArrayList<>();
			for(int i = 1; i < instance.getNumClusters(); i++) {
				for(int j = 1; j <= i; j++) {
					singleRestart.add(startEndLevelMap.get(new Pair<>(i, j)));
				}
			}
			pricingProblemStrategyOrder.add(singleRestart);
		}
		else if(pricingProblemStrategy==PricingProblemStrategy.EndFirstStartSecondWithRestart) {
			for(int i = 1; i < instance.getNumClusters(); i++) {
				List<Integer> singleRestart = new ArrayList<>();
				for(int j = 1; j <= i; j++) {
					singleRestart.add(startEndLevelMap.get(new Pair<>(i, j)));
				}
				pricingProblemStrategyOrder.add(singleRestart);
			}
		}
		else if(pricingProblemStrategy==PricingProblemStrategy.StartFirstEndSecond) {
			List<Integer> singleRestart = new ArrayList<>();
			for(int i = instance.getNumClusters()-1; i > 0; i--) {
				for(int j = i; j > 0; j--) {
					singleRestart.add(startEndLevelMap.get(new Pair<>(i, j)));
				}
			}
			pricingProblemStrategyOrder.add(singleRestart);
		}
		else if(pricingProblemStrategy==PricingProblemStrategy.StartFirstEndSecondWithRestart) {
			for(int i = instance.getNumClusters()-1; i > 0; i--) {
				List<Integer> singleRestart = new ArrayList<>();
				for(int j = i; j > 0; j--) {
					singleRestart.add(startEndLevelMap.get(new Pair<>(i, j)));
				}
				pricingProblemStrategyOrder.add(singleRestart);
			}
		}
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
		Map<Integer, Set<Cluster>> result = new LinkedHashMap<>();
		allRC = new LinkedHashMap<>();
		bestCliquePerLevel = new LinkedHashMap<>();
		bestRC = Double.POSITIVE_INFINITY;

		if(pricingProblemStrategy==PricingProblemStrategy.None || enumerate) {
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
							resultPotential.entrySet().stream()
							.sorted(Map.Entry.comparingByValue())
							.limit(limit)
							.collect(Collectors.toMap(
									Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new)); //Ascending order
					result.put(i, sorted.keySet());
				}
				else {
					Map<Cluster,Double> sorted =
							resultPotential.entrySet().stream()
							.sorted(Map.Entry.comparingByValue())
							.limit(columnLimit)
							.collect(Collectors.toMap(
									Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new)); //Ascending order
					result.put(i, sorted.keySet());
					allRC.put(i, sorted);
					
					for(Double obj: sorted.values()) {
						if(obj < bestRC) {
							bestRC = obj;
						}
					}
				}
			}
		}
		else if(GlobalParam.PP_RESTRICTED_GRAPH) {
			for(List<Integer> list: pricingProblemStrategyOrder) {
				// Do a restart
				boolean first = true;
				Set<Integer> restartResult = new LinkedHashSet<>();
				
				for(Integer level : list) {
					Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
					
					boolean performPP = false;
					if(first) {
						performPP = true;
					}
					else {
						
						try {
							generatePotentialColumns(level, enumerate);
						} catch (IloException e) {
							e.printStackTrace();
						}
						Map<Cluster,Double> sorted =
								resultPotential.entrySet().stream()
								.sorted(Map.Entry.comparingByValue())
								.limit(GlobalParam.PP_NUM_COLUMNS)
								.collect(Collectors.toMap(
										Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new)); //Ascending order
						result.put(level, sorted.keySet());
						for(Cluster cluster: sorted.keySet()) {
							restrictedGraph.addAll(Arrays.stream(cluster.getPointId()).boxed().collect(Collectors.toList()));
						}
						if(sorted.isEmpty()) {
							performPP = true;
						}
					}
					
					if(performPP) {
						restrictedGraph = new LinkedHashSet<>();
						try {
							generatePotentialColumns(level, enumerate);
						} catch (IloException e) {
							e.printStackTrace();
						}
						Map<Cluster,Double> sorted =
								resultPotential.entrySet().stream()
								.sorted(Map.Entry.comparingByValue())
								.limit(GlobalParam.PP_NUM_COLUMNS)
								.collect(Collectors.toMap(
										Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new)); //Ascending order
						result.put(level, sorted.keySet());
						if(first) {
							restartResult = new LinkedHashSet<>();
						}
						for(Cluster cluster: sorted.keySet()) {
							restrictedGraph.addAll(Arrays.stream(cluster.getPointId()).boxed().collect(Collectors.toList()));
						}
					}
					first = false;
				}
			}
		}
		else {
			for(List<Integer> list: pricingProblemStrategyOrder) {
				// Do a restart
				boolean first = true;
				Map<Cluster, Double> restartResult = new LinkedHashMap<>();

				for(Integer level : list) {
					Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
					double dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);

					// Reset BIX
					bix = null;

					boolean performPP = false;
					if(first) {
						performPP = true;
					}
					else {
						Map<Cluster, Double> temp = new LinkedHashMap<>();
						if(initReducedInstanceAndBranchMapping(level)) {
							// Check if restartResult contains negative RC columns that do not violate branching rules
							for(Entry<Cluster, Double> entry: restartResult.entrySet()) {
								// Calculate new RC
								double RC = - duals.getDualsMaxCluster()[level];
								for(int i: entry.getKey().getPointId()) {
									for(int j: entry.getKey().getPointId()) {
										if(i>=j) {
											continue;
										}
										RC += (dist * instance.getSquaredDistanceMatrix()[i][j] - duals.getDualsPoint()[level][i]);
									}
								}
								RC = RC / ((double) entry.getKey().getPointId().length);
								if(RC > GlobalParam.NEGATIVE_RC_THRESHOLD) {
									continue;
								}
	
								boolean violated = false;
	
								if(bix.isBranch()) {
									// Does it violate enforced?
									List<Integer> pointId = Arrays.stream(entry.getKey().getPointId()).boxed().collect(Collectors.toList());
									outerloop: for(Integer i: pointId) {
										for(List<Integer> enforced: bix.getBranchMapping()) {
											if(enforced.contains(i)) {
												if(!pointId.containsAll(enforced)) {
													violated = true;
													break outerloop;
												}
											}
										}
									}
	
									// Does it violate forbidden?
									outerloop: for(int i: entry.getKey().getPointId()) {
										for(int j: entry.getKey().getPointId()) {
											if(i>=j) {
												continue;
											}
											Pair<Integer, Integer> pair2 = new Pair<>(i, j);
											if(bi.getForbidden(level).contains(pair2)) {
												violated = true;
												break outerloop;
											}
										}
									}
								}
								if(!violated) {
									temp.put(entry.getKey(), RC);
								}
							}
						}
						if(!temp.isEmpty()) {
							result.put(level, temp.keySet());
						}
						else {
							performPP = true;
						}
					}

					if(performPP) {
						try {
							generatePotentialColumns(level, enumerate);
						} catch (IloException e) {
							e.printStackTrace();
						}
						Map<Cluster,Double> sorted =
								resultPotential.entrySet().stream()
								.sorted(Map.Entry.comparingByValue())
								.limit(GlobalParam.PP_NUM_COLUMNS)
								.collect(Collectors.toMap(
										Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new)); //Ascending order
						result.put(level, sorted.keySet());
						if(first) {
							restartResult = sorted;
						}
					}
					first = false;
				}

			}
		}
		log.info("Lowest RC {}", bestRC);
		return result;
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
	 * Return false if there are inconsistencies
	 * @param level
	 * @return
	 */
	private boolean initReducedInstanceAndBranchMapping(int level) {
		resultPotential = new LinkedHashMap<>();
		bestClique = null;

		if(bix==null) { //Do it once for each level
			if(!createReducedInstance(level)) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Algorithm 2 from Aloise (2012)
	 * @throws IloException 
	 */
	private void generatePotentialColumns(int level, boolean enumerate) throws IloException {
		if(!initReducedInstanceAndBranchMapping(level)) {
			return;
		}

		double dist = 0;
		if(!bix.isBranch()) {
			Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
			dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
			dist = Math.sqrt(dist);
		}

		// Initialise network
		boolean[][] network = new boolean[bix.getReducedInstance().getNumPoints()][bix.getReducedInstance().getNumPoints()];
		Map<Integer, Integer> vertexSet = new LinkedHashMap<>();
		
		for(int i = 0; i < bix.getReducedInstance().getNumPoints(); i++) {
			if(exclusionSet.contains(i)) {
				continue;
			}
			vertexSet.put(i, 0);
		}
		for(int i = 0; i < bix.getReducedInstance().getNumPoints(); i++) {
			for(int j = i+1 ; j < bix.getReducedInstance().getNumPoints(); j++) {
				if(exclusionSet.contains(i) || exclusionSet.contains(j)) {
					continue;
				}
				if(!bix.isBranch()) {
					if(enumerate) {
						if(dist * instance.getRegularDistanceMatrix()[i][j]<=duals.getDualsSqrtPoint()[level][i] + duals.getDualsSqrtPoint()[level][j] + calcClusterEnumerationIndividualBound(level, i) + calcClusterEnumerationIndividualBound(level, j)) { //TODO: Rounding errors may cause problems here
							network[i][j] = true;
							vertexSet.replace(i, vertexSet.get(i)+1);
							vertexSet.replace(j, vertexSet.get(j)+1);
						}	
					}
					else {
						if(dist * instance.getRegularDistanceMatrix()[i][j]<=duals.getDualsSqrtPoint()[level][i] + duals.getDualsSqrtPoint()[level][j]) { //TODO: Rounding errors may cause problems here
							network[i][j] = true;
							vertexSet.replace(i, vertexSet.get(i)+1);
							vertexSet.replace(j, vertexSet.get(j)+1);
						}
					}
				}
				else {
					if(bix.getReducedInstance().getRegularDistanceMatrix()[i][j]<=bix.getReducedDuals().getDualsSqrtPointSingle()[i] + bix.getReducedDuals().getDualsSqrtPointSingle()[j]) { //TODO: Rounding errors may cause problems here
						network[i][j] = true;
						vertexSet.replace(i, vertexSet.get(i)+1);
						vertexSet.replace(j, vertexSet.get(j)+1);
					}
				}
			}
		}

		if(enumerate) {
			findAllCliques(level, network, vertexSet);
		}
		else {
			findOptimalClique(level, network, vertexSet);
		}
	}

	/**
	 * Use cluster enumeration to find all cliques that have an objective lower than the optimality gap
	 * @param level
	 * @param network
	 */
	private void findAllCliques(int level, boolean[][] network, Map<Integer, Integer> vertexSet) {
		countEnumerate = 0;

		boolean cont = true;
		while(cont) {
			if(countEnumerate>=GlobalParam.ENUMERATE_PP_NUM_COLUMNS_MAXIMUM_IN_GRAPH) {
				break;
			}
			if(vertexSet.isEmpty()) {
				break;
			}

			// Find vertex with minimum degree
			int bestVertex = -1;
			int bestDegree = Integer.MAX_VALUE;
			for(Entry<Integer, Integer> entry: vertexSet.entrySet()) {
				int degree = entry.getValue();
				if(degree < bestDegree) {
					bestDegree = degree;
					bestVertex = entry.getKey();
				}
			}

			findAllCliquesContainingVertex(level, network, vertexSet, bestVertex);
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
	private void findAllCliquesContainingVertex(int level, boolean[][] network, Map<Integer, Integer> vertexSet, int bestVertex) {
		// Find neighbors (excluding current vertex)
		Set<Integer> Gi = findNeighbours(network, vertexSet, bestVertex);
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
		removeVertex(network, vertexSet, bestVertex);
	}
	
	private Set<Integer> findNeighbours(boolean[][] network, Map<Integer, Integer> vertexSet, int bestVertex) {
		Set<Integer> Gi = new LinkedHashSet<>();
		for(int i = 0; i < bix.getReducedInstance().getNumPoints(); i++) {
			if(i < bestVertex) {
				if(network[i][bestVertex]) {
					Gi.add(i);
				}
			}
			else {
				if(network[bestVertex][i]) {
					Gi.add(i);
				}
			}
		}
		return Gi;
	}
	
	private void removeVertex(boolean[][] network, Map<Integer, Integer> vertexSet, int removeVertex) {
		vertexSet.remove(removeVertex);
		for(int i = 0; i < bix.getReducedInstance().getNumPoints(); i++) {
			if(i < removeVertex) {
				network[i][removeVertex] = false;
			}
			else {
				network[removeVertex][i] = false;
			}
		}
	}
	
	private boolean networkContainsEdge(boolean[][] network, int i, int j) {
		if(i < j) {
			return network[i][j];
		}
		return network[j][i];
	}

	private void recursiveSearch(boolean[][] network, int level, boolean[] currentClique, List<Integer> remainingVertices, int currentVertex, int currentSize, double currentCost, double currentDualCost) {
		if(countEnumerate>=GlobalParam.ENUMERATE_PP_NUM_COLUMNS_MAXIMUM_IN_GRAPH) {
			return;
		}

		// Check if this is indeed a clique, only the last vertex has to be checked
		double nextCost = currentCost * currentSize;
		int nextSize = currentSize + 1;
		List<Integer> curPoints = new ArrayList<>();
		for(int i = 0; i < instance.getNumPoints(); i++) {
			if(currentClique[i] && i!=currentVertex) {
				if(!networkContainsEdge(network, i, currentVertex)) {
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

	private void findOptimalClique(int level, boolean[][] network, Map<Integer, Integer> vertexSet) {
		// Step 1
		double bestObj = Double.POSITIVE_INFINITY;
		boolean cont = true;
		while(cont) {
			if(vertexSet.isEmpty()) {
				break;
			}
		
			// Step 1a
			int bestVertex = -1;
			int bestDegree = Integer.MAX_VALUE;
			if(restrictedGraph.isEmpty()) {
				for(Entry<Integer, Integer> entry: vertexSet.entrySet()) {
					int degree = entry.getValue();
					if(degree < bestDegree) {
						bestDegree = degree;
						bestVertex = entry.getKey();
					}
				}
			}
			else {
				for(Entry<Integer, Integer> entry: vertexSet.entrySet()) {
					int degree = entry.getValue();
					if(degree < bestDegree && restrictedGraph.contains(entry.getKey())) {
						bestDegree = degree;
						bestVertex = entry.getKey();
					}
				}
				
				if(bestVertex==-1) {
					for(Entry<Integer, Integer> entry: vertexSet.entrySet()) {
						int degree = entry.getValue();
						if(degree < bestDegree) {
							bestDegree = degree;
							bestVertex = entry.getKey();
						}
					}
				}
			}
			
			// Step 1b
			Set<Integer> Gi = findNeighbours(network, vertexSet, bestVertex);
			Gi.add(bestVertex);

			// Step 1c
			List<Integer> Gi_list = new ArrayList<>(Gi);
			Collections.sort(Gi_list); //Sort
			CliqueFinder_Dinkelbach cliqueFinder = new CliqueFinder_Dinkelbach(instance, level, duals, Gi_list, bix); 
			cliqueFinder.solve();
			Triple<List<Integer>, Double, Double> pairQP = cliqueFinder.getClique();

			// Step 1d
			// Add all columns with actual negative RC 
			if(pairQP.third < GlobalParam.NEGATIVE_RC_THRESHOLD) {
				int[] pointId;
				if(!bix.isBranch()) {
					pointId = pairQP.first.stream().mapToInt(Integer::intValue).toArray();
				}
				else {
					// Convert to original instance
					List<Integer> oldIDs = new ArrayList<>();
					for(Integer i: pairQP.first) {
						oldIDs.addAll(bix.getBranchMapping().get(i));
					}
					Collections.sort(oldIDs);
					pointId = oldIDs.stream().mapToInt(Integer::intValue).toArray();
				}
				Cluster cluster = new Cluster(pointId, instance.getNumPoints());
				resultPotential.put(cluster, pairQP.third);
				if(GlobalParam.PP_STOP_WHEN_REACHING_POTENTIAL_COLUMNS && resultPotential.size()>GlobalParam.PP_NUM_POTENTIAL_COLUMNS) {
					GlobalParam.PP_STOP_EARLY = true;
					break;
				}
			}
			// Keep track of the optimal value of the pricing problem when setting kappa=0
			if(pairQP.second<bestObj) {
				bestObj = pairQP.second;
				bestClique = pairQP;
			}

			// Step 1e
			this.removeVertex(network, vertexSet, bestVertex);
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
		return allRC;
	}
}