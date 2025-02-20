package cl.hierarchical.model.extended;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
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
import cl.data.Solution;
import cl.data.type.DistanceType;
import cl.data.type.PricingProblemStrategy;
import cl.data.type.SolverType;
import cl.data.type.ValidInequalityType;
import cl.kmeans.util.GenerateAllClusters;
import cl.util.Pair;
import cl.util.Triple;
import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;

public class MIP_Colgen {
	private static Logger log = LoggerFactory.getLogger(MIP_Colgen.class);

	private Instance instance;
	private IloCplex model;
	private IloObjective objective;

	private boolean useVariablesDummy = false; // We only need dummy variables when branching

	private Map<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> variablesCluster; // (Start level, end level)
	private Map<Pair<Integer, Integer>, IloNumVar> variablesDummyPoint;
	private Map<Integer, IloNumVar> variablesDummy; // level
	private IloRange[][] constraintsPoint;
	private IloRange[] constraintsMerge;
	private IloRange[] constraintsStart;

	private SolverType solverType;

	private PricingProblemStrategy pricingProblemStrategy = PricingProblemStrategy.None;
	private PricingProblem pp;
	
	// Column management
	private Map<Pair<Integer, Integer>, Set<Cluster>> tabooList;
	private int numRemovedClusters;

	// Column generation
	private int numIter, numColumns;
	private double objectiveLB;

	// Branch and price
	private ArrayDeque<Triple<BranchInformation, List<Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean>>, Branch>> queue;
	private Map<Pair<Integer, Integer>, List<Entry<Cluster, IloNumVar>>> removedClusters;
	private List<Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean>> branchMap;
	private double bestObj = Double.MAX_VALUE;
	private double bestLP;
	private boolean firstBranch = true;
	private int numBranchAndPrice;
	private BinaryTreeLowerBound binaryTree;
	private Branch currentBranch;

	// Start solution
	private Solution bestStartSolution;
	private List<Double> startSolutionObjectives;

	// Solution
	private Solution bestSolution;

	// Computation time
	private long totalTime, rootNodeTime, pricingTime, masterTime;

	private boolean solveInteger = false;
	
	private boolean timeLimitReached = false;

	public MIP_Colgen(Instance instance, SolverType solverType) throws IloException {
		this(instance, solverType, PricingProblemStrategy.None, false, true);
	}
	
	public MIP_Colgen(Instance instance, SolverType solverType, PricingProblemStrategy pricingProblemStrategy) throws IloException {
		this(instance, solverType, pricingProblemStrategy, false, true);
	}

	public MIP_Colgen(Instance instance, SolverType solverType, PricingProblemStrategy pricingProblemStrategy, boolean useVariablesDummy, boolean resetGlobalCounters) throws IloException {
		this.instance = instance;
		this.instance.initStartEndLevelList();
		this.solverType = solverType;
		this.pricingProblemStrategy = pricingProblemStrategy;
		if(useVariablesDummy || solverType==SolverType.Colgen_BranchAndPrice || solverType==SolverType.Colgen_CuttingPlanes) {
			this.useVariablesDummy = true;
		}
		this.model = new IloCplex();
		this.model.setOut(null);	
		initModel(resetGlobalCounters);	
	}

	private void initModel(boolean resetGlobalCounters) throws IloException {
		initVariables();
		if(useVariablesDummy) {
			initDummyVariables();
		}

		// Constraints
		initConstraintsPoint();
		initConstraintsMaxEndCluster();
		initConstraintsStartCluster();

		// Objective
		initObjective();

		initPricingProblems();
		initEmptyLists(resetGlobalCounters);
		
		// Initialise first cluster
		initFirstCluster();
	}

	private void initFirstCluster() throws IloException {
		boolean[] pointId = new boolean[instance.getNumPoints()];
		Arrays.fill(pointId, Boolean.TRUE);
		
		Cluster cluster = new Cluster(pointId);
		addCluster(new Pair<>(0, 0), cluster, true);
	}
	
	private void initDummyVariables() throws IloException {
		variablesDummy = new LinkedHashMap<>();	
		for(int i = 0; i < instance.getNumClusters(); i++) {
			variablesDummy.put(i, model.numVar(0, instance.getNumClusters()));
		}

		variablesDummyPoint = new LinkedHashMap<>();
		for(int i = 0; i < instance.getNumClusters(); i++) {
			for(int j = 0; j < instance.getNumPoints(); j++) {
				variablesDummyPoint.put(new Pair<>(i, j), model.numVar(0, 1));
			}
		}
	}

	public void setSolveInteger(boolean solveInteger) {
		this.solveInteger = solveInteger;
	}

	public void solve() throws IloException {
		totalTime = System.currentTimeMillis();
		if(solverType==SolverType.Colgen_BranchAndPrice) {
			solveBranchAndPrice();	
		}
		else if(solverType==SolverType.Colgen_Enumerate) {
			solveBranchAndPrice();	
		}
		else if(solverType==SolverType.Colgen_CuttingPlanes) {
			solveRootNode();
		}
		else if(solverType==SolverType.Colgen_RootNode) {
			solveRootNode();
		}
		else if(solverType==SolverType.Colgen_Integer) {
			if(solveInteger) {
				model.solve();
				bestSolution = getCurrentSolution();
				bestSolution.setFeasible(true);
			}
			else {
				solveInteger();
			}
		}
		else if(solverType==SolverType.Colgen_GetDualInformation) {
			model.solve();
		}
		else {
			throw new IllegalArgumentException("SolverType "+ solverType + " is not implemented (yet)");
		}
		totalTime = System.currentTimeMillis() - totalTime;
	}

	private void solveInteger() throws IloException {
		log.info("Solve integer");
		generateColumns();
		MIP_Colgen intModel = new MIP_Colgen(instance, SolverType.Colgen_Integer);
		intModel.setSolveInteger(true);
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			for(Cluster c: entry.getValue().keySet()) {
				intModel.addCluster(entry.getKey(), c, true);
			}
		}
		intModel.solve();
		intModel.printSolution();
		intModel.printVariables();
		bestSolution = intModel.getBestSolution(false);
		intModel.cleanUp();
	}

	private void solveBranchAndPrice() throws IloException {
		// Solve root node
		log.info("Solve root node");
		if(bestStartSolution!=null) {
			bestSolution = bestStartSolution;
			bestSolution.setFeasible(true);
			bestObj = bestSolution.getObjective();
		}
		else {
			bestObj = Double.POSITIVE_INFINITY;
		}
		queue = new ArrayDeque<>();
		queue.add(new Triple<>(new BranchInformationColgen(instance), new ArrayList<>(), null));
		
		while(!queue.isEmpty()) {
			if(GlobalParam.USE_GLOBAL_TIME_LIMIT && (System.currentTimeMillis() - totalTime > GlobalParam.GLOBAL_TIME_LIMIT)) {
				timeLimitReached = true;
				return;
			}
			
			Triple<BranchInformation, List<Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean>>, Branch> branch = queue.pollLast();
			solveBranchAndPrice(branch);
			
			numBranchAndPrice++;
		}
	}

	/**
	 * Return the best LP bound
	 * @param bi
	 * @return
	 * @throws IloException
	 */
	private void solveBranchAndPrice(Triple<BranchInformation, List<Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean>>, Branch> triple) throws IloException {
		BranchInformation bi = triple.first;
		Branch parentBranch = triple.third;
		
		// Undo previous branch changes
		if(parentBranch!=null) {
			for(Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean> branch: branchMap) {
				bi.clear(branch.first);
			}
			undoRemoveViolatedClusters(removedClusters);
		}
		removedClusters = new LinkedHashMap<>();
		branchMap = triple.second;

		if(parentBranch!=null) {
			for(Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean> branch: branchMap) {
				if(branch.second) {
					bi.enforce(branch.first);
					removeViolatedClusters(branch.first, true);
				}
				else {
					bi.forbid(branch.first);
					removeViolatedClusters(branch.first, false);
				}
			}
			Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean> lastBranch = branchMap.get(branchMap.size()-1);
			currentBranch = new Branch(lastBranch.first, lastBranch.second);
			
			log.info("Current branch {}, {}", lastBranch.first, lastBranch.second);
		}
		
		model.solve();

		// Calculate the best gap. If the gap is small enough we stop branching.
		double bestGap = 100d * (bestObj - bestLP) / bestLP;
		log.info("BnP iter: {}, Best LP: {}, Best obj: {}, Best gap: {}%", numBranchAndPrice, bestLP, bestObj, bestGap);
		if(bestGap < GlobalParam.OPTIMALITY_GAP) { 
			log.info("Finished branch and price. Best LP : " + bestLP + ", Best obj: " + bestObj + ", Best gap: " + bestGap + "%");
			return;
		}

		if(firstBranch) {
			rootNodeTime = System.currentTimeMillis();
		}
		else {
			pp.setBranchInformation(bi);
		}
		generateColumns();
		
		if(GlobalParam.PRINT_DETAILS) {
			printVariables(); 
		}

		// Calculate current gap
		double curObj = model.getObjValue();
		double gap = 100d * (bestObj - curObj) / curObj;
		log.info("Finished generating columns. Current LP: " + curObj + ", Best obj: " + bestObj + ", Cur gap: " + gap + "%");

		if(firstBranch) { //This should only happen once after solving the rootNode
			firstBranch = false;
			bestLP = curObj;
			rootNodeTime = System.currentTimeMillis() - rootNodeTime;
			parentBranch = new Branch(new Pair<Integer, Pair<Integer, Integer>>(-1, new Pair<Integer, Integer>(-1, -1)), false);
			binaryTree = new BinaryTreeLowerBound(parentBranch, bestLP);
		}
		else {
			if(currentBranch.isEnforce()) {
				binaryTree.addRight(parentBranch, currentBranch, curObj);
			}
			else {
				binaryTree.addLeft(parentBranch, currentBranch, curObj);
				binaryTree.updateBestLP(currentBranch);
				log.info("Found a new bestLP of {}, previous value was {}", binaryTree.getBestLP(), bestLP);
				bestLP = binaryTree.getBestLP();
			}
			parentBranch = currentBranch;
		}

		Map<Pair<Integer, Cluster>, Double> cMap = getClusterNonIntegrality();
		if (cMap.isEmpty()) {
			// Integral solution found!!
			log.info("Integral solution found of value: {}", curObj);
			if (bestSolution == null || curObj < bestObj) {
				log.info("Storing better integral solution of value: {}", curObj);
				bestObj = curObj;
				bestSolution = getCurrentSolution();
				bestSolution.setFeasible(true);
			}
			return;
		}

		if(solverType==SolverType.Colgen_Enumerate) {
			GlobalParam.ENUMERATE_RC_COST_THRESHOLD = Math.max(bestObj - curObj, 0); // Absolute gap
			GlobalParam.ENUMERATE_GRAPH_THRESHOLD = GlobalParam.ENUMERATE_RC_COST_THRESHOLD;

			log.info("Find all clusters with RC below {}", GlobalParam.ENUMERATE_RC_COST_THRESHOLD);

			int totalCountEnumerate = 0;
			model.solve();
			updateDuals();
			Map<Integer, Set<Cluster>> cols = pp.generateColumns(true, null);
			for(Entry<Integer, Set<Cluster>> entry: cols.entrySet()) {
				log.info("Found {} clusters for level {}", entry.getValue().size(), entry.getKey());
			}
			for(Entry<Integer, Set<Cluster>> entry: cols.entrySet()) {
				log.info("Adding {} clusters to level {}", entry.getValue().size(), entry.getKey());
				for(Cluster c: entry.getValue()) {
					addCluster(instance.getStartEndLevelList().get(entry.getKey()), c, true);
					totalCountEnumerate++;
				}
			}
			GlobalParam.COUNTER_ENUMERATE += totalCountEnumerate;

			System.out.println("Counter "+GlobalParam.COUNTER_ENUMERATE);

			log.info("Found {} clusters with RC below {}", totalCountEnumerate, GlobalParam.ENUMERATE_RC_COST_THRESHOLD);
			if(GlobalParam.APPLY_CUTTING_PLANES_AFTER_ENUMERATE) {
				log.info("Apply cutting planes after enumerate");
				model.solve();
				printSolution();
				printVariables();
				applyCuttingPlanes();

				cMap = getClusterNonIntegrality();
				if(cMap.isEmpty()) {
					bestObj = model.getObjValue();
					bestSolution = getCurrentSolution();
					bestSolution.setFeasible(true);
					log.info("Storing better integral solution of value: {}", bestObj);
					return;
				}
			}
			else {
				log.info("Solve integer model after enumerate");
				MIP_Colgen intModel = new MIP_Colgen(instance, SolverType.Colgen_Integer, PricingProblemStrategy.None, false, false);
				intModel.setSolveInteger(true);
				for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
					for(Cluster c: entry.getValue().keySet()) {
						intModel.addCluster(entry.getKey(), c, true);
					}
				}
				intModel.solve();
				bestSolution = intModel.getBestSolution(false);
				intModel.cleanUp();
				return;
			}
		}

		if (bestSolution != null && curObj >= bestObj) {
			log.info("Bound - solution is worse than optimal");
			// Although this solution is not integer, it is worse than the upper bound.
			return;
		}
		
		if(clusterLargerThanOne()) {
			log.info("A cluster is selected more than once, this cannot be an optimal LP solution");
			return;
		}

		// Start the branching over here
		Pair<Integer, Pair<Integer, Integer>> branch = suggestHighestLevelLargestClusterBranch(cMap); 

		// First enforce then forbid (since it is a queue)
		log.info("Create branching on {}, forbid", branch); // zi + zj <= 1
		List<Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean>> branchMapForbid = new ArrayList<>(branchMap);
		branchMapForbid.add(new Pair<>(branch, false));
		queue.add(new Triple<>(bi, branchMapForbid, parentBranch));

		log.info("Create branching on {}, enforce", branch); // zi = zj
		List<Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean>> branchMapEnforce = new ArrayList<>(branchMap);
		branchMapEnforce.add(new Pair<>(branch, true));
		queue.add(new Triple<>(bi, branchMapEnforce, parentBranch));
	}

	private boolean levelsHaveOverlap(Pair<Integer, Integer> level1, Pair<Integer, Integer> level2) {
		if(level1.first >= level2.first && level1.second <= level2.first) {
			return true;
		}
		if(level2.first >= level1.first && level2.second <= level1.first) {
			return true;
		}
		return false;
	}

	/**
	 * When enforce that 2 points are in the same cluster, it must also hold for all higher levels
	 * 
	 * When forbidding that 2 points are in the same cluster, it must also hold for all lower levels
	 * 
	 * Enforce or forbid
	 * @param branch
	 * @param enforce
	 * @throws IloException 
	 */
	private void removeViolatedClusters(Pair<Integer, Pair<Integer, Integer>> branch, boolean enforce) throws IloException {
		if(enforce) {
			for(Pair<Integer, Integer> removeLevel: instance.getStartEndLevelList()) {
				if(removeLevel.second<= branch.first && removeLevel.first>= branch.first) {
					List<Entry<Cluster, IloNumVar>> clusterRemove = new ArrayList<>();
					for(Entry<Cluster, IloNumVar> entry: variablesCluster.get(removeLevel).entrySet()) {
						Cluster cluster = entry.getKey();
						Pair<Integer, Integer> pair = branch.second;
						// If a cluster only contains one of the points, the cluster is not valid
						if((cluster.getCluster()[pair.first]!=cluster.getCluster()[pair.second])) {
							clusterRemove.add(entry);
							entry.getValue().setUB(0);
						}
					}
					if(!removedClusters.containsKey(removeLevel)) {
						removedClusters.put(removeLevel, new ArrayList<>());
					}
					removedClusters.get(removeLevel).addAll(clusterRemove);
				}
			}
		}
		else {
			for(Pair<Integer, Integer> removeLevel: instance.getStartEndLevelList()) {
				if(removeLevel.second<= branch.first && removeLevel.first>= branch.first) {
					List<Entry<Cluster, IloNumVar>> clusterRemove = new ArrayList<>();
					for(Entry<Cluster, IloNumVar> entry: variablesCluster.get(removeLevel).entrySet()) {
						Cluster cluster = entry.getKey();
						Pair<Integer, Integer> pair = branch.second;
						// If a cluster contains both of the points, the cluster is not valid
						if((cluster.getCluster()[pair.first] && cluster.getCluster()[pair.second])) {
							clusterRemove.add(entry);
							entry.getValue().setUB(0);
						}
					}
					if(!removedClusters.containsKey(removeLevel)) {
						removedClusters.put(removeLevel, new ArrayList<>());
					}
					removedClusters.get(removeLevel).addAll(clusterRemove);
				}
			}
		}
	}
	
	private boolean clusterLargerThanOne() throws UnknownObjectException, IloException {
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			for(Entry<Cluster, IloNumVar> entry2: entry.getValue().entrySet()) {
				double val = model.getValue(entry2.getValue());
				if(val > 1 + GlobalParam.EPSILON3) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Set UB to 1 again
	 * @param removed
	 * @throws IloException 
	 */
	private void undoRemoveViolatedClusters(Map<Pair<Integer, Integer>, List<Entry<Cluster, IloNumVar>>> removed) throws IloException {
		for(Entry<Pair<Integer, Integer>, List<Entry<Cluster, IloNumVar>>> entry: removed.entrySet()) {
			for(Entry<Cluster, IloNumVar> entry2: entry.getValue()) {
				entry2.getValue().setUB(Double.POSITIVE_INFINITY);
			}
		}
	}

	private Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<IloNumVar, Double>> getClusterNonIntegralityOriginal() throws UnknownObjectException, IloException {
		Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<IloNumVar, Double>> result = new LinkedHashMap<>();
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			for(Entry<Cluster, IloNumVar> entry2: entry.getValue().entrySet()) {
				double val = model.getValue(entry2.getValue());
				if(isFractional(val)) {
					Pair<Pair<Integer, Integer>, Cluster> cluster = new Pair<>(entry.getKey(), entry2.getKey());
					result.put(cluster, new Pair<>(entry2.getValue(), val));
				}
			}
		}
		return result;
	}

	/**
	 * These fractional variables have the same cluster in common (at a different level)
	 * @return
	 * @throws UnknownObjectException
	 * @throws IloException
	 */
	private Map<Cluster, Map<Pair<Integer, Integer>, Pair<IloNumVar, Double>>> getClusterNonIntegralitySimilar() throws UnknownObjectException, IloException {
		Map<Cluster, Map<Pair<Integer, Integer>, Pair<IloNumVar, Double>>> result = new LinkedHashMap<>();
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			for(Entry<Cluster, IloNumVar> entry2: entry.getValue().entrySet()) {
				double val = model.getValue(entry2.getValue());
				if(isFractional(val)) {
					Cluster cluster = entry2.getKey();
					if(!result.containsKey(cluster)) {
						result.put(cluster, new LinkedHashMap<>());
					}
					result.get(cluster).put(entry.getKey(), new Pair<>(entry2.getValue(), val));
				}
			}
		}
		return result;
	}
	
	private Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, IloNumVar>>>> getClusterNonIntegralityPerConstraint() throws UnknownObjectException, IloException {
		Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, IloNumVar>>>> result = new LinkedHashMap<>();
		for(int i = 0; i < instance.getNumClusters(); i++) {
			Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, IloNumVar>>> levelResult = new LinkedHashMap<>();
			for(int j = 0; j < instance.getNumPoints(); j++) {
				Map<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> pointResult = new LinkedHashMap<>();
				levelResult.put(j, pointResult);
			}
			result.put(i, levelResult);
		}
		
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			for(Entry<Cluster, IloNumVar> entry2: entry.getValue().entrySet()) {
				double val = model.getValue(entry2.getValue());
				if(val>GlobalParam.EPSILON3) {
					for(int level = entry.getKey().second; level <= entry.getKey().first; level++) {
						for(int i: entry2.getKey().getPointId()) {
							if(!result.get(level).get(i).containsKey(entry.getKey())) {
								result.get(level).get(i).put(entry.getKey(), new LinkedHashMap<>());
							}
							result.get(level).get(i).get(entry.getKey()).put(entry2.getKey(), entry2.getValue());
						}
					}
				}
			}
		}
		return result;
	}

	private Map<Pair<Integer, Cluster>, Double> getClusterNonIntegrality() throws UnknownObjectException, IloException {
		Map<Pair<Integer, Cluster>, Double> result = new LinkedHashMap<>();
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			for(Entry<Cluster, IloNumVar> entry2: entry.getValue().entrySet()) {
				double val = model.getValue(entry2.getValue());
				for(int level = entry.getKey().second; level <= entry.getKey().first; level++) {
					Pair<Integer, Cluster> cluster = new Pair<>(level, entry2.getKey());
					if(!result.containsKey(cluster)) {
						result.put(cluster, 0d);
					}
					result.replace(cluster, result.get(cluster)+val);
				}
			}
		}

		Map<Pair<Integer, Cluster>, Double> finalResult = new LinkedHashMap<>();
		for(Entry<Pair<Integer, Cluster>, Double> entry: result.entrySet()) {
			if(isFractional(entry.getValue())) {
				finalResult.put(entry.getKey(), entry.getValue());
			}
		}
		return finalResult;
	}

	/**
	 * Fractional
	 * @param val
	 * @return
	 */
	private boolean isFractional(double val) {
		double rnd = Math.round(val);
		if (Math.abs(rnd - val) > GlobalParam.EPSILON3) {
			return true;
		}
		return false;
	}

	private boolean isFractional(IloNumVar var) throws UnknownObjectException, IloException {
		double val = model.getValue(var);
		return isFractional(val);
	}

	/**
	 * Also return the level
	 * Preferable, on the highest level T
	 * Find the largest fractional cluster, branch on 2 elements in that cluster that do not appear in another fractional cluster
	 * @param cMap 
	 * @return
	 */
	private Pair<Integer, Pair<Integer, Integer>> suggestHighestLevelLargestClusterBranch(Map<Pair<Integer, Cluster>, Double> cMap) {
		int highestLevel = -1;
		int bestNumElements = -1;
		Pair<Integer, Cluster> bestCluster = null;
		for (Entry<Pair<Integer, Cluster>, Double> entry: cMap.entrySet()) {
			double val = entry.getValue();
			if (val <= GlobalParam.EPSILON5 || val >= 1 - GlobalParam.EPSILON5) {
				continue;
			}
			if (entry.getKey().first>highestLevel && entry.getKey().second.getPointId().length>=2 ) { 
				bestCluster = entry.getKey();
				highestLevel = bestCluster.first;
				bestNumElements = bestCluster.second.getPointId().length;
			}
			else if(highestLevel==entry.getKey().first && entry.getKey().second.getPointId().length>=2 && entry.getKey().second.getPointId().length>bestNumElements) {
				bestCluster = entry.getKey();
				highestLevel = bestCluster.first;
				bestNumElements = bestCluster.second.getPointId().length;
			}
		}

		if(bestCluster.second.getPointId().length<2) {
			throw new IllegalArgumentException("This cluster contains less than 2 points, we cannot handle this situation (yet)");
		}

		// Find the second largest fractional cluster, with at least 1 element in common
		bestNumElements = -1;
		Pair<Integer, Cluster> secondCluster = null;
		for (Entry<Pair<Integer, Cluster>, Double> entry: cMap.entrySet()) {
			if(entry.getKey().first!=bestCluster.first) { // Only consider the level of the best cluster
				continue;
			}
			double val = entry.getValue();
			if (val <= GlobalParam.EPSILON5 || val >= 1 - GlobalParam.EPSILON5) {
				continue;
			}
			if (entry.getKey().second.getPointId().length>bestNumElements && !entry.getKey().second.equals(bestCluster.second)) { 
				// Check if they overlap at least 1 element
				boolean overlap = false;
				outerloop: for(Integer i: bestCluster.second.getPointId()) {
					for(Integer j: entry.getKey().second.getPointId()) {
						if(i==j) {
							overlap = true;
							break outerloop;
						}
					}
				}
				if(overlap) {
					secondCluster = entry.getKey();
					bestNumElements = secondCluster.second.getPointId().length;
				}
			}
		}

		// Select elements A and B from bestCluster, such that element A is in the second cluster but element B is not in the second cluster
		List<Integer> bestList = Arrays.stream(bestCluster.second.getPointId()).boxed().collect(Collectors.toList());
		List<Integer> secondList = Arrays.stream(secondCluster.second.getPointId()).boxed().collect(Collectors.toList());
		List<Integer> copyBestList = new ArrayList<>(bestList);
		bestList.removeAll(secondList);
		secondList.retainAll(copyBestList); // Intersection


		// Choose A from bestList, excluding secondList
		int idA = GlobalParam.RANDOM.nextInt(bestList.size());

		// Choose B from bestList, including secondList
		int idB = GlobalParam.RANDOM.nextInt(secondList.size());

		// First element should be the smallest
		int first = bestList.get(idA);
		int second = secondList.get(idB);
		if(second<first) {
			int temp = second;
			second = first;
			first = temp;
		}
		return new Pair<>(bestCluster.first, new Pair<>(first, second));
	}

	private void solveRootNode() throws IloException {
		generateColumns();
		bestSolution = getCurrentSolution();
		bestLP = model.getObjValue();
		Map<Pair<Integer, Cluster>, Double> cMap = getClusterNonIntegrality();
		if(cMap.isEmpty()) {
			bestSolution.setFeasible(true);
		}
	}

	public void initHierarchicalStartSolutions(List<Solution> startSolutions) throws IloException {
		log.info("Adding {} startsolutions", startSolutions.size());
		for(Solution startSolution : startSolutions) {
			initHierarchicalStartSolution(startSolution);
		}
	}

	public void initHierarchicalStartSolution(Solution startSolution) throws IloException {
		startSolutionObjectives.add(startSolution.getObjective());

		if(bestStartSolution==null || startSolution.getObjective() < bestStartSolution.getObjective()) {
			bestStartSolution = startSolution; 
			log.info("Found better start solution with objective {}", bestStartSolution.getObjective());
		}
		for(Entry<Integer, List<Cluster>> entry: startSolution.getHierarchicalClusters().entrySet()) {
			for(Cluster cluster: entry.getValue()) {
				addClusterToEachLevel(cluster, true);
			}
		}
	}
	
	private void performColumnManagement() throws UnknownObjectException, IloException {
		int count = 0;
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			Set<Cluster> tabooSet = tabooList.get(entry.getKey());
			for(Entry<Cluster, IloNumVar> cluster: entry.getValue().entrySet()) {
				if(removedClusters!=null && removedClusters.containsKey(entry.getKey()) && removedClusters.get(entry.getKey()).contains(cluster)) {
					continue;
				}
				
				if(tabooSet.contains(cluster.getKey())) {
					continue;
				}
				
				if(model.getReducedCost(cluster.getValue())>GlobalParam.COLUMN_MANAGEMENT_THRESHOLD && model.getValue(cluster.getValue())<GlobalParam.EPSILON3) {
					cluster.getKey().setInactiveIter(cluster.getKey().getInactiveIter()+1);
				}
				else {
					cluster.getKey().setInactiveIter(0);
				}
				
				if(cluster.getKey().getInactiveIter()>=GlobalParam.COLUMN_MANAGEMENT_NUM_ITER) {
					count++;
				}
			}
		}
		
		IloNumVar[] removeColumns = new IloNumVar[count];
		int i = 0;
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			Iterator<Entry<Cluster, IloNumVar>> iter = entry.getValue().entrySet().iterator();
			
			while (iter.hasNext()) {
				Entry<Cluster, IloNumVar> cluster = iter.next();
				if(cluster.getKey().getInactiveIter()>=GlobalParam.COLUMN_MANAGEMENT_NUM_ITER) {
					removeColumns[i] = cluster.getValue();
					iter.remove();
					tabooList.get(entry.getKey()).add(cluster.getKey());
					i++;
				}
			}
		}
		
		model.end(removeColumns);
		numRemovedClusters += count;
		log.info("Perform column management, we removed {} columns. In total we removed {} columns", count, numRemovedClusters);
	}

	private void generateColumns() throws IloException {
		boolean cont = true;
		GlobalParam.ONLY_USE_VNS = true;
		while(cont) {
			if(GlobalParam.USE_GLOBAL_TIME_LIMIT && (System.currentTimeMillis() - totalTime > GlobalParam.GLOBAL_TIME_LIMIT)) {
				timeLimitReached = true;
				return;
			}
			
			numIter++;
			if(GlobalParam.PRINT_DETAILS) {
				log.info("NumIter {}", numIter);
			}

			long startMasterTime = System.currentTimeMillis(); 
			model.solve();
			masterTime += (System.currentTimeMillis() - startMasterTime);
			
			double gap = (model.getObjValue() - objectiveLB) / objectiveLB * 100;
			if(GlobalParam.PRINT_DETAILS) {
				int numVars = 0;
				for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
					numVars += entry.getValue().size();
				}
				log.info("Actual number of variables in the model {}", numVars);
				log.info("Obj {}, objectiveLB {}, gap {}%, feasible {}", model.getObjValue(), objectiveLB, gap, model.isPrimalFeasible());
			}
			if(!model.isPrimalFeasible())
			{
				throw new IllegalArgumentException("Not feasible");
			}
			if(Math.abs(model.getObjValue() - objectiveLB) < GlobalParam.EPSILON5) {
				log.info("Obj {} is close enough to objectiveLB {}", model.getObjValue(), objectiveLB);
				break;
			}
			updateDuals();

			boolean addColumn = false;

			long startPricingTime = System.currentTimeMillis();
			Map<Integer, Set<Cluster>> cols = pp.generateColumns();
			pricingTime += (System.currentTimeMillis() - startPricingTime);

			if(GlobalParam.USE_COLUMN_MANAGEMENT) {
				performColumnManagement();
			}
			for(Entry<Integer, Set<Cluster>> entry: cols.entrySet()) {
				for(Cluster c: entry.getValue()) {
					boolean added = addCluster(instance.getStartEndLevelList().get(entry.getKey()), c, true);
					if(added) {
						numColumns++;
						addColumn = true;
					}
				}
			}
			if(GlobalParam.PRINT_DETAILS) {
				log.info("TotalNumColumns {}", numColumns);
			}

			if(!addColumn) {
				if(GlobalParam.ONLY_USE_VNS) {
					GlobalParam.ONLY_USE_VNS = false;
					continue;
				}
				if(solverType==SolverType.Colgen_CuttingPlanes && applyCuttingPlanes()) {
					GlobalParam.ONLY_USE_VNS = true;
					continue;
				}
				break;
			}
			// Reset these parameters
			GlobalParam.ONLY_USE_VNS = true;
		}			
		if(GlobalParam.USE_COLUMN_MANAGEMENT) { // Because we changed the model
			model.solve();
		}
	}

	private boolean applyCuttingPlanes() throws IloException {
		log.info("Start cutting planes procedure");

		// Find fractional variables
		boolean success = false;
		int count = 0;
		while(true) {
			model.solve();
			boolean added = addValidInequalities();
			if(!added) {
				break;
			}
			success = true;
			count++;
		}
		log.info("Finish cutting planes procedure after {} iterations", count);
		updateDuals();
		return success;
	}

	/**
	 * First use Type1, if no violations are found use Type2
	 * @return
	 * @throws IloException
	 */
	private boolean addValidInequalities() throws IloException {
		int numCuts = 0;

		// Retrieve current values of variables
		Map<Cluster, Map<Pair<Integer, Integer>, Pair<IloNumVar, Double>>> cMap = getClusterNonIntegralitySimilar();
		Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<IloNumVar, Double>> cMapOriginal = getClusterNonIntegralityOriginal();
		Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, IloNumVar>>>> cMapConstraint = getClusterNonIntegralityPerConstraint();
		Map<Pair<Integer, Integer>, Map<Cluster, Double>> valuesVariablesCluster = new LinkedHashMap<>();
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			Map<Cluster, Double> map = new LinkedHashMap<>();
			for(Entry<Cluster, IloNumVar> entry2: entry.getValue().entrySet()) {
				map.put(entry2.getKey(), model.getValue(entry2.getValue()));
			}
			valuesVariablesCluster.put(entry.getKey(), map);
		}

		if(GlobalParam.VALID_INEQUALITIES.contains(ValidInequalityType.Type1)) {
			numCuts += addValidInequalities1(cMap, valuesVariablesCluster, true);
		}
		if(numCuts==0 && GlobalParam.VALID_INEQUALITIES.contains(ValidInequalityType.Type2)) {
			numCuts += addValidInequalities2(cMapOriginal, cMapConstraint);
		}
		log.info("Added {} cuts", numCuts);
		GlobalParam.COUNTER_CUTTING_PLANES += numCuts;
		return numCuts!=0;
	}

	/**
	 * e.g. x_7^A + \sum_B x_6^B <= 1
	 * where clusters A and B both contain item i
	 * but cluster A is not a subset of B
	 * 
	 * @param cMap
	 * @return
	 * @throws UnknownObjectException
	 * @throws IloException
	 */
	private int addValidInequalities1(Map<Cluster, Map<Pair<Integer, Integer>, Pair<IloNumVar, Double>>> cMap, Map<Pair<Integer, Integer>, Map<Cluster, Double>> valuesVariablesCluster, boolean addOneCutPerViolation) throws IloException {
		log.info("Find valid inequalities of type 1");
		int numCuts = 0;
		for(Entry<Cluster, Map<Pair<Integer, Integer>, Pair<IloNumVar, Double>>> entry: cMap.entrySet()) {
			int highestLevel = Integer.MAX_VALUE;	
			for(Entry<Pair<Integer, Integer>, Pair<IloNumVar, Double>> entry2: entry.getValue().entrySet()) {
				if(entry2.getKey().first < highestLevel) {
					highestLevel = entry2.getKey().first;
				}
			}

			Cluster originalCluster = entry.getKey();
			for(int level = highestLevel; level > 1; level--) {
				// Add one cut per level
				for(int i: originalCluster.getPointId()) {
					// Add one cut per point
					double check = 0;
					IloNumExpr expr = model.constant(0);
					for(Entry<Pair<Integer, Integer>, Pair<IloNumVar, Double>> entry2: entry.getValue().entrySet()) {
						expr = model.sum(expr, entry2.getValue().first);
						check += valuesVariablesCluster.get(entry2.getKey()).get(originalCluster);
					}

					int numElements = 0;
					for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry2: variablesCluster.entrySet()) {
						if(entry2.getKey().first >= level && entry2.getKey().second <= level) {
							for(Entry<Cluster, IloNumVar> entry3: entry2.getValue().entrySet()) {
								Cluster cluster = entry3.getKey();
								if(Cluster.checkValidInequalityType1(originalCluster, cluster, i)) {
									numElements++;
									expr = model.sum(expr, entry3.getValue());
									check += valuesVariablesCluster.get(entry2.getKey()).get(cluster);
								}
							}
						}
					}
					if(numElements>0 && check > 1+GlobalParam.EPSILON3) {
						numCuts++;
						IloRange newConstraint = model.addLe(expr, 1);
						if(addOneCutPerViolation) {
							break;
						}
					}
				}
			}
		}
		return numCuts;
	}

	/**
	 * Create a graph containing all variables with positive coefficient
	 * An edge exists if two variables cannot be selected at the same time
	 * Find the largest clique
	 * 
	 * @param cMap
	 * @param valuesVariablesCluster
	 * @param addOneCutPerViolation
	 * @return
	 * @throws IloException
	 */
	private int addValidInequalities2(Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<IloNumVar, Double>> cMapOriginal,
			Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, IloNumVar>>>> cMapConstraint) throws IloException {
		log.info("Find valid inequalities of type 2");
		int numCuts = 0;

		CuttingPlanes<IloNumVar> cuttingPlanes = new CuttingPlanes<>(cMapOriginal, cMapConstraint);
		List<List<IloNumVar>> cuts = cuttingPlanes.findValidInequalities2();
		for(List<IloNumVar> cut: cuts) {
			IloNumExpr expr = model.constant(0);
			for(IloNumVar var: cut) {
				expr = model.sum(expr, var);
			}
			IloRange newConstraint = model.addLe(expr, 1);
			numCuts++;
		}
		
		return numCuts;
	}

	private void initEmptyLists(boolean resetGlobalCounters) {
		startSolutionObjectives = new ArrayList<>();

		if(GlobalParam.USE_COLUMN_MANAGEMENT) {
			tabooList = new LinkedHashMap<>();
			for(Pair<Integer, Integer> pair: instance.getStartEndLevelList()) {
				tabooList.put(pair, new LinkedHashSet<>());
			}
		}
		
		// Init time
		totalTime = 0;
		rootNodeTime = 0;
		pricingTime = 0;
		masterTime = 0;

		if(resetGlobalCounters) {
			// Reset counters
			GlobalParam.COUNTER_VNS = 0;
			GlobalParam.COUNTER_QP = 0;
			GlobalParam.COUNTER_KNAPSACK = 0;
			GlobalParam.COUNTER_CUTTING_PLANES = 0;
			GlobalParam.COUNTER_ENUMERATE = 0;
		}
		numBranchAndPrice = 0;
	}

	private void initAllClusters() throws IloException {
		GenerateAllClusters gen = new GenerateAllClusters(instance);
		List<Cluster> clusters = gen.generate();

		for(Pair<Integer, Integer> pair: instance.getStartEndLevelList()) {
			for(Cluster cluster: clusters) {
				addCluster(pair, cluster, false);
			}
		}
	}

	private boolean addClusterToEachLevel(Cluster cluster, boolean startSolution) throws IloException {
		boolean added = false;
		for(Pair<Integer, Integer> pair: instance.getStartEndLevelList()) {
			// Check if we are allowed to add this cluster to the given level
			if(pair.first==0 && instance.getNumPoints()!=cluster.getPointId().length) {
				continue;
			}

			if(addCluster(pair, cluster, startSolution)) {
				added = true;
			}
		}
		return added;
	}

	public boolean addCluster(Pair<Integer, Integer> pair, Cluster cluster, boolean startSolution) throws IloException {
		if (cluster == null) {
			return false;
		}
		if(cluster.getPointId().length==instance.getNumPoints() && !(pair.first==0 && pair.second==0)) {
			// Do not add the cluster containing all elements, except at level 0
			return false;
		}
		if(pair.second==0 && pair.first!=0) {
			return false;
		}
		if(variablesCluster.get(pair).containsKey(cluster)) {
			if (!startSolution && GlobalParam.DEBUG) {
				throw new IllegalArgumentException("Cluster "+ cluster + " was already added to level "+ pair);
			}
			return false;
		}
		if(GlobalParam.USE_COLUMN_MANAGEMENT) {
			cluster = cluster.copy();
		}
		double dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
		double clusterCost = instance.computeClusterCost(cluster) * dist;
		IloColumn col = model.column(objective, clusterCost);
		for(Integer j: cluster.getPointId()) {
			for(int level = pair.first; level >= pair.second; level--) {
				col = col.and(model.column(constraintsPoint[level][j], 1));
			}
		}
		col = col.and(model.column(constraintsMerge[pair.second], 1));
		col = col.and(model.column(constraintsStart[pair.first], 1));

		IloNumVar var;
		if(solveInteger) {
			var = model.boolVar(col);
		}
		else {
			var = model.numVar(col, 0, Double.POSITIVE_INFINITY);
		}
		variablesCluster.get(pair).put(cluster, var);

		return true;
	}

	private void initPricingProblems() {
		if(GlobalParam.DISTANCE_TYPE==DistanceType.Centroid) {
			pp = new PricingProblemColgen(instance, pricingProblemStrategy);
		}
		else {
			pp = new PricingProblemColgenMedoid(instance);
		}
	}
	
	public Triple<double[][], double[], double[]> getOriginalDuals() throws UnknownObjectException, IloException {
		double[][] dualsPointTemp = new double[instance.getNumClusters()][instance.getNumPoints()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			for(int j = 0; j < instance.getNumPoints(); j++) {
				dualsPointTemp[i][j] = model.getDual(constraintsPoint[i][j]);
			}
		}
		
		double[] dualsMergeTemp = new double[instance.getNumClusters()]; 
		for(int i = 0; i < instance.getNumClusters(); i++) {
			dualsMergeTemp[i] = model.getDual(constraintsMerge[i]);
		}
		double[] dualsStartTemp = new double[instance.getNumClusters()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			dualsStartTemp[i] = model.getDual(constraintsStart[i]); // These duals have a different sign
		}
		
		return new Triple<>(dualsPointTemp, dualsMergeTemp, dualsStartTemp);
	}

	private DualInformation getDuals() throws IloException {
		double[][] dualsPointTemp = new double[instance.getNumClusters()][instance.getNumPoints()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			for(int j = 0; j < instance.getNumPoints(); j++) {
				dualsPointTemp[i][j] = model.getDual(constraintsPoint[i][j]);
			}
		}

		double[][] dualsPoint = new double[instance.getStartEndLevelList().size()][instance.getNumPoints()];
		double[][] dualsSqrtPoint = new double[instance.getStartEndLevelList().size()][instance.getNumPoints()];
		for(int i = 0; i < instance.getStartEndLevelList().size(); i++) {
			Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(i);
			for(int j = 0; j < instance.getNumPoints(); j++) {
				double temp = 0;
				for(int m = pair.second; m <= pair.first; m++) {
					temp += dualsPointTemp[m][j];
				}
				dualsPoint[i][j] = temp;
				dualsSqrtPoint[i][j] = Math.sqrt(dualsPoint[i][j]);
			}
		}

		double[] dualsMergeTemp = new double[instance.getNumClusters()]; 
		for(int i = 0; i < instance.getNumClusters(); i++) {
			dualsMergeTemp[i] = model.getDual(constraintsMerge[i]);
		}
		double[] dualsStartTemp = new double[instance.getNumClusters()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			dualsStartTemp[i] = model.getDual(constraintsStart[i]); // These duals have a different sign
		}

		double[] dualsMerge = new double[instance.getStartEndLevelList().size()];
		for(int i = 0; i < instance.getStartEndLevelList().size(); i++) {
			Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(i);
			dualsMerge[i] = dualsMergeTemp[pair.second] + dualsStartTemp[pair.first];
		}
		return new DualInformation(dualsPoint, dualsSqrtPoint, dualsMerge, null, null);
	}

	private void updateDuals() throws IloException {
		pp.setDuals(getDuals());
	}

	/**
	 * Constraint 1
	 * Each point is in at least one cluster
	 * @throws IloException 
	 */
	private void initConstraintsPoint() throws IloException {
		constraintsPoint = new IloRange[instance.getNumClusters()][instance.getNumPoints()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			for(int j = 0; j < instance.getNumPoints(); j++) {
				IloNumExpr expr = model.constant(0);
				if(useVariablesDummy) {
					expr = model.sum(expr, variablesDummyPoint.get(new Pair<>(i, j)));
				}

				IloRange constraint = model.addGe(expr, 1); //TODO: Eq or Ge
				if(GlobalParam.DEBUG) {
					constraint.setName("constraintsPoint_"+i+"_"+j);
				}
				constraintsPoint[i][j] = constraint;
			}
		}
	}

	/**
	 * Constraint 2
	 * Each level exactly 2 clusters merge (except at level 0)
	 * sum x_c <= 2
	 * @throws IloException 
	 */
	private void initConstraintsMaxEndCluster() throws IloException {
		constraintsMerge = new IloRange[instance.getNumClusters()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			IloNumExpr expr = model.constant(0);
			int maxCluster = 2;
			if(i==0) {
				maxCluster = 1;
			}
			constraintsMerge[i] = model.addLe(expr, maxCluster); //TODO: Le or Eq
		}
	}

	/**
	 * Constraint 3
	 * At each level, exactly 1 cluster can start
	 * @throws IloException
	 */
	private void initConstraintsStartCluster() throws IloException {
		constraintsStart = new IloRange[instance.getNumClusters()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			IloNumExpr expr = model.constant(0);
			if(useVariablesDummy) {
				expr = model.sum(expr, variablesDummy.get(i));
			}
			int maxCluster = 1;
			if(i==instance.getNumClusters()-1) {
				maxCluster = instance.getNumClusters();
			}
			constraintsStart[i] = model.addGe(expr, maxCluster); //TODO: Le or Eq
		}
	}

	private void initObjective() throws IloException {
		IloNumExpr expr = model.constant(0);
		if(useVariablesDummy) {
			for(int i = 0; i < instance.getNumClusters(); i++) {
				IloNumExpr temp = model.prod(GlobalParam.BIGM6, variablesDummy.get(i));
				expr = model.sum(expr, temp);
			}
			
			for(int i = 0; i < instance.getNumClusters(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					IloNumExpr temp = model.prod(GlobalParam.BIGM6, variablesDummyPoint.get(new Pair<>(i, j)));
					expr = model.sum(expr, temp);
				}
			}
		}
		objective = model.addMinimize(expr);
	}

	public void cleanUp() throws IloException {
		model.clearModel();
		model.end();
	}

	private void initVariables() {
		variablesCluster = new LinkedHashMap<>();
		for(Pair<Integer, Integer> pair: instance.getStartEndLevelList()) {
			variablesCluster.put(pair, new LinkedHashMap<>());
		}
	}

	public Solution getBestSolution(boolean setDetailedInformation) {
		bestSolution.setSolverType(solverType);
		bestSolution.setPricingProblemStrategy(pricingProblemStrategy);
		bestSolution.setTimeLimitReached(timeLimitReached);

		if(setDetailedInformation) {
			// Time
			bestSolution.setTotalTime(totalTime);
			bestSolution.setMasterTime(masterTime);
			bestSolution.setPricingTime(pricingTime);
			bestSolution.setRootNodeTime(rootNodeTime);

			// Solutions
			double bestStartObjective = Double.POSITIVE_INFINITY; 
			if(bestStartSolution!=null) {
				bestStartObjective = bestStartSolution.getObjective();
			}
			bestSolution.setBestStartObjective(bestStartObjective);
			bestSolution.setStartObjectives(startSolutionObjectives);
			bestSolution.setBestLP(bestLP);

			// Counters
			bestSolution.setNumIter(numIter);
			bestSolution.setNumColumns(numColumns);
			bestSolution.setNumBranchAndPrice(numBranchAndPrice);
			bestSolution.setNumBranchAndBound(GlobalParam.COUNTER_BRANCH_AND_BOUND);
			bestSolution.setNumVNS(GlobalParam.COUNTER_VNS);
			bestSolution.setNumQP(GlobalParam.COUNTER_QP);
			bestSolution.setNumKnapsack(GlobalParam.COUNTER_KNAPSACK);
			bestSolution.setNumCuttingPlanes(GlobalParam.COUNTER_CUTTING_PLANES);
			bestSolution.setNumEnumerate(GlobalParam.COUNTER_ENUMERATE);
		}

		return bestSolution;
	}

	private Solution getCurrentSolution() throws UnknownObjectException, IloException {
		Map<Integer, List<Cluster>> hierarchicalClusters = new LinkedHashMap<>();
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			List<Cluster> clusters = new ArrayList<>();
			for(Entry<Cluster, IloNumVar> entry2: entry.getValue().entrySet()) {
				if(model.getValue(entry2.getValue())>0.5) {
					clusters.add(entry2.getKey());
				}
			}
			for(int i = entry.getKey().second; i <= entry.getKey().first; i++) {
				if(!hierarchicalClusters.containsKey(i)) {
					hierarchicalClusters.put(i, new ArrayList<>());
				}
				hierarchicalClusters.get(i).addAll(clusters);
			}
		}
		return new Solution(instance, hierarchicalClusters, model.getObjValue());
	}

	public void printSolution() throws IloException {
		System.out.println("is feaisble? "+model.isPrimalFeasible());
		System.out.println("objective: "+ model.getObjValue());
	}

	public void printVariables() throws IloException {
		System.out.println("Variables");
		for(Entry<Pair<Integer, Integer>, Map<Cluster, IloNumVar>> entry: variablesCluster.entrySet()) {
			System.out.println("Level: "+entry.getKey());
			for(Entry<Cluster, IloNumVar> entry2: entry.getValue().entrySet()) {
				if(model.getValue(entry2.getValue())>GlobalParam.EPSILON3) {
					System.out.println(model.getValue(entry2.getValue()) +  " "+ entry2.getKey());
				}
			}
		}
	}

}
