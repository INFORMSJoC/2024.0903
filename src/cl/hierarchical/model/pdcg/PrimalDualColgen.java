package cl.hierarchical.model.pdcg;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.fasterxml.jackson.databind.ObjectMapper;

import cl.data.Cluster;
import cl.data.GlobalParam;
import cl.data.Instance;
import cl.data.ObjectiveWeight;
import cl.data.ObjectiveWeightFunction;
import cl.data.Solution;
import cl.data.type.DistanceType;
import cl.data.type.PricingProblemStrategy;
import cl.data.type.SolverType;
import cl.data.type.ValidInequalityType;
import cl.hierarchical.model.extended.BinaryTreeLowerBound;
import cl.hierarchical.model.extended.Branch;
import cl.hierarchical.model.extended.BranchInformation;
import cl.hierarchical.model.extended.BranchInformationColgen;
import cl.hierarchical.model.extended.CuttingPlanes;
import cl.hierarchical.model.extended.DualInformation;
import cl.hierarchical.model.extended.MIP_Colgen;
import cl.hierarchical.model.extended.PricingProblem;
import cl.hierarchical.model.extended.PricingProblemColgen;
import cl.hierarchical.model.extended.PricingProblemColgenMedoid;
import cl.util.MatrixOperations;
import cl.util.Pair;
import cl.util.Triple;
import cl.util.UtilFunctions;
import cl.util.run.SpreadsheetGenerator;
import ilog.concert.IloException;
import ilog.cplex.IloCplex.UnknownObjectException;

public class PrimalDualColgen {
	private static Logger log = LoggerFactory.getLogger(PrimalDualColgen.class);

	private Instance instance;
	private RealMatrix matrixA;
	private RealMatrix vectorB;
	private RealMatrix vectorC;
	
	private boolean useVariablesDummy;

	private Map<Pair<Integer, Integer>, Map<Cluster, PDCGM_ClusterInfo>> variablesClusterMap; // (Start level, end level)
	private List<PDCGM_ClusterInfo> variablesClusterList;

	// Start solution
	private List<Solution> startSolutions;
	private Solution bestStartSolution;
	private List<Double> startSolutionObjectives;

	private int numVars, numSlackVars, numCons, numCoverConstraints;

	// Colgen
	private SolverType solverType;
	private int numIter, numColumns;

	private PricingProblemStrategy pricingProblemStrategy;
	private PricingProblem pp;
	
	// Best solution
	private double bestObj;
	private Solution bestSolution;

	// Current solution
	private double curObj;
	private Solution currentSolution;
	
	// PDCGM
	private PrimalDualInteriorPointMethod PDIPM;
	private PrimalDualPoint optimalPoint;

	// Branch and price
	private ArrayDeque<Triple<BranchInformation, List<Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean>>, Branch>> queue;
	private List<PDCGM_ClusterInfo> removedClusters;
	private List<Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean>> branchMap;
	private double bestLP;
	private boolean firstBranch = true;
	private BinaryTreeLowerBound binaryTree;
	private Branch currentBranch;
	
	// Keep track of time
	private long totalTime, rootNodeTime, pricingTime, masterTime;
	private boolean timeLimitReached;
	private int numBranchAndPrice;
	
	private ObjectMapper om = new ObjectMapper();
	
	public PrimalDualColgen(Instance instance, SolverType solverType, PricingProblemStrategy pricingProblemStrategy) {
		this(instance, solverType, pricingProblemStrategy, false, true);
	}

	public PrimalDualColgen(Instance instance, SolverType solverType, PricingProblemStrategy pricingProblemStrategy, boolean useVariablesDummy, boolean resetGlobalCounters) {
		this.instance = instance;
		this.solverType = solverType;
		this.pricingProblemStrategy = pricingProblemStrategy;
		if(useVariablesDummy || solverType==SolverType.Colgen_BranchAndPrice || solverType==SolverType.Colgen_CuttingPlanes) {
			this.useVariablesDummy = true;
		}
		this.instance.initStartEndLevelList();

		initModel(resetGlobalCounters);
	}

	private void initModel(boolean resetGlobalCounters) {
		initPricingProblems();
		initEmptyLists(resetGlobalCounters);
		initConstraints();

		initFirstCluster();
	}

	private void initConstraints() {
		numVars = 0;

		// Constraints are:
		// objective (1)
		// cover constraints n*k
		// merge constraints k
		// start constraints k
		// var constraints

		numCons = (instance.getNumPoints()+2) * instance.getNumClusters();
		numCoverConstraints = instance.getNumPoints() * instance.getNumClusters();
		vectorB = new Array2DRowRealMatrix(numCons, 1);
		int count = 0;
		for(int i = 0; i < instance.getNumPoints(); i++) {
			for(int level = 0; level < instance.getNumClusters(); level++) {
				vectorB.setEntry(count, 0, 1);
				count++;
			}
		}
		for(int level = 0; level < instance.getNumClusters(); level++) {
			// Constraints merge
			if(level==0) {
				vectorB.setEntry(count, 0, 1);
			}
			else {
				vectorB.setEntry(count, 0, 2);
			}
			count++;
		}
		for(int level = 0; level < instance.getNumClusters(); level++) {
			// Constraints start
			if(level==instance.getNumClusters()-1) {
				vectorB.setEntry(count, 0, instance.getNumClusters());
			}
			else
			{
				vectorB.setEntry(count, 0, 1);
			}
			count++;
		}

		// Add slack vars 
		numSlackVars = 0; 
	}

	public void solve() {
		totalTime = System.currentTimeMillis();
		if(bestStartSolution!=null) {
			bestSolution = bestStartSolution;
			bestSolution.setFeasible(true);
			bestObj = bestSolution.getObjective();
			
			if(GlobalParam.WRITE_BETTER_INTEGER_SOLUTION) {
				writeBetterSolution();
			}
		}
		else {
			bestObj = Double.POSITIVE_INFINITY;
		}
		
		if(solverType==SolverType.Colgen_RootNode) {
			solveRootNode();
		}
		else if(solverType==SolverType.Colgen_CuttingPlanes) {
			solveRootNode();
		}
		else if(solverType==SolverType.Colgen_BranchAndPrice) {
			solveBranchAndPrice();
		}
		else if(solverType==SolverType.Colgen_Enumerate) {
			solveBranchAndPrice();	
		}
		else if(solverType==SolverType.Colgen_BranchPriceCutEnumerate) {
			solveBranchAndPrice();
		}
		else {
			throw new IllegalArgumentException("SolverType "+ solverType + " is not implemented (yet)");
		}
		totalTime = System.currentTimeMillis() - totalTime;
	}
	
	private void solveBranchAndPrice() {
		// Solve root node
		log.info("Solve root node");
				
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
	
	private void solveBranchAndPrice(Triple<BranchInformation, List<Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean>>, Branch> triple) {
		BranchInformation bi = triple.first;
		Branch parentBranch = triple.third;
		
		// Undo previous branch changes
		if(parentBranch!=null) {
			for(Pair<Pair<Integer, Pair<Integer, Integer>>, Boolean> branch: branchMap) {
				bi.clear(branch.first);
			}
			undoRemoveViolatedClusters(removedClusters);
		}
		removedClusters = new ArrayList<>();
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
			printClustersFromIPM(optimalPoint);
		}

		// Calculate current gap
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
				
				if(GlobalParam.WRITE_BETTER_INTEGER_SOLUTION) {
					writeBetterSolution();
				}
			}
			parentBranch = currentBranch;
		}

		Map<Pair<Integer, Cluster>, Double> cMap = getClusterNonIntegralityPerLevel(optimalPoint);
		if (cMap.isEmpty()) {
			// Integral solution found!!
			Solution currentSolution = getCurrentIntegerSolution(optimalPoint);
			log.info("Integral solution found of value: {}", currentSolution.getObjective());
			if (bestSolution == null || curObj < bestObj) {
				log.info("Storing better integral solution of value: {}", currentSolution.getObjective());
				bestObj = currentSolution.getObjective();
				bestSolution = currentSolution;
				
				if(GlobalParam.WRITE_BETTER_INTEGER_SOLUTION) {
					writeBetterSolution();
				}
			}
			return;
		}
		
		if(solverType==SolverType.Colgen_Enumerate) {
			GlobalParam.ENUMERATE_RC_COST_THRESHOLD = Math.max(bestObj - curObj, 0); // Absolute gap
			GlobalParam.ENUMERATE_GRAPH_THRESHOLD = GlobalParam.ENUMERATE_RC_COST_THRESHOLD;
		
			log.info("Find all clusters with RC below {}", GlobalParam.ENUMERATE_RC_COST_THRESHOLD);
			
			int totalCountEnumerate = 0;
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
				PDIPM = new PrimalDualInteriorPointMethod(matrixA, vectorB, vectorC, GlobalParam.EPSILON3, true);
				PDIPM.solve();
				applyCuttingPlanes();
				PrimalDualPoint curPoint = PDIPM.getSolutionPrimalDualPoint();
				Map<Pair<Integer, PDCGM_ClusterInfo>, Double> cMapEnumerate = getClusterNonIntegrality(curPoint);
				if(cMapEnumerate.isEmpty()) {
					bestSolution = getCurrentIntegerSolution(curPoint);
					bestObj = bestSolution.getObjective();
					log.info("Storing better integral solution of value: {}", bestObj);
					return;
				}
			}
			else {
				log.info("Solve integer model after enumerate");
				MIP_Colgen intModel;
				try {
					intModel = new MIP_Colgen(instance, SolverType.Colgen_Integer, PricingProblemStrategy.None, false, false);
					intModel.setSolveInteger(true);
					for(PDCGM_ClusterInfo clusterInfo: variablesClusterList) {
						if(clusterInfo!=null) {
							intModel.addCluster(clusterInfo.getLevel(), clusterInfo.getCluster(), true);
						}
					}
					intModel.solve();
					bestSolution = intModel.getBestSolution(false);
					intModel.cleanUp();
				} catch (IloException e) {
					e.printStackTrace();
				}
				return;
			}
		}
			
		if (bestSolution != null && curObj >= bestObj) {
			log.info("Bound - solution is worse than optimal");
			// Although this solution is not integer, it is worse than the upper bound.
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
	
	private void writeBetterSolution() {
		String date_file = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date());
		File outputFile = new File(GlobalParam.WRITE_BETTER_INTEGER_SOLUTION_FILE_NAME+"/java_"+date_file+".xlsx");
		List<Solution> tempList = new ArrayList<>();
		tempList.add(getBestSolution(true, System.currentTimeMillis() - totalTime));
		Map<File, List<Solution>> tempMap = new LinkedHashMap<>();
		tempMap.put(outputFile, tempList);
		try {
			SpreadsheetGenerator.writeSpreadsheet(outputFile, tempMap, true);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Set UB to 1 again
	 * @param removed
	 * @throws IloException 
	 */
	private void undoRemoveViolatedClusters(List<PDCGM_ClusterInfo> removed) {
		for(PDCGM_ClusterInfo clusterInfo: removed) {
			vectorC.setEntry(numSlackVars+clusterInfo.getId(), 0, clusterInfo.getClusterCost());
		}
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
	private void removeViolatedClusters(Pair<Integer, Pair<Integer, Integer>> branch, boolean enforce) {
		if(enforce) {
			for(PDCGM_ClusterInfo clusterInfo: variablesClusterList) {
				if(clusterInfo==null) {
					continue;
				}
				if(clusterInfo.getLevel().second <= branch.first && clusterInfo.getLevel().first >= branch.first) {
					Cluster cluster = clusterInfo.getCluster();
					Pair<Integer, Integer> pair = branch.second;
					// If a cluster only contains one of the points, the cluster is not valid
					if((cluster.getCluster()[pair.first]!=cluster.getCluster()[pair.second])) {
						removedClusters.add(clusterInfo);
						vectorC.setEntry(numSlackVars + clusterInfo.getId(), 0, GlobalParam.BIGM3);
					}
				}
			}
		}
		else {
			for(PDCGM_ClusterInfo clusterInfo: variablesClusterList) {
				if(clusterInfo==null) {
					continue;
				}
				if(clusterInfo.getLevel().second <= branch.first && clusterInfo.getLevel().first >= branch.first) {
					Cluster cluster = clusterInfo.getCluster();
					Pair<Integer, Integer> pair = branch.second;
					// If a cluster contains both of the points, the cluster is not valid
					if((cluster.getCluster()[pair.first] && cluster.getCluster()[pair.second])) {
						removedClusters.add(clusterInfo);
						vectorC.setEntry(numSlackVars + clusterInfo.getId(), 0, GlobalParam.BIGM3);
					}
				}
			}
		}
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
	
	private void solveRootNode() {
		generateColumns();
		
		Map<Pair<Integer, PDCGM_ClusterInfo>, Double> nonIntegrality = getClusterNonIntegrality(optimalPoint);
		printClustersFromIPM(optimalPoint);
		bestLP = curObj;
		if(nonIntegrality.isEmpty()) {
			bestSolution = getCurrentIntegerSolution(optimalPoint);
			log.info("Found an integral solution with value: {}", bestSolution.getObjective());
		}
		else {
			bestSolution = getCurrentSolution(optimalPoint);
			log.info("Found a fractional solution with value: {}", bestSolution.getObjective());
		}
	}

	/**
	 * Based on Large-scale optimization with the primal-dual column generation method
	 * Gondzio (2016)
	 */
	private void generateColumns() {
		// Intialise
		PrimalDualPoint curPoint = null;
		PrimalDualPoint startPrimalDualPoint = null;
		PDCGM_matrixAAT[][] matrixAAT = null;
		Pair<Map<Integer, List<Cluster>>, Double> solution = null;
		int previousIndex = -1;
		double epsilon = 0.5;
		double epsilonMax = 0.4;
		double D = 10+GlobalParam.EPSILON3;
		double UB = Double.POSITIVE_INFINITY;
		double LB = Double.NEGATIVE_INFINITY;
		double relGap = Double.POSITIVE_INFINITY;
		
		// BranchPriceCutEnumerate
		boolean useCuttingPlane = true;
		Map<Integer, Double> bounds = new LinkedHashMap<>();
		double lowestBound = Double.POSITIVE_INFINITY;

		GlobalParam.ONLY_USE_VNS = true;
		boolean primalFeasible = true;
		if(useVariablesDummy) {
			primalFeasible = false;
		}
		while(true) {
			if(relGap < GlobalParam.EPSILON3 &&
					(solverType==SolverType.Colgen_BranchAndPrice || solverType==SolverType.Colgen_Enumerate || solverType==SolverType.Colgen_RootNode  || solverType==SolverType.Colgen_BranchPriceCutEnumerate)) {
				break;
			}
			log.info("PDCGM iter {}", numIter);
			numIter++;
			
			// Step 2. Solve IPM
			PDIPM = new PrimalDualInteriorPointMethod(matrixA, vectorB, vectorC, epsilon, primalFeasible);
			if(startPrimalDualPoint!=null) {
				PDIPM.setStartPrimalDualPoint(startPrimalDualPoint);
			}
			if(matrixAAT!=null) {
				PDIPM.setMatrixAAT(matrixAAT, previousIndex);
			}
			long startMasterTime = System.currentTimeMillis(); 
			PDIPM.solve();
			masterTime += (System.currentTimeMillis() - startMasterTime);
			
			// Get IPM data
			curPoint = PDIPM.getSolutionPrimalDualPoint();
			matrixAAT = PDIPM.getMatrixAAT();
			previousIndex = matrixA.getColumnDimension();

			solution = getClustersFromIPM(curPoint);
			double obj = Double.POSITIVE_INFINITY;
			if(PDIPM.isPrimalFeasible()) {
				obj = solution.second;
			}
			double LB_dual = getDualLowerBound(curPoint);
			if(GlobalParam.PRINT_DETAILS) {
				log.info("Obj {}, LB dual {}", obj, LB_dual);
			}

			// Step 3. Update UB
			UB = Math.min(UB, obj);

			// Step 4. Generate columns
			pp.setDuals(getDuals(curPoint));
			GlobalParam.PP_STOP_EARLY = false;
			long startPricingTime = System.currentTimeMillis();
			Map<Integer, Set<Cluster>> cols = pp.generateColumns();
			pricingTime += (System.currentTimeMillis() - startPricingTime);
			boolean addColumn = false;
			for(Entry<Integer, Set<Cluster>> entry: cols.entrySet()) {
				for(Cluster c: entry.getValue()) {
					boolean added = addCluster(instance.getStartEndLevelList().get(entry.getKey()), c, true);
					if(added) {
						numColumns++;
						addColumn = true;
					}
					else {
						if(bounds.get(entry.getKey())==null) {
							bounds.put(entry.getKey(), pp.getRC().get(entry.getKey()).get(c));
						}
						else if(pp.getRC().get(entry.getKey()).get(c) < bounds.get(entry.getKey())) {
							bounds.replace(entry.getKey(), pp.getRC().get(entry.getKey()).get(c));
						}
						if(pp.getRC().get(entry.getKey()).get(c) < lowestBound) {
							lowestBound = pp.getRC().get(entry.getKey()).get(c); 
						}
					}
				}
			}
			if(GlobalParam.PRINT_DETAILS) {
				log.info("TotalNumColumns {}", numColumns);
			}

			// Step 5. Update LB
			if(!GlobalParam.PP_STOP_EARLY && !GlobalParam.ONLY_USE_VNS) {
				// The LB is only valid, if we used an exact QP and we did not stop the PP early
				double RC = pp.getBestRC();
				if(Double.isInfinite(RC)) {
					RC = 0;
				}
				LB = Math.max(LB, LB_dual + RC); 
			}

			if(PDIPM.isPrimalFeasible()) {
				// Step 6. Update gap
				relGap = (UB - LB) / (GlobalParam.EPSILON5 + Math.abs(UB)); 
	
				// Step 7. Update epsilon
				epsilon = Math.max(Math.min(epsilonMax, relGap / D), GlobalParam.EPSILON5);
			}
			if(GlobalParam.PRINT_DETAILS) {
				log.info("UB {}, LB {}, relGap {}, epsilon {}", UB, LB, relGap, epsilon);
			}

			// Step 8. Update PrimalDualPoint for next iteration
			if(addColumn) {
				startPrimalDualPoint = updateWarmStart(curPoint);
			}
			else {
				// Use an exact solver
				if(GlobalParam.ONLY_USE_VNS) {
					log.info("Switching to exact QP");
					GlobalParam.ONLY_USE_VNS = false;
					epsilon = Math.min(epsilon, GlobalParam.EPSILON5);
					if(solverType==SolverType.Colgen_CuttingPlanes) {
						epsilon = Math.min(epsilon, GlobalParam.EPSILON7);
					}
					primalFeasible = true;
					
					// Reset bounds
					bounds = new LinkedHashMap<>();
					lowestBound = Double.POSITIVE_INFINITY;
					continue;
				}
				if(solverType==SolverType.Colgen_CuttingPlanes && applyCuttingPlanes()) {
					curPoint = PDIPM.getSolutionPrimalDualPoint();
					solution = getClustersFromIPM(curPoint);
					UB = solution.second;
					break;
				}
				if(solverType==SolverType.Colgen_BranchPriceCutEnumerate && firstBranch && useCuttingPlane && applyCuttingPlanes(GlobalParam.BRANCH_PRICE_CUT_ENUMERATE_NUM_CUTTING_PLANES)) {
					startPrimalDualPoint = PDIPM.getSolutionPrimalDualPoint();
					matrixAAT = null; // Recompute from scratch
					GlobalParam.ONLY_USE_VNS = true;
					primalFeasible = false;
					relGap = Double.POSITIVE_INFINITY;
					UB = Double.POSITIVE_INFINITY;
					useCuttingPlane = false;
					continue;
				}
				if(solverType==SolverType.Colgen_BranchPriceCutEnumerate) {
					if(lowestBound < -GlobalParam.EPSILON2) {
						log.info("Find all clusters with RC below a level specific threshold");
						
						int totalCountEnumerate = 0;
						Map<Integer, Set<Cluster>> colsEnum = pp.generateColumns(true, bounds);
						for(Entry<Integer, Set<Cluster>> entry: colsEnum.entrySet()) {
							log.info("Found {} clusters for level {}", entry.getValue().size(), entry.getKey());
						}
						
						for(Entry<Integer, Set<Cluster>> entry: colsEnum.entrySet()) {
							log.info("Adding {} clusters to level {}", entry.getValue().size(), entry.getKey());
							
							if(!pp.getRC().containsKey(entry.getKey())) {
								continue;
							}
							Map<Cluster, Double> clusters = pp.getRC().get(entry.getKey());
							Map<Cluster, Double> sortedClusters =
									clusters.entrySet().stream()
								       .sorted(Map.Entry.comparingByValue(Comparator.reverseOrder()))
								       .collect(Collectors.toMap(
								          Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new));
							
							int count = 0;
							for(Cluster c: sortedClusters.keySet()) {
								if(count>GlobalParam.PP_NUM_COLUMNS) {
									break;
								}
								
								if(addCluster(instance.getStartEndLevelList().get(entry.getKey()), c, true)) {
									totalCountEnumerate++;
									count++;
								}
							}
						}
						
						if(totalCountEnumerate>0) {
							startPrimalDualPoint = updateWarmStart(curPoint);
							GlobalParam.ONLY_USE_VNS = true;
							GlobalParam.COUNTER_ENUMERATE += totalCountEnumerate;
							log.info("Found {} clusters with RC below {}", totalCountEnumerate, GlobalParam.ENUMERATE_RC_COST_THRESHOLD);
							continue;	
						}
					}
				}
				break;
			}

			// Reset these parameters
			GlobalParam.ONLY_USE_VNS = true;
		}

		// Update current objective
		curObj = UB;
		optimalPoint = curPoint;
		
		log.info("Finished generating columns with PDCGM, solution with obj "+ curObj);
	}
	
	private boolean applyCuttingPlanes() {
		return applyCuttingPlanes(Integer.MAX_VALUE);
	}
	
	private boolean applyCuttingPlanes(int maxNumIter) {
		log.info("Start cutting planes procedure");

		// Find fractional variables
		boolean success = false;
		int count = 0; 
		this.printClustersFromIPM(PDIPM.getSolutionPrimalDualPoint());
		while(true) {
			if(GlobalParam.USE_GLOBAL_TIME_LIMIT && (System.currentTimeMillis() - totalTime > GlobalParam.GLOBAL_TIME_LIMIT)) {
				timeLimitReached = true;
				return false;
			}
			if(count>=maxNumIter) {
				break;
			}
			PrimalDualPoint curPoint = PDIPM.getSolutionPrimalDualPoint();
			int numCuts = addValidInequalities(curPoint);
			if(numCuts==0) {
				break;
			}
			success = true;
			count++;
			
			PDIPM = new PrimalDualInteriorPointMethod(matrixA, vectorB, vectorC, GlobalParam.EPSILON5, true);
			PDIPM.setFeasibleEpsilon(GlobalParam.EPSILON3);
		
			RealMatrix newX = MatrixOperations.addRow(curPoint.getPointX(), 0.45, numCuts);
			RealMatrix newPi = MatrixOperations.addRow(curPoint.getPointPi(), 0.01, numCuts);
			RealMatrix newS = MatrixOperations.addRow(curPoint.getPointS(), 0.45, numCuts);
			PDIPM.setStartPrimalDualPoint(new PrimalDualPoint(newX, newPi, newS));
			PDIPM.solve();
		}
		log.info("Finish cutting planes procedure after {} iterations", count);
		return success;
	}
	
	/**
	 * These fractional variables have the same cluster in common (at a different level)
	 * @return
	 * @throws UnknownObjectException
	 * @throws IloException
	 */
	private Map<Cluster, Map<Pair<Integer, Integer>, Pair<Integer, Double>>> getClusterNonIntegralitySimilar(PrimalDualPoint curPoint) {
		Map<Cluster, Map<Pair<Integer, Integer>, Pair<Integer, Double>>> result = new LinkedHashMap<>();
		for(int i = numSlackVars; i < curPoint.getPointX().getRowDimension(); i++) {
			double val = curPoint.getPointX().getEntry(i, 0);
			if(isFractional(val)) {
				int actualIdx = i - numSlackVars;
				PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
				if(clusterInfo==null) {
					continue;
				}
				Cluster cluster = clusterInfo.getCluster();
				if(!result.containsKey(cluster)) {
					result.put(cluster, new LinkedHashMap<>());
				}
				result.get(cluster).put(clusterInfo.getLevel(), new Pair<>(i, val));
			}
		}
		return result;
	}
	
	private Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<Integer, Double>> getClusterNonIntegralityOriginal(PrimalDualPoint curPoint) {
		Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<Integer, Double>> result = new LinkedHashMap<>();
		for(int i = numSlackVars; i < curPoint.getPointX().getRowDimension(); i++) {
			double val = curPoint.getPointX().getEntry(i, 0);
			if(isFractional(val)) {
				int actualIdx = i - numSlackVars;
				PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
				if(clusterInfo==null) {
					continue;
				}
				Pair<Pair<Integer, Integer>, Cluster> cluster = new Pair<>(clusterInfo.getLevel(), clusterInfo.getCluster());
				result.put(cluster, new Pair<>(i, val));
			}
		}
		return result;
	}
	
	private Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, Integer>>>> getClusterNonIntegralityPerConstraint(PrimalDualPoint curPoint) {
		Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, Integer>>>> result = new LinkedHashMap<>();
		for(int i = 0; i < instance.getNumClusters(); i++) {
			Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, Integer>>> levelResult = new LinkedHashMap<>();
			for(int j = 0; j < instance.getNumPoints(); j++) {
				Map<Pair<Integer, Integer>, Map<Cluster, Integer>> pointResult = new LinkedHashMap<>();
				levelResult.put(j, pointResult);
			}
			result.put(i, levelResult);
		}
		
		for(int idx = numSlackVars; idx < curPoint.getPointX().getRowDimension(); idx++) {
			double val = curPoint.getPointX().getEntry(idx, 0);
			if(val>GlobalParam.EPSILON3) {
				int actualIdx = idx - numSlackVars;
				PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
				if(clusterInfo==null) {
					continue;
				}
				for(int level = clusterInfo.getLevel().second; level <= clusterInfo.getLevel().first; level++) {
					for(int i: clusterInfo.getCluster().getPointId()) {
						if(!result.get(level).get(i).containsKey(clusterInfo.getLevel())) {
							result.get(level).get(i).put(clusterInfo.getLevel(), new LinkedHashMap<>());
						}
						result.get(level).get(i).get(clusterInfo.getLevel()).put(clusterInfo.getCluster(), idx);
					}
				}
			}
		}
		return result;
	}
	
	/**
	 * First use Type1, if no violations are found use Type2
	 * @return
	 * @throws IloException
	 */
	private int addValidInequalities(PrimalDualPoint curPoint) {
		int numCuts = 0;

		// Retrieve current values of variables
		Map<Cluster, Map<Pair<Integer, Integer>, Pair<Integer, Double>>> cMap = getClusterNonIntegralitySimilar(curPoint);
		Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<Integer, Double>> cMapOriginal = getClusterNonIntegralityOriginal(curPoint);
		Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, Integer>>>> cMapConstraint = getClusterNonIntegralityPerConstraint(curPoint);
		Map<Pair<Integer, Integer>, Map<Cluster, Double>> valuesVariablesCluster = new LinkedHashMap<>();
		for(int idx = numSlackVars; idx < curPoint.getPointX().getRowDimension(); idx++) {
			double val = curPoint.getPointX().getEntry(idx, 0);
			int actualIdx = idx - numSlackVars;
			PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
			if(clusterInfo==null) {
				continue;
			}
			if(!valuesVariablesCluster.containsKey(clusterInfo.getLevel())) {
				valuesVariablesCluster.put(clusterInfo.getLevel(), new LinkedHashMap<>());
			}
			valuesVariablesCluster.get(clusterInfo.getLevel()).put(clusterInfo.getCluster(), val);
		}

		if(GlobalParam.VALID_INEQUALITIES.contains(ValidInequalityType.Type1)) {
			numCuts += addValidInequalities1(cMap, valuesVariablesCluster, curPoint, true);
		}
		if(numCuts==0 && GlobalParam.VALID_INEQUALITIES.contains(ValidInequalityType.Type2)) {
			numCuts += addValidInequalities2(cMapOriginal, cMapConstraint);
		}
		log.info("Added {} cuts", numCuts);
		GlobalParam.COUNTER_CUTTING_PLANES += numCuts;
		return numCuts;
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
	private int addValidInequalities1(Map<Cluster, Map<Pair<Integer, Integer>, Pair<Integer, Double>>> cMap, Map<Pair<Integer, Integer>, Map<Cluster, Double>> valuesVariablesCluster, PrimalDualPoint curPoint, boolean addOneCutPerViolation) {
		log.info("Find valid inequalities of type 1");
		int numCuts = 0;
		for(Entry<Cluster, Map<Pair<Integer, Integer>, Pair<Integer, Double>>> entry: cMap.entrySet()) {
			if(numCuts>=GlobalParam.CUTTING_PLANES_NUM_CUTS) {
				break;
			}
			int highestLevel = Integer.MAX_VALUE;	
			for(Entry<Pair<Integer, Integer>, Pair<Integer, Double>> entry2: entry.getValue().entrySet()) {
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
					List<Integer> cut = new ArrayList<>();
					for(Entry<Pair<Integer, Integer>, Pair<Integer, Double>> entry2: entry.getValue().entrySet()) {
						cut.add(entry2.getValue().first);
						check += valuesVariablesCluster.get(entry2.getKey()).get(originalCluster);
					}

					int numElements = 0;
					for(int idx = numSlackVars; idx < curPoint.getPointX().getRowDimension(); idx++) {
						double val = curPoint.getPointX().getEntry(idx, 0);
						int actualIdx = idx - numSlackVars;
						PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
						if(clusterInfo==null) {
							continue;
						}
						if(clusterInfo.getLevel().first >= level && clusterInfo.getLevel().second <= level) {
								Cluster cluster = clusterInfo.getCluster();
								if(Cluster.checkValidInequalityType1(originalCluster, cluster, i)) {
									numElements++;
									cut.add(idx);
									check += val;
								}
						}
					}
					if(numElements>0 && check > 1+GlobalParam.EPSILON2) {
						//TODO: make sure that we do not add the same constraint twice
						
						numCuts++;
						// Add a constraint
						RealMatrix rowA = new Array2DRowRealMatrix(numSlackVars+numVars, 1);
						List<PDCGM_ClusterInfo> clusters = new ArrayList<>();
						for(int idx: cut) {
							clusters.add(variablesClusterList.get(idx));
							rowA.setEntry(idx, 0, 1);
						}
						matrixA = MatrixOperations.addRow(matrixA, rowA);
						vectorB = MatrixOperations.addRow(vectorB, 1+GlobalParam.EPSILON5);
						numCons++;
						
						// Add a variable
						RealMatrix colA = new Array2DRowRealMatrix(numCons, 1);
						colA.setEntry(numCons-1, 0, 1);
						matrixA = MatrixOperations.addColumn(matrixA, colA);
						vectorC = MatrixOperations.addRow(vectorC, 0);
						numVars++;
						variablesClusterList.add(null);
						
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
	private int addValidInequalities2(Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<Integer, Double>> cMapOriginal,
			Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, Integer>>>> cMapConstraint) {
		log.info("Find valid inequalities of type 2");
		int numCuts = 0;

		CuttingPlanes<Integer> cuttingPlanes = new CuttingPlanes<>(cMapOriginal, cMapConstraint);
		cuttingPlanes.setMaxNumCuts(GlobalParam.CUTTING_PLANES_NUM_CUTS);
		List<List<Integer>> cuts = cuttingPlanes.findValidInequalities2();
		Set<Set<Integer>> setCuts = new LinkedHashSet<>();
		for(List<Integer> cut: cuts) {
			setCuts.add(new LinkedHashSet<>(cut));
		}
		for(Set<Integer> cut: setCuts) {
			numCuts++;
			
			// Add a constraint
			RealMatrix rowA = new Array2DRowRealMatrix(numSlackVars+numVars, 1);
			for(int idx: cut) {
				rowA.setEntry(idx, 0, 1);
			}
			matrixA = MatrixOperations.addRow(matrixA, rowA);
			vectorB = MatrixOperations.addRow(vectorB, 1+GlobalParam.EPSILON5);
			numCons++;
			
			// Add a variable
			RealMatrix colA = new Array2DRowRealMatrix(numCons, 1);
			colA.setEntry(numCons-1, 0, 1);
			matrixA = MatrixOperations.addColumn(matrixA, colA);
			vectorC = MatrixOperations.addRow(vectorC, 0);
			numVars++;
			
			variablesClusterList.add(null);
		}
		
		return numCuts;
	}

	private PrimalDualPoint updateWarmStart(PrimalDualPoint curPoint) {
		RealMatrix newX = new Array2DRowRealMatrix(numSlackVars + numVars, 1);
		MatrixOperations.setElementsMatrix(newX, curPoint.getPointX(), 0.01);

		RealMatrix newS = new Array2DRowRealMatrix(numSlackVars + numVars, 1);
		MatrixOperations.setElementsMatrix(newS, curPoint.getPointS(), 0.01);

		return new PrimalDualPoint(newX, curPoint.getPointPi(), newS);
	}

	private double getDualLowerBound(PrimalDualPoint curPoint) {
		return vectorB.transpose().multiply(curPoint.getPointPi()).getEntry(0, 0);
	}

	private DualInformation getDuals(PrimalDualPoint curPoint) {
		RealMatrix vectorPi = curPoint.getPointPi();

		double[][] dualsPointTemp = new double[instance.getNumClusters()][instance.getNumPoints()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			for(int j = 0; j < instance.getNumPoints(); j++) {
				dualsPointTemp[i][j] = vectorPi.getEntry(i*instance.getNumPoints() + j, 0);
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

		int startIndex = instance.getNumClusters() * instance.getNumPoints();
		double[] dualsMergeTemp = new double[instance.getNumClusters()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			dualsMergeTemp[i] = vectorPi.getEntry(startIndex+i, 0);
		}
		startIndex += instance.getNumClusters();
		double[] dualsStartTemp = new double[instance.getNumClusters()];
		for(int i = 0; i < instance.getNumClusters(); i++) {
			dualsStartTemp[i] = vectorPi.getEntry(startIndex+i, 0); // These duals have a different sign
		}

		double[] dualsMerge = new double[instance.getStartEndLevelList().size()];
		for(int i = 0; i < instance.getStartEndLevelList().size(); i++) {
			Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(i);
			dualsMerge[i] = dualsMergeTemp[pair.second] + dualsStartTemp[pair.first];
		}
		return new DualInformation(dualsPoint, dualsSqrtPoint, dualsMerge, null, null);
	}
	
	private Solution getCurrentSolution(PrimalDualPoint curPoint) {
		Map<Integer, List<Cluster>> hierarchicalClusters = new LinkedHashMap<>();
		for(int i = 0; i < instance.getNumClusters(); i++) {
			hierarchicalClusters.put(i, new ArrayList<>());
		}
		for(int i = numSlackVars; i < curPoint.getPointX().getRowDimension(); i++) {
			double val = curPoint.getPointX().getEntry(i, 0);
			if(val > 0.5) { 
				int actualIdx = i - numSlackVars;
				PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
				if(clusterInfo==null) {
					continue;
				}
				for(int level = clusterInfo.getLevel().second; level <= clusterInfo.getLevel().first; level++) {
					hierarchicalClusters.get(level).add(clusterInfo.getCluster());
				}
			}
		}
		Solution solution = new Solution(instance, hierarchicalClusters, curObj);
		return solution;
	}

	private Solution getCurrentIntegerSolution(PrimalDualPoint curPoint) {
		Map<Integer, List<Cluster>> hierarchicalClusters = new LinkedHashMap<>();
		for(int i = 0; i < instance.getNumClusters(); i++) {
			hierarchicalClusters.put(i, new ArrayList<>());
		}
		for(int i = numSlackVars; i < curPoint.getPointX().getRowDimension(); i++) {
			double val = curPoint.getPointX().getEntry(i, 0);
			if(val > GlobalParam.EPSILON1) { 
				int actualIdx = i - numSlackVars;
				PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
				if(clusterInfo==null) {
					continue;
				}
				for(int level = clusterInfo.getLevel().second; level <= clusterInfo.getLevel().first; level++) {
					hierarchicalClusters.get(level).add(clusterInfo.getCluster());
				}
			}
		}
		Solution solution = new Solution(instance, hierarchicalClusters, instance.computeHierarchicalClusterCost(hierarchicalClusters));
		solution.setFeasible(true);
		return solution;
	}

	private Pair<Map<Integer, List<Cluster>>, Double> getClustersFromIPM(PrimalDualPoint curPoint) {
		Map<Integer, List<Cluster>> hierarchicalClusters = new LinkedHashMap<>();
		double obj = 0;
		for(int i = numSlackVars; i < curPoint.getPointX().getRowDimension(); i++) {
			double val = curPoint.getPointX().getEntry(i, 0);
			if(val > GlobalParam.EPSILON10) { // We want a high precision when calculating UB
				int actualIdx = i - numSlackVars;
				PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
				if(clusterInfo==null) {
					continue;
				}
				double dist = ObjectiveWeightFunction.getWeight(clusterInfo.getLevel().first, clusterInfo.getLevel().second);
				obj += dist * instance.computeClusterCost(clusterInfo.getCluster()) * curPoint.getPointX().getEntry(i, 0);
				for(int level = clusterInfo.getLevel().second; level <= clusterInfo.getLevel().first; level++) {
					if(!hierarchicalClusters.containsKey(level)) {
						hierarchicalClusters.put(level, new ArrayList<>());
					}
					hierarchicalClusters.get(level).add(clusterInfo.getCluster());
				}
			}
		}
		return new Pair<>(hierarchicalClusters, obj);
	}
	
	private void printClustersFromIPM(PrimalDualPoint curPoint) {
		for(int i = numSlackVars; i < curPoint.getPointX().getRowDimension(); i++) {
			double val = curPoint.getPointX().getEntry(i, 0);
			if(val > GlobalParam.EPSILON2) { 
				int actualIdx = i - numSlackVars;
				PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
				if(clusterInfo==null) {
					continue;
				}
				System.out.println(UtilFunctions.roundDecimals(val, 5) + " "+clusterInfo.getLevel() + " "+ clusterInfo.getCluster());
			}
		}
	}
	
	/**
	 * Return which level it was k=1,..., K
	 * @param curPoint
	 * @return
	 */
	private Map<Pair<Integer, Cluster>, Double> getClusterNonIntegralityPerLevel(PrimalDualPoint curPoint) {
		Map<Pair<Integer, Cluster>, Double> result = new LinkedHashMap<>();
		for(int i = numSlackVars; i < curPoint.getPointX().getRowDimension(); i++) {
			double val = curPoint.getPointX().getEntry(i, 0);
			int actualIdx = i - numSlackVars;
			PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
			if(clusterInfo==null) {
				continue;
			}
			for(int level = clusterInfo.getLevel().second; level <= clusterInfo.getLevel().first; level++) {
				Pair<Integer, Cluster> cluster = new Pair<>(level, clusterInfo.getCluster());
				if(!result.containsKey(cluster)) {
					result.put(cluster, 0d);
				}
				result.replace(cluster, result.get(cluster)+val);
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
	 * Return the index in numVars
	 * @param curPoint
	 * @return
	 */
	private Map<Pair<Integer, PDCGM_ClusterInfo>, Double> getClusterNonIntegrality(PrimalDualPoint curPoint) {
		Map<Pair<Integer, PDCGM_ClusterInfo>, Double> result = new LinkedHashMap<>();
		for(int i = numSlackVars; i < curPoint.getPointX().getRowDimension(); i++) {
			double val = curPoint.getPointX().getEntry(i, 0);
			if(isFractional(val)) {
				int actualIdx = i - numSlackVars;
				PDCGM_ClusterInfo clusterInfo = variablesClusterList.get(actualIdx);
				if(clusterInfo==null) {
					continue;
				}
				result.put(new Pair<>(i, clusterInfo), val);
			}
		}
		return result;
	}

	private boolean isFractional(double val) {
		double rnd = Math.round(val);
		if (Math.abs(rnd - val) > 5*GlobalParam.EPSILON2) {
			return true;
		}
		return false;
	}

	private void initFirstCluster() {
		boolean[] pointId = new boolean[instance.getNumPoints()];
		Arrays.fill(pointId, Boolean.TRUE);

		Cluster cluster = new Cluster(pointId);
		addCluster(new Pair<>(0, 0), cluster, true);


	}

	private void initEmptyLists(boolean resetGlobalCounters) {
		variablesClusterMap = new LinkedHashMap<>();
		for(Pair<Integer, Integer> pair: instance.getStartEndLevelList()) {
			variablesClusterMap.put(pair, new LinkedHashMap<>());
		}
		variablesClusterList = new ArrayList<>();
		startSolutionObjectives = new ArrayList<>();
		
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

	private void initPricingProblems() {
		if(GlobalParam.DISTANCE_TYPE==DistanceType.Centroid) {
			pp = new PricingProblemColgen(instance, pricingProblemStrategy);
		}
		else {
			pp = new PricingProblemColgenMedoid(instance);
		}
		if(solverType==SolverType.Colgen_CuttingPlanes) {
			pp.setColumnLimit(GlobalParam.PP_NUM_POTENTIAL_COLUMNS);
		}
	}

	public void initHierarchicalStartSolutions(List<Solution> startSolutions) {
		this.startSolutions = startSolutions;
		log.info("Adding {} startsolutions", startSolutions.size());
		for(Solution startSolution : startSolutions) {
			initHierarchicalStartSolution(startSolution);
		}
	}

	public void initHierarchicalStartSolution(Solution startSolution) {
		if(GlobalParam.OBJECTIVE_WEIGHT==ObjectiveWeight.proportionalToLevel) {
			ObjectiveWeightFunction.setWeightsProportionalToLevel(instance, startSolution.getHierarchicalClusters());
			startSolution = new Solution(instance, startSolution.getHierarchicalClusters(), instance.computeHierarchicalClusterCost(startSolution.getHierarchicalClusters()));
		}
		
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

	private boolean addClusterToEachLevel(Cluster cluster, boolean startSolution) {
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

	/**
	 * If dummyCluster, ignore all checks
	 * @param pair
	 * @param cluster
	 * @param startSolution
	 * @param dummyCluster
	 * @return
	 */
	private boolean addCluster(Pair<Integer, Integer> pair, Cluster cluster, boolean startSolution) {
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
		if(variablesClusterMap.get(pair).containsKey(cluster)) {
			if (!startSolution && GlobalParam.DEBUG) {
				throw new IllegalArgumentException("Cluster "+ cluster + " was already added to level "+ pair);
			}
			return false;
		}

		RealMatrix colA = new Array2DRowRealMatrix(numCons, 1);
		for(int level = pair.second; level <= pair.first; level++) {
			// Cover constraints
			int idx = level*instance.getNumPoints();
			for(int i: cluster.getPointId()) {
				colA.setEntry(idx+i, 0, 1);
			}
		}
		// Merge constraints
		colA.setEntry(numCoverConstraints + pair.second, 0, 1);
		// Start constraints
		colA.setEntry(numCoverConstraints + instance.getNumClusters() + pair.first, 0, 1);
		double dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
		double clusterCost = instance.computeClusterCost(cluster) * dist;

		if(matrixA==null) {
			matrixA = colA;
			vectorC = new Array2DRowRealMatrix(1, 1);
			vectorC.setEntry(0, 0, clusterCost);
		}
		else {
			matrixA = MatrixOperations.addColumn(matrixA, colA);
			vectorC = MatrixOperations.addRow(vectorC, clusterCost);
		}
		PDCGM_ClusterInfo clusterInfo = new PDCGM_ClusterInfo(pair, cluster, clusterCost, numVars);
		numVars++;

		variablesClusterMap.get(pair).put(cluster, clusterInfo);
		variablesClusterList.add(clusterInfo);
		return true;
	}
	
	public Solution getBestSolution(boolean setDetailedInformation) {
		return getBestSolution(setDetailedInformation, totalTime);
	}

	public Solution getBestSolution(boolean setDetailedInformation, long totalTime) {
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
}
