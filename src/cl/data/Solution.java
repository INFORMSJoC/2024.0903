package cl.data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonProperty;

import cl.data.type.HierarchicalHeuristicType;
import cl.data.type.PricingProblemStrategy;
import cl.data.type.SolverType;
import cl.util.AdjustedRandIndex;

public class Solution {
	private Instance instance;
	private double objective;
	private long totalTime;
	private boolean feasible; //If the solution is integer valued
	private List<Cluster> clusters;
	private Map<Integer, List<Cluster>> hierarchicalClusters;
	private SolverType solverType;
	private PricingProblemStrategy pricingProblemStrategy;
	private HierarchicalHeuristicType hierarchicalHeuristicType;
	private boolean timeLimitReached;
	
	// Column generation
	private int numIter, numColumns;
	private long heuristicTime, rootNodeTime, pricingTime, masterTime;
	private double bestStartObjective, bestLP;
	private List<Double> startObjectives;
	private int numBranchAndPrice, numVNS, numQP, numBranchAndBound, numKnapsack, numCuttingPlanes, numEnumerate;
	
	@JsonCreator
	public Solution(@JsonProperty("instance") Instance instance, @JsonProperty("objective") double objective, @JsonProperty("totalTime") long totalTime,
			@JsonProperty("feasible") boolean feasible, @JsonProperty("clusters") List<Cluster> clusters,
			@JsonProperty("hierarchicalClustersList") List<List<Cluster>> hierarchicalClustersList, @JsonProperty("solverType") SolverType solverType,
			@JsonProperty("pricingProblemStrategy") PricingProblemStrategy pricingProblemStrategy, @JsonProperty("hierarchicalHeuristicType") HierarchicalHeuristicType hierarchicalHeuristicType,
			@JsonProperty("timeLimitReached") boolean timeLimitReached, @JsonProperty("numIter") int numIter, @JsonProperty("numColumns") int numColumns,
			@JsonProperty("heuristicTime") long heuristicTime, @JsonProperty("rootNodeTime") long rootNodeTime,
			@JsonProperty("pricingTime") long pricingTime, @JsonProperty("masterTime") long masterTime, @JsonProperty("bestStartObjective") double bestStartObjective,
			@JsonProperty("bestLP") double bestLP, @JsonProperty("startObjectives") List<Double> startObjectives,
			@JsonProperty("numBranchAndPrice") int numBranchAndPrice, @JsonProperty("numVNS") int numVNS, @JsonProperty("numQP") int numQP,
			@JsonProperty("numBranchAndBound") int numBranchAndBound, @JsonProperty("numKnapsack") int numKnapsack, @JsonProperty("numCuttingPlanes") int numCuttingPlanes,
			@JsonProperty("numEnumerate") int numEnumerate) {
		this.instance = instance;
		this.objective = objective;
		this.totalTime = totalTime;
		this.feasible = feasible;
		this.clusters = clusters;
		this.hierarchicalClusters = new LinkedHashMap<>(); 
		for(int i = 0; i < hierarchicalClustersList.size(); i++) {
			this.hierarchicalClusters.put(i, hierarchicalClustersList.get(i));
		}
		this.solverType = solverType;
		this.pricingProblemStrategy = pricingProblemStrategy;
		this.hierarchicalHeuristicType = hierarchicalHeuristicType;
		this.timeLimitReached = timeLimitReached;
		this.numIter = numIter;
		this.numColumns = numColumns;
		this.heuristicTime = heuristicTime;
		this.rootNodeTime = rootNodeTime;
		this.pricingTime = pricingTime;
		this.masterTime = masterTime;
		this.bestStartObjective = bestStartObjective;
		this.bestLP = bestLP;
		this.startObjectives = startObjectives;
		this.numBranchAndPrice = numBranchAndPrice;
		this.numVNS = numVNS;
		this.numQP = numQP;
		this.numBranchAndBound = numBranchAndBound;
		this.numKnapsack = numKnapsack;
		this.numCuttingPlanes = numCuttingPlanes;
		this.numEnumerate = numEnumerate;
	}

	public Solution(Instance instance, double objective) {
		this.instance = instance;
		this.objective = objective;
	}
	
	public Solution(Instance instance, List<Cluster> clusters, double objective) {
		this.instance = instance;
		this.clusters = clusters;
		this.objective = objective;
	}
	
	public Solution(Instance instance, Map<Integer, List<Cluster>> hierarchicalClusters, double objective) {
		this.instance = instance;
		this.hierarchicalClusters = hierarchicalClusters;
		this.objective = objective;
	}

	public HierarchicalHeuristicType getHierarchicalHeuristicType() {
		return hierarchicalHeuristicType;
	}

	public void setHierarchicalHeuristicType(HierarchicalHeuristicType hierarchicalHeuristicType) {
		this.hierarchicalHeuristicType = hierarchicalHeuristicType;
	}

	public PricingProblemStrategy getPricingProblemStrategy() {
		return pricingProblemStrategy;
	}

	public void setPricingProblemStrategy(PricingProblemStrategy pricingProblemStrategy) {
		this.pricingProblemStrategy = pricingProblemStrategy;
	}

	public boolean isTimeLimitReached() {
		return timeLimitReached;
	}

	public void setTimeLimitReached(boolean timeLimitReached) {
		this.timeLimitReached = timeLimitReached;
	}

	public int getNumEnumerate() {
		return numEnumerate;
	}

	public void setNumEnumerate(int numEnumerate) {
		this.numEnumerate = numEnumerate;
	}

	public int getNumKnapsack() {
		return numKnapsack;
	}

	public void setNumKnapsack(int numKnapsack) {
		this.numKnapsack = numKnapsack;
	}

	public int getNumBranchAndPrice() {
		return numBranchAndPrice;
	}

	public void setNumBranchAndPrice(int numBranchAndPrice) {
		this.numBranchAndPrice = numBranchAndPrice;
	}

	public int getNumVNS() {
		return numVNS;
	}

	public void setNumVNS(int numVNS) {
		this.numVNS = numVNS;
	}

	public int getNumQP() {
		return numQP;
	}

	public void setNumQP(int numQP) {
		this.numQP = numQP;
	}

	public int getNumBranchAndBound() {
		return numBranchAndBound;
	}

	public void setNumBranchAndBound(int numBranchAndBound) {
		this.numBranchAndBound = numBranchAndBound;
	}

	public int getNumCuttingPlanes() {
		return numCuttingPlanes;
	}

	public void setNumCuttingPlanes(int numCuttingPlanes) {
		this.numCuttingPlanes = numCuttingPlanes;
	}

	public SolverType getSolverType() {
		return solverType;
	}

	public void setSolverType(SolverType solverType) {
		this.solverType = solverType;
	}

	public Instance getInstance() {
		return instance;
	}

	@JsonIgnore
	public Map<Integer, List<Cluster>> getHierarchicalClusters() {
		return hierarchicalClusters;
	}
	
	public List<List<Cluster>> getHierarchicalClustersList() {
		List<List<Cluster>> hierarchicalClustersList = new ArrayList<>();
		for(Entry<Integer, List<Cluster>> entry: hierarchicalClusters.entrySet()) {
			hierarchicalClustersList.add(entry.getValue());
		}
		return hierarchicalClustersList;
	}

	public List<Cluster> getClusters() {
		return clusters;
	}

	public double getObjective() {
		return objective;
	}

	public long getHeuristicTime() {
		return heuristicTime;
	}

	public void setHeuristicTime(long heuristicTime) {
		this.heuristicTime = heuristicTime;
	}

	public void setTotalTime(long totalTime) {
		this.totalTime = totalTime;
	}

	public long getTotalTime() {
		return totalTime;
	}

	public boolean isFeasible() {
		return feasible;
	}

	public void setFeasible(boolean feasible) {
		this.feasible = feasible;
	}

	public int getNumIter() {
		return numIter;
	}

	public void setNumIter(int numIter) {
		this.numIter = numIter;
	}

	public long getRootNodeTime() {
		return rootNodeTime;
	}

	public void setRootNodeTime(long rootNodeTime) {
		this.rootNodeTime = rootNodeTime;
	}

	public long getPricingTime() {
		return pricingTime;
	}

	public void setPricingTime(long pricingTime) {
		this.pricingTime = pricingTime;
	}

	public long getMasterTime() {
		return masterTime;
	}

	public void setMasterTime(long masterTime) {
		this.masterTime = masterTime;
	}

	public int getNumColumns() {
		return numColumns;
	}

	public void setNumColumns(int numColumns) {
		this.numColumns = numColumns;
	}

	public double getBestStartObjective() {
		return bestStartObjective;
	}

	public void setBestStartObjective(double bestStartObjective) {
		this.bestStartObjective = bestStartObjective;
	}
	
	public List<Double> getStartObjectives() {
		return startObjectives;
	}

	public void setStartObjectives(List<Double> startObjectives) {
		this.startObjectives = startObjectives;
	}

	public double getBestLP() {
		return bestLP;
	}

	public void setBestLP(double bestLP) {
		this.bestLP = bestLP;
	}
	
	/**
	 * Return objective before standardisation
	 * @return
	 */
	@JsonIgnore
	public double getOriginalObjective() {
		if(hierarchicalClusters==null) {
			return 0;
		}
		return instance.getOriginalInstance().computeHierarchicalClusterCost(hierarchicalClusters);
	}

	/**
	 * Return ARI of the largest level available 
	 * @return
	 */
	@JsonIgnore
	public double getARI() {
		if(instance.getTrueNumClusters()!=0) {
			int ind = instance.getNumClusters()-1;
			List<Integer> predCluster = AdjustedRandIndex.convertClusterToPartition(instance, hierarchicalClusters.get(ind)); 
			List<Integer> trueCluster = AdjustedRandIndex.convertClusterToPartition(instance, instance.getTrueClusters());
			return AdjustedRandIndex.calcARI(predCluster, trueCluster);
		}
		return Double.NaN;
	}
	
	public void printHierarchy() {
		System.out.println("Hierarchy, obj "+ objective);
		for(Entry<Integer, List<Cluster>> entry: hierarchicalClusters.entrySet()) {
			System.out.println("Level: "+entry.getKey());
			for(Cluster cluster: entry.getValue()) {
				System.out.println(Arrays.toString(cluster.getPointId())+ " "+ instance.computeClusterCost(cluster));
			}
		}
		
		if(instance.getTrueNumClusters()!=0) {
			System.out.println();
			System.out.println("True clusters (with costs)");
			for(Cluster cluster: instance.getTrueClusters()) {
				System.out.println(Arrays.toString(cluster.getPointId()) + " "+ instance.computeClusterCost(cluster));
			}
			System.out.println("ARI "+ getARI());
		}
	}
}
