package cl.data;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonProperty;

import cl.data.type.DistanceType;
import cl.util.Pair;

public class Instance {
	private String instanceName;
	private double[][] coordinates, coordinatesOriginal;
	private int numPoints;
	private int numDimensions;
	private int numClusters;
	private double[][] squaredDistanceMatrix, regularDistanceMatrix;	
	
	// True cluster
	private int trueNumClusters = 0;
	private int[] trueClusterLabels;
	
	private List<Pair<Integer, Integer>> startEndLevelList;

	@JsonCreator
	public Instance(@JsonProperty("instanceName") String instanceName, @JsonProperty("coordinates") double[][] coordinates, @JsonProperty("coordinatesOriginal") double[][] coordinatesOriginal,
			@JsonProperty("numPoints") int numPoints,
			@JsonProperty("numDimensions") int numDimensions, @JsonProperty("numClusters") int numClusters, @JsonProperty("squaredEuclideanDistanceMatrix") double[][] squaredEuclideanDistanceMatrix,
			@JsonProperty("regularDistanceMatrix") double[][] regularDistanceMatrix, @JsonProperty("trueNumClusters") int trueNumClusters,
			@JsonProperty("trueClusterLabels") int[] trueClusterLabels) {
		this.instanceName = instanceName;
		this.coordinates = coordinates;
		this.coordinatesOriginal = coordinatesOriginal;
		this.numPoints = numPoints;
		this.numDimensions = numDimensions;
		this.numClusters = numClusters;
		this.squaredDistanceMatrix = squaredEuclideanDistanceMatrix;
		this.regularDistanceMatrix = regularDistanceMatrix;
		this.trueNumClusters = trueNumClusters;
		this.trueClusterLabels = trueClusterLabels;
	}

	public Instance(String instanceName, int trueNumClusters, int[] trueClusterLabels, double[][] coordinates, double[][] coordinatesOriginal, int numPoints, int numDimensions, int numClusters) {
		this.instanceName = instanceName;
		this.trueNumClusters = trueNumClusters;
		this.trueClusterLabels = trueClusterLabels;
		this.coordinates = coordinates;
		this.coordinatesOriginal = coordinatesOriginal;
		this.numPoints = numPoints;
		this.numDimensions = numDimensions;
		this.numClusters = numClusters;
	}

	public Instance(String instanceName, int numPoints, int numDimensions, int numClusters,
			double[][] squaredEuclideanDistanceMatrix, double[][] regularDistanceMatrix) {
		this.instanceName = instanceName;
		this.numPoints = numPoints;
		this.numDimensions = numDimensions;
		this.numClusters = numClusters;
		this.squaredDistanceMatrix = squaredEuclideanDistanceMatrix;
		this.regularDistanceMatrix = regularDistanceMatrix;
	}

	@JsonIgnore
	public List<Pair<Integer, Integer>> getStartEndLevelList() {
		return startEndLevelList;
	}

	public void initStartEndLevelList() {
		startEndLevelList = new ArrayList<>();
		for(int startLevel = 0; startLevel < numClusters; startLevel++) {
			for(int endLevel = startLevel; endLevel >= 0; endLevel--) {
				startEndLevelList.add(new Pair<>(startLevel, endLevel));
			}
		}
	}

	public int getTrueNumClusters() {
		return trueNumClusters;
	}

	public int[] getTrueClusterLabels() {
		return trueClusterLabels;
	}
	
	@JsonIgnore
	public List<Cluster> getTrueClusters() {
		List<Cluster> trueClusters = new ArrayList<>();
		Map<Integer, List<Integer>> temp = new LinkedHashMap<>();
		for(int i = 0; i < numPoints; i++) {
			int clusterID = trueClusterLabels[i];
			
			if(!temp.containsKey(clusterID)) {
				temp.put(clusterID, new ArrayList<>());
			}
			temp.get(clusterID).add(i);
		}
		
		for(Entry<Integer, List<Integer>> entry: temp.entrySet()) {
			trueClusters.add(new Cluster(entry.getValue().stream().mapToInt(j -> j).toArray(), numPoints));
		}
		return trueClusters;
	}

	public String getInstanceName() {
		return instanceName;
	}

	public int getNumPoints() {
		return numPoints;
	}

	public double[][] getCoordinates() {
		return coordinates;
	}

	public int getNumDimensions() {
		return numDimensions;
	}

	public int getNumClusters() {
		return numClusters;
	}

	public void setNumClusters(int numClusters) {
		this.numClusters = numClusters;
	}

	public double[][] getSquaredDistanceMatrix() {
		if(squaredDistanceMatrix==null) {
			throw new IllegalArgumentException("Distance matrix is not initialised, use computeSquaredEuclideanDistanceMatrix.");
		}
		return squaredDistanceMatrix;
	}
	
	/**
	 * Coordinates before standardisation
	 * @return
	 */
	@JsonIgnore
	public Instance getOriginalInstance() {
		Instance instance = new Instance(instanceName, trueNumClusters, trueClusterLabels, coordinatesOriginal, coordinatesOriginal, numPoints, numDimensions, numClusters);
		instance.computeDistanceMatrices();
		return instance;
	}

	/**
	 * Symmetric matrix containing d_ij^2
	 * @param numPoints
	 * @param coordinates
	 * @return
	 */
	public void computeDistanceMatrices() {
		squaredDistanceMatrix = new double[numPoints][numPoints];
		regularDistanceMatrix = new double[numPoints][numPoints]; 
		for(int i = 0; i < numPoints; i++) {
			for(int j = i; j < numPoints; j++) {
				if(i==j) {
					squaredDistanceMatrix[i][j] = 0;
					regularDistanceMatrix[i][j] = 0;
				}
				else {
					double squaredDist = GlobalParam.DISTANCE_FUNCTION.getSquaredDistance(coordinates[i], coordinates[j]); //EuclideanDistance.getSquaredEuclideanDistance(coordinates[i], coordinates[j]);
					squaredDistanceMatrix[i][j] = squaredDist;
					squaredDistanceMatrix[j][i] = squaredDist;
					
					double regularDist = GlobalParam.DISTANCE_FUNCTION.getRegularDistance(coordinates[i], coordinates[j]); // EuclideanDistance.getRegularEuclideanDistance(coordinates[i], coordinates[j]);
					regularDistanceMatrix[i][j] = regularDist;
					regularDistanceMatrix[j][i] = regularDist;
				}
			}
		}
	}

	public double[][] getRegularDistanceMatrix() {
		return regularDistanceMatrix;
	}

	/**
	 * Sum of squared distances from all points of a given set to its centroid
	 * Is equal to the sum of squared distances between pairs of points of this set, divided by its cardinality
	 * @param cluster
	 * @return
	 */
	public double computeClusterCost(Cluster cluster) {
		if(GlobalParam.DISTANCE_TYPE==DistanceType.Centroid) {
			double cost = 0;
			for(int i: cluster.getPointId()) {
				for(int j: cluster.getPointId()) {
					if(i>=j) {
						continue;
					}
					cost += squaredDistanceMatrix[i][j]; 
				}
			}
			cost = cost / ((double) cluster.getPointId().length);
			return cost;
		}
		else {
			double cost = Double.POSITIVE_INFINITY;
			List<Integer> poindId = Arrays.stream(cluster.getPointId()).boxed().collect(Collectors.toList());
			Medoid medoid = new Medoid(poindId, null, -1);
			for(Integer i: cluster.getPointId()) { // Since we do not know what the medoid is at the moment
				medoid.setMedoidCoord(coordinates[i]);
				medoid.setMedoidId(i);
				double candidateCost = GlobalParam.DISTANCE_FUNCTION.getDistanceOfMedoid(this, medoid);
				if(candidateCost<cost) {
					cost = candidateCost;
				}
			}
			return cost;
		}
	}
	
	/**
	 * We assume the cluster at level 1 is also present
	 * @param hierarchicalClusters
	 * @return
	 */
	public double computeHierarchicalClusterCost(Map<Integer, List<Cluster>> hierarchicalClusters) {
		double cost = 0;
		if(GlobalParam.OBJECTIVE_WEIGHT==ObjectiveWeight.noWeight) {
			for(Entry<Integer, List<Cluster>> entry: hierarchicalClusters.entrySet()) {
				for(Cluster cluster: entry.getValue()) {
					cost += computeClusterCost(cluster);
				}
			}
			return cost;
		}
		else if(GlobalParam.OBJECTIVE_WEIGHT==ObjectiveWeight.noDoubleCounting) {
			return hierarchicalClusters.values().stream()
		    .flatMap(List::stream)
		    .distinct()
		    .mapToDouble(this::computeClusterCost)
		    .sum();
		}		
		else if(GlobalParam.OBJECTIVE_WEIGHT==ObjectiveWeight.proportionalToLevel) {
			for(Entry<Integer, List<Cluster>> entry: hierarchicalClusters.entrySet()) {
				for(Cluster cluster: entry.getValue()) {
					cost += ObjectiveWeightFunction.getProportionalToLevel(entry.getKey()) * computeClusterCost(cluster); 
				}
			}
			return cost;
		}
		else {
			throw new IllegalArgumentException("OBJECTIVE_WEIGHT not defined");
		}
	}

	public double computeFirstClusterCost() {
		int[] pointId = IntStream.range(0, numPoints).toArray();
		Cluster cluster = new Cluster(pointId, numPoints);
		return computeClusterCost(cluster);
	}
	
	public static Instance readFromTxt(File f, int numClusters) {
		return readFromTxt(f, numClusters, true, true);
	}
	
	public static Instance readFromTxt(File f, int numClusters, boolean standardise, boolean computeDistanceMatrix) {
		try {
			Scanner sc = new Scanner(f);

			String header = sc.nextLine();
			String[] headers = header.split(" ");
			int numPoints = Integer.parseInt(headers[0]);
			int numDimensions = Integer.parseInt(headers[1]);
			
			boolean trueClusters = false;
			int trueNumClusters = 0;
			int[] trueClusterLabels = null;
			if(headers.length==3) {
				trueClusters = true;
				trueNumClusters = Integer.parseInt(headers[2]);
				trueClusterLabels = new int[numPoints];
			}
			
			if(numClusters>numPoints) {
				sc.close();
				throw new IllegalArgumentException("NumClusters should be less than or equal to numPoints");
			}

			double[][] coordinatesOriginal = new double[numPoints][numDimensions];
			for(int i = 0; i < numPoints; i++) {
				for(int j = 0; j < numDimensions; j++) {
					String s = sc.next();
					coordinatesOriginal[i][j] = Double.parseDouble(s);
				}
				if(trueClusters) {
					trueClusterLabels[i] = sc.nextInt();
				}
			}
			double[][] coordinates = new double[numPoints][numDimensions];
			if(standardise) {
				for(int j = 0; j < numDimensions; j++) {
					// Calculate mean and std
					SummaryStatistics ss = new SummaryStatistics();
					for(int i = 0; i < numPoints; i++) {
						ss.addValue(coordinatesOriginal[i][j]);
					}
					
					// New coordinates
					for(int i = 0; i < numPoints; i++) {
						coordinates[i][j] = (coordinatesOriginal[i][j] - ss.getMean()) / ss.getStandardDeviation();
					}
				}
			}
			else {
				for(int j = 0; j < numDimensions; j++) {
					for(int i = 0; i < numPoints; i++) {
						coordinates[i][j] = coordinatesOriginal[i][j];
					}
				}
			}

			sc.close();
			Instance instance = new Instance(f.getName(), trueNumClusters, trueClusterLabels, coordinates, coordinatesOriginal, numPoints, numDimensions, numClusters);
			if(computeDistanceMatrix) {
				// TODO: Do we need to calculate this more efficiently?
				instance.computeDistanceMatrices();
			}
			return instance;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return null;
	}


}