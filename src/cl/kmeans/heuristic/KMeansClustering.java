package cl.kmeans.heuristic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import cl.data.GlobalParam;
import cl.data.Instance;
import cl.util.distance.EuclideanDistance;

public class KMeansClustering {

	private Instance instance;

	public KMeansClustering(Instance instance) {
		this.instance = instance;
	}

	public Map<CentroidKMeans, List<Integer>> fit() { 
		int maxIterations = 1000;
		List<CentroidKMeans> centroids = randomCentroids();
		Map<CentroidKMeans, List<Integer>> clusters = new HashMap<>();
		Map<CentroidKMeans, List<Integer>> lastState = new HashMap<>();

		// iterate for a pre-defined number of times
		for (int i = 0; i < maxIterations; i++) {
			boolean isLastIteration = i == maxIterations - 1;
			// in each iteration we should find the nearest centroid for each record
			for (int j = 0; j < instance.getNumPoints(); j++) {
				CentroidKMeans centroid = nearestCentroid(j, centroids);
				assignToCluster(clusters, j, centroid);
			}

			// if the assignments do not change, then the algorithm terminates
			boolean shouldTerminate = isLastIteration || clusters.equals(lastState);
			lastState = clusters;
			if (shouldTerminate) { 
				break; 
			}
			// at the end of each iteration we should relocate the centroids
			centroids = relocateCentroids(clusters);
			clusters = new HashMap<>();
		}

		return lastState;
	}

	private CentroidKMeans Average(CentroidKMeans centroid, List<Integer> records) {
		if (records == null || records.isEmpty()) { 
			return centroid;
		}

		double[] average = centroid.getCoordinates();
		for (int record : records) {
			for(int j = 0; j < instance.getNumDimensions(); j++) {
				average[j] += instance.getCoordinates()[record][j];
			}
		}
		for(int j = 0; j < instance.getNumDimensions(); j++) {
			average[j] = average[j] / ((double) records.size());
		}
		return new CentroidKMeans(average);
	}

	private List<CentroidKMeans> relocateCentroids(Map<CentroidKMeans, List<Integer>> clusters) {
		return clusters.entrySet().stream().map(e -> Average(e.getKey(), e.getValue())).collect(Collectors.toList());
	}

	private List<CentroidKMeans> randomCentroids() {
		List<CentroidKMeans> centroids = new ArrayList<>();
		double[] maxs = new double[instance.getNumDimensions()];
		double[] mins = new double[instance.getNumDimensions()];

		for (int i = 0; i < instance.getNumPoints(); i++) {
			for(int j = 0; j < instance.getNumDimensions(); j++) {
				if(j==0 || instance.getCoordinates()[i][j] > maxs[j]) {
					maxs[j] = instance.getCoordinates()[i][j];
				}
				if(j== 0 || instance.getCoordinates()[i][j] < mins[j]) {
					mins[j] = instance.getCoordinates()[i][j];
				}
			}
		}
		for (int i = 0; i < instance.getNumClusters(); i++) {
			double[] coordinates = new double[instance.getNumDimensions()];
			for (int j = 0; j < instance.getNumDimensions(); j++) {
				double max = maxs[j];
				double min = mins[j];
				coordinates[j] = GlobalParam.RANDOM.nextDouble() * (max - min) + min;
			}
			centroids.add(new CentroidKMeans(coordinates));
		}
		return centroids;
	}

	private void assignToCluster(Map<CentroidKMeans, List<Integer>> clusters, int i, CentroidKMeans centroid) {
		clusters.compute(centroid, (key, list) -> {
			if (list == null) {
				list = new ArrayList<>();
			}

			list.add(i);
			return list;
		});
	}

	private CentroidKMeans nearestCentroid(int i, List<CentroidKMeans> centroids) {
		double minimumDistance = Double.MAX_VALUE;
		CentroidKMeans nearest = null;

		for (CentroidKMeans centroid : centroids) {
			double currentDistance = EuclideanDistance.getSquaredEuclideanDistance(instance.getCoordinates()[i], centroid.getCoordinates());

			if (currentDistance < minimumDistance) {
				minimumDistance = currentDistance;
				nearest = centroid;
			}
		}

		return nearest;
	}

	public static class CentroidKMeans {

		private final double[] coordinates;

		public CentroidKMeans(double[] coordinates) {
			this.coordinates = coordinates;
		}

		public double[] getCoordinates() {
			return coordinates;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + Arrays.hashCode(coordinates);
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			CentroidKMeans other = (CentroidKMeans) obj;
			if (!Arrays.equals(coordinates, other.coordinates))
				return false;
			return true;
		}


	}

	public static class Record {
		private final String description;
		private final Map<String, Double> features;

		public Record(String description, Map<String, Double> features) {
			this.description = description;
			this.features = features;
		}

		public String getDescription() {
			return description;
		}

		public Map<String, Double> getFeatures() {
			return features;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((description == null) ? 0 : description.hashCode());
			result = prime * result + ((features == null) ? 0 : features.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Record other = (Record) obj;
			if (description == null) {
				if (other.description != null)
					return false;
			} else if (!description.equals(other.description))
				return false;
			if (features == null) {
				if (other.features != null)
					return false;
			} else if (!features.equals(other.features))
				return false;
			return true;
		}
	}
}
