package cl.data;

import java.util.ArrayList;
import java.util.List;

public class Centroid {

	private List<Integer> cluster;
	private double[] coord;
	
	public Centroid(Instance instance) {
		this(new ArrayList<>(), new double[instance.getNumDimensions()]);
	}
	
	public Centroid(List<Integer> cluster, double[] coord) {
		this.cluster = cluster;
		this.coord = coord;
	}

	public List<Integer> getCluster() {
		return cluster;
	}

	public double[] getCoord() {
		return coord;
	}

	public void updateCentroid(Instance instance, Integer ind)
	{
		double curClusterSize = cluster.size();
		for(int i = 0; i < instance.getNumDimensions(); i++)
		{
			// TODO: this is only for Euclidean distance
			coord[i] = (coord[i] * curClusterSize + instance.getCoordinates()[ind][i]) / (curClusterSize + 1);
		}
		cluster.add(ind);
	}
	
	public Centroid copy()
	{
		return new Centroid(new ArrayList<>(cluster), coord.clone());
	}
}
