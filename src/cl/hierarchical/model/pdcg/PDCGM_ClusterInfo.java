package cl.hierarchical.model.pdcg;

import cl.data.Cluster;
import cl.util.Pair;

public class PDCGM_ClusterInfo {
	private Pair<Integer, Integer> level;
	private Cluster cluster;
	private double clusterCost;
	private int id;
	
	public PDCGM_ClusterInfo(Pair<Integer, Integer> level, Cluster cluster, double clusterCost, int id) {
		this.level = level;
		this.cluster = cluster;
		this.clusterCost = clusterCost;
		this.id = id;
	}

	public Pair<Integer, Integer> getLevel() {
		return level;
	}

	public Cluster getCluster() {
		return cluster;
	}

	public double getClusterCost() {
		return clusterCost;
	}

	public int getId() {
		return id;
	}
	
	
}
