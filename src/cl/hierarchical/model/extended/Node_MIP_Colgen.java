package cl.hierarchical.model.extended;

import cl.data.Cluster;
import cl.util.Pair;
import ilog.concert.IloNumVar;

public class Node_MIP_Colgen {
	private Pair<Integer, Integer> level;
	private Cluster cluster;
	private Pair<IloNumVar, Double> var;
	
	
	
	public Node_MIP_Colgen(Pair<Integer, Integer> level, Cluster cluster, Pair<IloNumVar, Double> var) {
		this.level = level;
		this.cluster = cluster;
		this.var = var;
	}
	public Pair<Integer, Integer> getLevel() {
		return level;
	}
	public Cluster getCluster() {
		return cluster;
	}
	
	public Pair<IloNumVar, Double> getVar() {
		return var;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((cluster == null) ? 0 : cluster.hashCode());
		result = prime * result + ((level == null) ? 0 : level.hashCode());
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
		Node_MIP_Colgen other = (Node_MIP_Colgen) obj;
		if (cluster == null) {
			if (other.cluster != null)
				return false;
		} else if (!cluster.equals(other.cluster))
			return false;
		if (level == null) {
			if (other.level != null)
				return false;
		} else if (!level.equals(other.level))
			return false;
		return true;
	}
	
	
}
