package cl.hierarchical.model.extended;

import java.util.Map;
import java.util.Set;

import cl.data.Cluster;

public interface PricingProblem {
	void setDuals(DualInformation duals);
	void setBranchInformation(BranchInformation bi);
	Map<Integer, Set<Cluster>> generateColumns();
	Map<Integer, Set<Cluster>> generateColumns(boolean enumerate, Map<Integer, Double> bounds);
	Map<Integer, Map<Cluster, Double>> getRC();
	double getBestRC();
	void setColumnLimit(int columnLimit);
}
