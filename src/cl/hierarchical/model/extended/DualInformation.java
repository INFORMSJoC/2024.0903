package cl.hierarchical.model.extended;

import java.util.Map;

import cl.data.Cluster;

public class DualInformation {
	private double[][] dualsPoint;
	private double[][] dualsSqrtPoint;
	private double[] dualsMaxCluster;
	private Map<Integer, Map<Cluster, Double>> dualsHierarchy;
	private double[] dualsOneChangePerHierarchy;

	// When we are branching
	private double[] dualsPointSingle, dualsSqrtPointSingle;

	public DualInformation(double[] dualsPointSingle, double[] dualsSqrtPointSingle) {
		this.dualsPointSingle = dualsPointSingle;
		this.dualsSqrtPointSingle = dualsSqrtPointSingle;
	}

	public DualInformation(double[][] dualsPoint, double[][] dualsSqrtPoint, double[] dualsMaxCluster,
			Map<Integer, Map<Cluster, Double>> dualsHierarchy, double[] dualsOneChangePerHierarchy) {
		this.dualsPoint = dualsPoint;
		this.dualsSqrtPoint = dualsSqrtPoint;
		this.dualsMaxCluster = dualsMaxCluster;
		this.dualsHierarchy = dualsHierarchy;
		this.dualsOneChangePerHierarchy = dualsOneChangePerHierarchy;
	}

	public double[][] getDualsSqrtPoint() {
		return dualsSqrtPoint;
	}
	public double[][] getDualsPoint() {
		return dualsPoint;
	}
	public double[] getDualsMaxCluster() {
		return dualsMaxCluster;
	}
	public Map<Integer, Map<Cluster, Double>> getDualsHierarchy() {
		return dualsHierarchy;
	}
	public double[] getDualsOneChangePerHierarchy() {
		return dualsOneChangePerHierarchy;
	}

	public double[] getDualsPointSingle() {
		return dualsPointSingle;
	}

	public void setDualsPointSingle(double[] dualsPointSingle) {
		this.dualsPointSingle = dualsPointSingle;
	}

	public double[] getDualsSqrtPointSingle() {
		return dualsSqrtPointSingle;
	}

	public void setDualsSqrtPointSingle(double[] dualsSqrtPointSingle) {
		this.dualsSqrtPointSingle = dualsSqrtPointSingle;
	}	


}