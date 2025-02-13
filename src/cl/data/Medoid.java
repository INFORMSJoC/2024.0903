package cl.data;

import java.util.List;

public class Medoid {

	private List<Integer> cluster;
	private double[] medoidCoord;
	private int medoidId;

	public Medoid(List<Integer> cluster, double[] medoidCoord, int medoidId) {
		this.cluster = cluster;
		this.medoidCoord = medoidCoord;
		this.medoidId = medoidId; // MedoidID corresponds with the same medoidCoord in the Instance
	}

	public List<Integer> getCluster() {
		return cluster;
	}

	public double[] getMedoidCoord() {
		return medoidCoord;
	}

	public int getMedoidId() {
		return medoidId;
	}

	public void setMedoidCoord(double[] medoidCoord) {
		this.medoidCoord = medoidCoord;
	}

	public void setMedoidId(int medoidId) {
		this.medoidId = medoidId;
	}
}