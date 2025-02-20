package cl.data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonProperty;

public class Cluster {
	private boolean[] cluster;
	private int[] pointId;
	private int numPoints;
	private double inactiveIter = 0; // For column management
		
	public Cluster(boolean[] cluster) {
		this.cluster = cluster;
		this.numPoints = cluster.length;
		this.pointId = computePointId();
	}
	
	public Cluster(int[] pointId, int numPoints) {
		this.pointId = pointId;
		this.numPoints = numPoints;
		this.cluster = computeCluster();
	}
	
	@JsonCreator
	public Cluster(@JsonProperty("poindId") int[] pointId, @JsonProperty("numPoints") int numPoints, @JsonProperty("cluster") boolean[] cluster) {
		this.pointId = pointId;
		this.numPoints = numPoints;
		this.cluster = cluster;
	}
	
	@JsonIgnore
	public double getInactiveIter() {
		return inactiveIter;
	}

	public void setInactiveIter(double inactiveIter) {
		this.inactiveIter = inactiveIter;
	}
	
	public Cluster copy() {
		return new Cluster(pointId, numPoints, cluster);
	}

	public Cluster deepCopy() {
		return new Cluster(pointId.clone(), numPoints, cluster.clone());
	}

	/**
	 * E.g. return [1, 3, 5, 9, ...]
	 * @return
	 */
	public int[] getPointId() {
		return pointId;
	}
	
	public boolean[] getCluster() {
		return cluster;
	}

	/**
	 * Get number of points of instance
	 * @return
	 */
	public int getNumPoints() {
		return numPoints;
	}

	private int[] computePointId() {
		List<Integer> temp = new ArrayList<>();
		for(int i = 0; i < cluster.length; i++) {
			if(cluster[i]) {
				temp.add(i);
			}
		}
		return temp.stream().mapToInt(i -> i).toArray();
	}
	
	private boolean[] computeCluster() {
		boolean[] temp = new boolean[numPoints];
		for(int i: pointId) {
			temp[i] = true;
		}
		return temp;
	}
	
	public static Cluster merge(Cluster c1, Cluster c2) {
		boolean[] bool1 = c1.getCluster();
		boolean[] bool2 = c2.getCluster();
		
		boolean[] boolNew = new boolean[c1.getNumPoints()];
		for(int i = 0; i < c1.getNumPoints(); i++) {
			if(bool1[i] || bool2[i]) {
				boolNew[i] = true;
			}
		}
		return new Cluster(boolNew);
	}
	
	/**
	 * Return true if this cluster is a subset of a larger cluster
	 * @param superSet
	 * @return
	 */
	public boolean isSubsetOf(Cluster otherCluster) {
		for(int i: pointId) {
			if(!otherCluster.getCluster()[i]) {
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Return true if cluster 1 and cluster 2 do not have any objects in common
	 * @return
	 */
	public static boolean intersectionIsEmpty(Cluster c1, Cluster c2) {
		for(int i: c1.getPointId()) {
			for(int j: c2.getPointId()) {
				if(i>j) {
					continue;
				}
				if(i==j) {
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * We already know that originalCluster contains object i
	 * 
	 * Return true, if
	 * - otherCluster contain object i
	 * - originalCluster is not a subset of otherCluster
	 * - originalCluster is not equal to otherCluster
	 * 
	 * In other words, there must be at least one object j which is in originalCluster, but not in otherCluster
	 * 
	 * @param c1
	 * @param c2
	 * @param object
	 * @return
	 */
	public static boolean checkValidInequalityType1(Cluster originalCluster, Cluster otherCluster, int object) {
		boolean check = false;
		for(int i: originalCluster.getPointId()) {
			if(!otherCluster.getCluster()[i]) {
				check = true;
			}
		}
		
		if(!check) {
			return false;
		}
		
		for(int i: otherCluster.getPointId()) {
			if(i==object) {
				return true;
			}
			if(i>object) {
				return false;
			}
		}
		return false;
	}
	
	/**
	 * Return true if this cluster contains element i
	 * @param i
	 * @return
	 */
	public boolean contains(int i) {
		for(int j: pointId) {
			if(j==i) {
				return true;
			}
			if(j>i) {
				return false;
			}
		}
		return false;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(cluster);
		result = prime * result + numPoints;
		result = prime * result + Arrays.hashCode(pointId);
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
		Cluster other = (Cluster) obj;
		if (!Arrays.equals(cluster, other.cluster))
			return false;
		if (numPoints != other.numPoints)
			return false;
		if (!Arrays.equals(pointId, other.pointId))
			return false;
		return true;
	}

	@Override
	public String toString() {
		return Arrays.toString(pointId);
	}
}