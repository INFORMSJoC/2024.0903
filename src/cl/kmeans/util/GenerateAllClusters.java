package cl.kmeans.util;

import java.util.ArrayList;
import java.util.List;

import cl.data.Cluster;
import cl.data.Instance;
/**
 * Generates actual clusters
 * @author rickw
 *
 */
public class GenerateAllClusters {
	private Instance instance;
	private List<Cluster> clusters;

	public GenerateAllClusters(Instance instance) {
		this.instance = instance;
	}

	public List<Cluster> generate() {
		clusters = new ArrayList<>();
		for(int i = 0; i < instance.getNumPoints(); i++) {
			// Initialise cluster
			boolean[] currentCluster = new boolean[instance.getNumPoints()];
			currentCluster[i] = true;
			List<Integer> remainingVertices = new ArrayList<>();
			for(int j = i+1; j < instance.getNumPoints(); j++) {
				remainingVertices.add(j);
			}

			// Find all cliques that contain the vertex 
			recursiveSearch(currentCluster, remainingVertices);		
		}
		return clusters;
	}
	
	private void recursiveSearch(boolean[] currentCluster, List<Integer> remainingVertices) {
		clusters.add(new Cluster(currentCluster));
		
		if(remainingVertices.isEmpty()) {
			return;
		}
		
		List<Integer> nextRemainingVertices = new ArrayList<>(remainingVertices);

		// Add new vertex
		while(!nextRemainingVertices.isEmpty()) {
			int nextVertex = nextRemainingVertices.remove(0);
			boolean[] nextCluster = currentCluster.clone();
			nextCluster[nextVertex] = true;
			recursiveSearch(nextCluster, nextRemainingVertices);
		}
	}
}
