package cl.hierarchical.model.hp.qp;

import java.util.ArrayList;
import java.util.List;

import cl.hierarchical.model.hp.qp.FibonacciHeap.EntryHeap;
/**
 * A fast implementation of Dijkstra's algorithm, based on the description in:
 * https://www.boost.org/doc/libs/1_39_0/libs/graph/doc/dijkstra_shortest_paths.html
 */
public class DijkstraAlgorithm {

	private double[][] network;
	private int n;

	private int[] shortestPathCovered;
	private int[] prev;
	private EntryHeap[] entries;
	
	private FibonacciHeap<Integer> queue;

	public DijkstraAlgorithm(double[][] network, int n) {
		this.network = network;
		this.n = n;

		initNetwork();
	}

	private void initNetwork() {
		// Initialise
		shortestPathCovered = new int[n];
		prev = new int[n];
		entries = new EntryHeap[n];
	}

	private void reset() {
		queue = new FibonacciHeap<>();
		for (int i = 0; i < n; i++) { 
			shortestPathCovered[i] = 0;
			prev[i] = -1;
			entries[i] = new EntryHeap<Integer>(i, Double.POSITIVE_INFINITY);
		} 
	}

	/**
	 * Find shortest path from source to sink
	 * 
	 * If entirePath: return entire path
	 * Else only return the first and last node
	 * 
	 * @param source
	 * @param sink
	 */
	public List<Integer> solve(int source, int sink) {
		// Initialise
		reset();
		entries[source].setPriority(0);
		queue.enqueue(source, 0);
			
		// Find the shortestPath for all other vertices
		while(!queue.isEmpty()) {
			// Closest vertex
			int u = minDistance(); 
			shortestPathCovered[u] = 2;
			if(u==sink) {
				// We found the sink, so we can stop
				break;
			}

			// process adjacent nodes of the current vertex
			if(Double.isFinite(entries[u].getPriority())) {
				for(int v = 0; v < n; v++) {
					EntryHeap entry = entries[v];
					// if vertex v not covered then update it  
					if (shortestPathCovered[v]!=2 && network[u][v] != 0 &&  entries[u].getPriority() + network[u][v] < entries[v].getPriority()) {
						double shortestPathsV = entries[u].getPriority() + network[u][v]; 
						prev[v] = u;
						if(shortestPathCovered[v]==0) {
							shortestPathCovered[v] = 1;
							EntryHeap<Integer> entryHeap = entries[v];
							entryHeap.setPriority(shortestPathsV);
							queue.enqueue(entryHeap);
						}
						else if(shortestPathCovered[v]==1) {
							queue.decreaseKey(entries[v], shortestPathsV);
						}
					}
				}
			}
		} 

		// get path
		List<Integer> path = new ArrayList<>();
		if(prev[sink]!=-1) {
			int u = sink;
			while(true) {
				if(u==source || u==-1) {
					break;
				}
				path.add(u);
				u = prev[u];
			}
			path.add(source);

			// Reverse direction
			List<Integer> reversePath = new ArrayList<>();
				for(int i = 0; i < path.size(); i++) {
					reversePath.add(path.get(path.size()-i-1));
				}
			path = reversePath;
		}
		return path;
	}

	private int minDistance()   { 
		// Initialize min value 
		return queue.dequeueMin().getValue();
	} 
}
