package cl.hierarchical.model.hp.qp;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import cl.data.GlobalParam;
import cl.util.Pair;

/**
 * RoofDualBound procedure based on the description in:
 * - A simple enumerative algorithm for unconstrained 0-1 quadratic programming by Pierre Hansen, Brigitte Jaumard, and Christophe Meyer
 */
public class Dinkelbach_QP_Hansen_RoofDualBound {
	
	private List<Integer> Gi;

	// Posiform
	private double[][] matrixW;
	private double[] vectorC;
	private double psiZero; // Current lowerbound

	// Graph
	private double[][] graph;
	private int sourceID;
	private int sinkID;

	public Dinkelbach_QP_Hansen_RoofDualBound(List<Integer> Gi,
			double[][] matrixW, double[] vectorC, double psiZero) {
		this.Gi = Gi;
		this.matrixW = matrixW;
		this.vectorC = vectorC;
		this.psiZero = psiZero;
	}

	/**
	 * Note that this implementation is different from Hansen 2000.
	 * In Hansen 2000: r is updated inside the if-statement
	 */
	public void solve()
	{
		int n = Gi.size();

		// Construct graph
		sourceID = 2*n; // Source is 2*n
		sinkID = 2*n+1;	// Sink is 2*n+1
		graph = new double[2*n+2][2*n+2]; // Upper triangular
		
		// Add edges
		// Note that Wi=0 (see page 4, Hansen 2000)
		for(int i = 0; i < 2*n; i++) {
			for(int j = i+1; j < 2*n; j++) {
				updateGraph(i, j, false);
			}
		}
		
		for(int i = 0; i < 2*n; i++) {
			updateGraph(i, false);
		}	

		DijkstraAlgorithm shortestPathSolver = new DijkstraAlgorithm(graph, 2*n+2);
		int r = 0;
		while(r<2*n) {
			if(vectorC[r]<=GlobalParam.EPSILON5) {
				r++;
				continue;
			}
			// Look for a SP
			List<Integer> path = shortestPathSolver.solve(r, sinkID);
			
			if(path.isEmpty()) {
				r++;
				continue;
			}
			path.add(0, sourceID);
			
			// Determine start and end node of the path
			int firstNode = path.get(1);
			int lastNode = path.get(path.size()-2);
			if(lastNode%2==0) {
				lastNode = lastNode + 1;
			}
			else {
				lastNode = lastNode -1;
			}
			List<Integer> firstAndLastNode = new ArrayList<>();
			firstAndLastNode.add(firstNode);
			firstAndLastNode.add(lastNode);

			// We have to check if c_r > epsilon and c_s > epsilon
			if(vectorC[firstNode]<=GlobalParam.EPSILON5
					|| vectorC[lastNode]<=GlobalParam.EPSILON5) {
				continue;
			}
	
			// Create new posiform, check for each element of the posiform what delta can maximally be
			double[][] newPosiformW = new double[2*n][2*n]; //New "matrixW"
			Set<Pair<Integer, Integer>> posiformID = new LinkedHashSet<>();

			int prev = -1; //null;
			for(int cur: path) {
				if(cur==sourceID || cur==sinkID) {
					continue;
				}
				if(prev==-1) {
					prev = cur;
					continue;
				}

				int prevId;
				if(prev%2==0) {
					prevId = prev+1;
				}
				else {
					prevId = prev-1;
				}

				// New literal
				int a = cur;
				int b = prevId;
				if(prevId<cur) {
					a = prevId;
					b = cur;
				}
				posiformID.add(new Pair<>(a, b));
				newPosiformW[a][b] += 1;
				prev = cur;
			}

			// Determine max feasible delta for W
			double maxDelta = Double.POSITIVE_INFINITY;
			for(Pair<Integer, Integer> pair: posiformID) {
				int a = pair.first;
				int b = pair.second;
				double delta = matrixW[a][b] / newPosiformW[a][b];
				if(delta < maxDelta) {
					maxDelta = delta;
				}
			}

			// Determine max feasible delta for C
			if(firstNode==lastNode) {
				double delta = vectorC[firstNode]/2d;
				if(delta<maxDelta) {
					maxDelta = delta;
				}
			}
			else {
				for(int node: firstAndLastNode) {
					if(vectorC[node]<maxDelta) {
						maxDelta = vectorC[node];
					}
				}
			}

			// Update posiform

			// Update W
			Set<Pair<Integer, Integer>> updateEdges = new LinkedHashSet<>(posiformID);
			for(Pair<Integer, Integer> pair: posiformID) {
				int a = pair.first;
				int b = pair.second;
				double coeff = maxDelta * newPosiformW[a][b];

				matrixW[a][b] -= coeff;

				if(a%2==0) {
					a++;
				}
				else {
					a--;
				}
				if(b%2==0) {
					b++;
				}
				else {
					b--;
				}

				updateEdges.add(new Pair<>(a, b));
				matrixW[a][b] += coeff;
			}

			// Update C
			if(firstNode==lastNode) {
				vectorC[firstNode] -= (maxDelta * 2);
			}
			else {
				for(int node: firstAndLastNode) {
					vectorC[node] -= maxDelta;
				}
			}

			// Update graph
			for(Pair<Integer, Integer> pair: updateEdges) {
				int a = pair.first;
				int b = pair.second;
				updateGraph(a, b, true);
			}

			for(int node: firstAndLastNode) {
				updateGraph(node, true);
			}

			// Update psiZero
			psiZero += maxDelta;
		}
	}

	private void updateGraph(int i, int j, boolean removePrevEdge) {
		// yk_bar -> yl
		int nodeStart;
		if(i%2==0) {
			nodeStart = i+1;
		}
		else {
			nodeStart = i-1;
		}
		int nodeEnd = j;
		double weight = matrixW[i][j]/2d;
		if(removePrevEdge) {
			graph[nodeStart][nodeEnd] = 0;
		}
		if(weight>GlobalParam.EPSILON5) {
			graph[nodeStart][nodeEnd] = weight;
		}

		// yl_bar -> yk
		if(j%2==0) {
			nodeStart = j+1;
		}
		else {
			nodeStart = j-1;	
		}
		nodeEnd = i;
		weight = matrixW[i][j]/2d;
		if(removePrevEdge) {
			graph[nodeStart][nodeEnd] = 0;
		}
		if(weight>GlobalParam.EPSILON5) {
			graph[nodeStart][nodeEnd] = weight;
		}
	}

	/**
	 * removePrevEdge=True if we want to remove edges already in the graph
	 * 
	 * @param i
	 * @param removePrevEdge
	 */
	private void updateGraph(int i, boolean removePrevEdge) {
		// Edges from source
		int node = i;
		double weight = vectorC[i]/2d; //Use this weight for all arcs
		if(removePrevEdge) {
			graph[sourceID][node] = 0;
		}
		if(weight>GlobalParam.EPSILON5) {
			graph[sourceID][node] = weight;
		}

		// Edges to sink
		if(i%2==0) {
			node = i+1;
		}
		else {
			node = i-1;
		}
		if(removePrevEdge) {
			graph[node][sinkID] = 0;
		}
		if(weight>GlobalParam.EPSILON5) {
			graph[node][sinkID] = weight;
		}
	}

	public double[][] getMatrixW() {
		return matrixW;
	}

	public double[] getVectorC() {
		return vectorC;
	}

	public double getPsiZero() {
		return psiZero;
	}	
}