package cl.hierarchical.model.compact;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import cl.data.Cluster;
import cl.data.GlobalParam;
import cl.data.Instance;
import cl.data.Solution;
import cl.data.type.SolverType;
import ilog.concert.IloException;
import ilog.concert.IloIntExpr;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;
/**
 * Compact formulation 2
 * @author rickw
 *
 */
public class MIP_Compact_AMSS_Hierarchical {

	private Instance instance;
	private IloCplex model;
	private IloObjective objective;

	private int tStart, tEnd;

	private IloNumVar[][][] variablesZ;
	private IloIntVar[][][] variablesGamma;
	
	private Solution bestSolution = new Solution(instance, Double.POSITIVE_INFINITY);

	public MIP_Compact_AMSS_Hierarchical(Instance instance) throws IloException {
		this.instance = instance;
		this.model = new IloCplex();
		this.model.setOut(null);	
		this.model.setParam(IloCplex.Param.TimeLimit, GlobalParam.GLOBAL_TIME_LIMIT / 1000d); // In seconds
		initModel();	
	}

	private void initModel() throws IloException {
		tStart = 0;
		tEnd = instance.getNumClusters();

		initVariablesZ();
		initVariablesGamma();
		initConstraints2();
		initConstraints3();
		initConstraints4();
		initConstraints14();
		initConstraints15();
		initConstraints16();
		initConstraintsHierarchical();
		initObjective();
	}
	
	public void initStartSolution(Solution startSolution) throws IloException {
		IloNumVar[] startVar = new IloNumVar[instance.getNumPoints() * instance.getNumClusters() * instance.getNumClusters()];
		double[] startVal = new double[instance.getNumPoints() * instance.getNumClusters() * instance.getNumClusters()];
		int idx = 0;
		for(int t = tStart; t < tEnd; t++) {
			List<Cluster> clusterList = startSolution.getHierarchicalClustersList().get(t);
			for(int i = 0; i < instance.getNumClusters(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					startVar[idx] = variablesGamma[j][i][t];
					
					// Check if point j is in cluster i
					if(i<clusterList.size() && clusterList.get(i).getCluster()[j]) {
						startVal[idx] = 1;
					}
					
					idx++;
				}
			}
		}
		model.addMIPStart(startVar, startVal);
	}

	/**
	 * sum_j z_ij = 1
	 * @throws IloException 
	 */
	private void initConstraints2() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				IloNumExpr expr = model.constant(0);
				for(int j = 0; j < instance.getNumPoints(); j++) {
					expr = model.sum(expr, variablesZ[i][j][t]);
				}
				model.addEq(expr, 1);
			}
		}
	}

	/**
	 * z_ij <= z_ii
	 * @throws IloException 
	 */
	private void initConstraints3() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					if(i==j) {
						continue;
					}
					model.addLe(variablesZ[i][j][t], variablesZ[i][i][t]); 
				}
			}
		}
	}

	private void initConstraints4() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			IloNumExpr expr = model.constant(0);
			for(int i = 0; i < instance.getNumPoints(); i++) {
				expr = model.sum(expr, variablesZ[i][i][t]);
			}
			model.addEq(expr, t+1);
		}
	}


	/**
	 * sum_k gamma_ik = 1
	 * @throws IloException
	 */
	private void initConstraints14() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				IloIntExpr expr = model.constant(0);
				for(int j = 0; j < t+1; j++) {
					expr = model.sum(expr, variablesGamma[i][j][t]);
				}
				model.addEq(expr, 1);
			}
		}
	}

	/**
	 * sum_i gamma_ik >= 1
	 * @throws IloException
	 */
	private void initConstraints15() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < (t+1); i++) {
				IloIntExpr expr = model.constant(0);
				for(int j = 0; j < instance.getNumPoints(); j++) {
					expr = model.sum(expr, variablesGamma[j][i][t]);
				}
				model.addGe(expr, 1);
			}
		}
	}

	/**
	 * z_ij <= 1 + gamma_ik - gamma_jk
	 * @throws IloException 
	 */
	private void initConstraints16() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					if(i==j) {
						continue;
					}
					for(int k = 0; k < (t+1); k++) {
						IloIntExpr expr = model.sum(1, variablesGamma[i][k][t]);
						expr = model.diff(expr, variablesGamma[j][k][t]);
						model.addLe(variablesZ[i][j][t], expr);
					}
				}
			}
		}
	}

	/**
	 * z_ijt >= 1/N (gamma_ik(t+1) + gamma_jk(t+1) - 1)
	 * @throws IloException 
	 */
	private void initConstraintsHierarchical() throws IloException {
		for(int i = 0; i < instance.getNumPoints(); i++) {
			for(int j = 0; j < instance.getNumPoints(); j++) {
				for(int t = tStart; t < tEnd-1; t++) {
					for(int k = 0; k < (t+1)+1; k++) { // (t+1)+1: because of the number of clusters in the previous period
						IloNumExpr temp = model.sum(variablesGamma[i][k][t+1], variablesGamma[j][k][t+1]);
						temp = model.diff(temp, 1);
						temp = model.prod(1d/instance.getNumPoints(), temp);
						model.addGe(variablesZ[i][j][t], temp);
					}
				}
			}
		}
	}

	private void initObjective() throws IloException {
		IloNumExpr expr = model.constant(0);
		for(int t = tStart; t < tEnd; t++) {
			IloNumExpr tempT = model.constant(0);
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					if(i<=j) {
						continue;
					}
					IloNumExpr temp = model.prod(instance.getSquaredDistanceMatrix()[i][j], variablesZ[i][j][t]); 
					tempT = model.sum(tempT, temp);
				}
			}
			expr = model.sum(expr, tempT);
		}
		objective = model.addMinimize(expr);
	}

	private void initVariablesZ() throws IloException {
		variablesZ = new IloNumVar[instance.getNumPoints()][instance.getNumPoints()][instance.getNumClusters()];
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					variablesZ[i][j][t] = model.numVar(0, 1);
				}
			}
		}
	}

	private void initVariablesGamma() throws IloException {
		variablesGamma = new IloIntVar[instance.getNumPoints()][instance.getNumClusters()][instance.getNumClusters()];
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumClusters(); j++) {
					variablesGamma[i][j][t] = model.boolVar();
					
					model.addLe(variablesGamma[i][j][t], 1); // TODO: maybe redundant
				}
			}
		}
	}

	public void solve() throws IloException {
		long computationTime = System.currentTimeMillis();
		model.solve();
		computationTime = System.currentTimeMillis() - computationTime;
		
		if(model.isPrimalFeasible()) {
			bestSolution = getCurrentSolution();
			bestSolution.setFeasible(true);
		}
		else {
			bestSolution = new Solution(instance, new LinkedHashMap<>(), Double.POSITIVE_INFINITY);
		}
		bestSolution.setTotalTime(computationTime);
		bestSolution.setTimeLimitReached(computationTime>GlobalParam.GLOBAL_TIME_LIMIT);
	}

	public void cleanUp() throws IloException {
		model.clearModel();
		model.end();
	}
	
	private Solution getCurrentSolution() throws UnknownObjectException, IloException {
		Map<Integer, List<Cluster>> hierarchicalClusters = new LinkedHashMap<>();
		Map<Integer, Map<Integer, Set<Integer>>> dendogram = getDendogram();
		for(Entry<Integer, Map<Integer, Set<Integer>>> entry: dendogram.entrySet()) {
			List<Cluster> clusters = new ArrayList<>();
			for(Set<Integer> entry2: entry.getValue().values()) {
				boolean[] cluster = new boolean[instance.getNumPoints()];
				for(Integer i: entry2) {
					cluster[i] = true;
				}
				clusters.add(new Cluster(cluster));
			}
			hierarchicalClusters.put(entry.getKey()-1, clusters); //TOOD: level starts at 1
		}
		return new Solution(instance, hierarchicalClusters, model.getObjValue());
	}
	
	public Solution getBestSolution() {
		bestSolution.setSolverType(SolverType.MIP2);
		return bestSolution;
	}
	
	public Map<Integer, Map<Integer, Set<Integer>>> getDendogram() throws UnknownObjectException, IloException {
		Map<Integer, Map<Integer, Set<Integer>>> dendrogram = new LinkedHashMap<>();
		for(int t = tStart; t < tEnd; t++) {
			Set<Integer> notCovered =  Arrays.stream(IntStream.range(0, instance.getNumPoints()).toArray()).boxed().collect(Collectors.toSet());

			Map<Integer, Set<Integer>> mapT = new LinkedHashMap<>();

			boolean cont = true;
			while(cont) {
				if(notCovered.isEmpty()) {
					break;
				}

				Set<Integer> covered = new LinkedHashSet<>();
				int i = notCovered.iterator().next();
				for(Integer j: notCovered) {
					if ((model.getValue(variablesZ[i][j][t]))>=(1d/(instance.getNumPoints()+0.1))) { //Because of rounding errors
						covered.add(j);
					}
				}

				mapT.put(i, covered);
				notCovered.removeAll(covered);
			}
			dendrogram.put((t+1), mapT);
		}
		return dendrogram;
	}

	public void printSolution() throws IloException {
		System.out.println("is feasible? "+model.isPrimalFeasible());
		System.out.println("objective: "+ model.getObjValue());
	}

	public double getObjectiveValue() throws IloException {
		return model.getObjValue();
	}
	
	public void printVariablesZ() throws UnknownObjectException, IloException {
		System.out.println("Z");
		for(int t = tStart; t < tEnd; t++) {
			System.out.println("T: "+(t+1));
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					System.out.print(model.getValue(variablesZ[i][j][t]) + " ");
				}
				System.out.println();
			}
		}
	}
	
	/**
	 * Cluster ID is equal to the (lowest) member ID
	 * @throws UnknownObjectException
	 * @throws IloException
	 */
	public void printDendrogram() throws UnknownObjectException, IloException {
		Map<Integer, Map<Integer, Set<Integer>>> dendrogram = getDendogram();
		for(Entry<Integer, Map<Integer, Set<Integer>>> mapT: dendrogram.entrySet()) {
			System.out.println("T: "+mapT.getKey());

			for(Entry<Integer, Set<Integer>> cluster: mapT.getValue().entrySet()) {
				System.out.println("Cluster ID: "+ cluster.getKey() + ", " + Arrays.toString(cluster.getValue().toArray()));
			}
		}
	}
}