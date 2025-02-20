package cl.hierarchical.model.compact;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import cl.data.Cluster;
import cl.data.Instance;
import cl.data.Solution;
import cl.data.type.SolverType;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;
import cl.data.GlobalParam;
/**
 * Compact formulation 1
 * @author rickw
 *
 */
public class MIP_Compact_MSS_Hierarchical {

	private Instance instance;
	private IloCplex model;
	private IloObjective objective;

	private int tStart, tEnd;

	private final double epsilon = 0.00000001;

	private IloNumVar[][][] variablesZ;
	private IloNumVar[][][] variablesZeta;
	
	private Solution bestSolution = new Solution(instance, Double.POSITIVE_INFINITY);

	private boolean solveLPRelaxation = false;
	private boolean useValidInequalities = true;
	
	public MIP_Compact_MSS_Hierarchical(Instance instance) throws IloException {
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
		initVariablesZeta();
		initConstraints2();
		initConstraints3();
		initConstraints4();
		initConstraints9();
		initConstraints10();
		if(useValidInequalities) {
			initConstraints6();
			initConstraints7();
			initConstraints11();
			initConstraints12();
		}
		initConstraintsHierarchical();
		initObjective();
	}
	
	/**
	 * https://www.ibm.com/docs/en/icos/22.1.0?topic=mip-starting-from-solution-starts
	 * @param startSolution
	 * @throws IloException 
	 */
	public void initStartSolution(Solution startSolution) throws IloException {
		IloNumVar[] startVar = new IloNumVar[instance.getNumPoints() * instance.getNumPoints() * instance.getNumClusters()];
		double[] startVal = new double[instance.getNumPoints() * instance.getNumPoints() * instance.getNumClusters()];
		int idx = 0;
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					startVar[idx] = variablesZeta[i][j][t];
					
					// Check if i and j are in the same cluster
					for(Cluster cluster: startSolution.getHierarchicalClustersList().get(t)) {
						if(cluster.getCluster()[i] && cluster.getCluster()[j]) {
							startVal[idx] = 1;
							break;
						}
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
	 * z_ij = z_ji
	 * @throws IloException
	 */
	private void initConstraints6() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = i+1; j < instance.getNumPoints(); j++) {
					model.addEq(variablesZ[i][j][t], variablesZ[j][i][t]);				
				}
			}
		}
	}

	/**
	 * z_ij + z_il - z_jl <= z_ii
	 * @throws IloException 
	 */
	private void initConstraints7() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					for(int k = 0; k < instance.getNumPoints(); k++) {
						if(i==j) {
							continue;
						}
						IloNumExpr expr = model.constant(0);
						expr = model.sum(expr, variablesZ[i][j][t]);
						expr = model.sum(expr, variablesZ[i][k][t]);
						expr = model.diff(expr, variablesZ[j][k][t]);
						model.addLe(expr, variablesZ[i][i][t]);
					}
				}
			}
		}
	}

	/**
	 * z_ij <= zeta_ij
	 * @throws IloException 
	 */
	private void initConstraints9() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					model.addLe(variablesZ[i][j][t], variablesZeta[i][j][t]);
				}
			}
		}
	}

	private void initConstraints10() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					if(i==j) {
						continue;
					}

					IloNumExpr lhs = model.diff(variablesZ[i][i][t], variablesZ[i][j][t]);
					IloNumExpr rhs = model.diff(1, variablesZeta[i][j][t]);
					model.addLe(lhs, rhs);
				}
			}
		}
	}

	/**
	 * (N-K+1)z_ij >= zeta_ij
	 * @throws IloException 
	 */
	private void initConstraints11() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					IloNumExpr expr = model.prod((instance.getNumPoints()-(t+1)+1), variablesZ[i][j][t]);
					model.addGe(expr, variablesZeta[i][j][t]);
				}
			}
		}
	}

	/**
	 * (N-K+1)(z_ii-z_ij) >= 1-zeta_ij
	 * @throws IloException 
	 */
	private void initConstraints12() throws IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					if(i==j) { 
						continue;
					}
					IloNumExpr temp = model.diff(variablesZ[i][i][t], variablesZ[i][j][t]);
					IloNumExpr lhs = model.prod((instance.getNumPoints()-(t+1)+1), temp);
					IloNumExpr rhs = model.diff(1, variablesZeta[i][j][t]);
					model.addGe(lhs, rhs);
				}
			}
		}
	}

	/**
	 * zeta_ijt >= zeta_ij(t+1)
	 * @throws IloException 
	 */
	private void initConstraintsHierarchical() throws IloException {
		for(int i = 0; i < instance.getNumPoints(); i++) {
			for(int j = 0; j < instance.getNumPoints(); j++) {
				for(int t = tStart; t < tEnd-1; t++) {
					model.addGe(variablesZeta[i][j][t], variablesZeta[i][j][t+1]);
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

	private void initVariablesZeta() throws IloException {
		variablesZeta = new IloIntVar[instance.getNumPoints()][instance.getNumPoints()][instance.getNumClusters()];
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					variablesZeta[i][j][t] = model.boolVar();
					if(solveLPRelaxation) {
						variablesZeta[i][j][t] = model.numVar(0, 1);
					}
				}
			}
		}
	}

	public void solve() throws IloException {
		long computationTime = System.currentTimeMillis();
		model.solve();
		computationTime = System.currentTimeMillis() - computationTime;
		
		if(solveLPRelaxation) {
			bestSolution = new Solution(instance, model.getObjValue());
		}
		else if(model.isPrimalFeasible()) {
			bestSolution = getCurrentSolution();
			bestSolution.setFeasible(true);
		}
		else {
			bestSolution = new Solution(instance, new LinkedHashMap<>(), Double.POSITIVE_INFINITY);
		}
		bestSolution.setTotalTime(computationTime);
		bestSolution.setTimeLimitReached(computationTime>GlobalParam.GLOBAL_TIME_LIMIT);
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

	public void cleanUp() throws IloException {
		model.clearModel();
		model.end();
	}

	public Solution getBestSolution() {
		bestSolution.setSolverType(SolverType.MIP1);
		return bestSolution;
	}

	public double getObjectiveValue() throws IloException {
		return model.getObjValue();
	}

	public void printSolution() throws IloException {
		System.out.println("is feasible? "+model.isPrimalFeasible());
		System.out.println("objective: "+ model.getObjValue());
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
					if (Math.round(model.getValue(variablesZeta[i][j][t]))==1) {
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

	public void checkVariablesZ() throws UnknownObjectException, IloException {
		for(int t = tStart; t < tEnd; t++) {
			for(int i = 0; i < instance.getNumPoints(); i++) {
				System.out.print(i + ": ");
				for(int j = 0; j < instance.getNumPoints(); j++) {
					if(model.getValue(variablesZ[i][j][t]) > epsilon) {
						System.out.print(j + " ");
					}
				}
				System.out.println();
			}
		}
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

	public void printVariablesZeta() throws UnknownObjectException, IloException {
		System.out.println("Zeta");
		for(int t = tStart; t < tEnd; t++) {
			System.out.println("T: "+(t+1));
			for(int i = 0; i < instance.getNumPoints(); i++) {
				for(int j = 0; j < instance.getNumPoints(); j++) {
					System.out.print(Math.round(model.getValue(variablesZeta[i][j][t])) + " ");
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