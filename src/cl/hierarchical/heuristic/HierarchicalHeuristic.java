package cl.hierarchical.heuristic;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cl.data.GlobalParam;
import cl.data.Instance;
import cl.data.Solution;
import cl.data.type.HierarchicalHeuristicType;
import cl.data.type.SolverType;
import cl.util.run.SpreadsheetGenerator;

public class HierarchicalHeuristic {
	private Instance instance;
	private HierarchicalHeuristicType hhType;
	
	private Solution bestSolution;
	private List<Solution> solutions;

	private static Logger log = LoggerFactory.getLogger(HierarchicalHeuristic.class);
	
	private long totalTime;

	public HierarchicalHeuristic(Instance instance, HierarchicalHeuristicType hhType) {
		this.instance = instance;
		this.hhType = hhType;
	}

	public void solve() {
		totalTime = System.currentTimeMillis();
		solutions = new ArrayList<>();
		bestSolution = new Solution(instance, Double.POSITIVE_INFINITY);
		log.info("Starting {} heuristic with {} restarts.", hhType, GlobalParam.HEURISTIC_RESTARTS);
		
		for(int i = 0; i < GlobalParam.HEURISTIC_RESTARTS; i++) {
			if(hhType==HierarchicalHeuristicType.hmc) {
				for(int k = 2; k <= instance.getNumClusters(); k++) {
					HMC hmc = new HMC(instance, k, instance.getNumClusters(), SolverType.Heuristic_HMC);
					hmc.solve();
					Solution sol = hmc.getSolution();
					solutions.add(sol);
					if(sol.getObjective() < bestSolution.getObjective()) {
						log.info("Using HMC with k {}, we found a better hierarchical solution with objective {}.", k, sol.getObjective());
						bestSolution = sol;
						if(GlobalParam.WRITE_BETTER_INTEGER_SOLUTION) {
							writeBestSolution();
						}
					}
				}
			}
			else if(hhType==HierarchicalHeuristicType.agglomerative) {
				HMC hmc = new HMC(instance, instance.getNumClusters(), instance.getNumClusters(), SolverType.Heuristic_Agg);
				hmc.solve();
				Solution sol = hmc.getSolution();
				solutions.add(sol);
				if(sol.getObjective() < bestSolution.getObjective()) {
					log.info("Using agglomerative heuristic, we found a better hierarchical solution with objective {}.", sol.getObjective());
					bestSolution = sol;
					if(GlobalParam.WRITE_BETTER_INTEGER_SOLUTION) {
						writeBestSolution();
					}
				}
			}
			else if(hhType==HierarchicalHeuristicType.divisive) {
				HMC hmc = new HMC(instance, 2, instance.getNumClusters(), SolverType.Heuristic_Div);
				hmc.solve();
				Solution sol = hmc.getSolution();
				solutions.add(sol);
				if(sol.getObjective() < bestSolution.getObjective()) {
					log.info("Using divisive heuristic, we found a better hierarchical solution with objective {}.", sol.getObjective());
					bestSolution = sol;
					if(GlobalParam.WRITE_BETTER_INTEGER_SOLUTION) {
						writeBestSolution();
					}
				}
			}
			else {
				throw new IllegalArgumentException("This method does not have an implementation yet");
			}
		}
		log.info("Best objective using {} heuristic is {}", hhType, bestSolution.getObjective());
	}
	
	private void writeBestSolution() {
		String date_file = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date());
		File outputFile = new File(GlobalParam.WRITE_BETTER_INTEGER_SOLUTION_FILE_NAME+"/java_"+date_file+".xlsx");
		List<Solution> tempList = new ArrayList<>();
		bestSolution = getBestSolution();
		bestSolution.setTotalTime((System.currentTimeMillis() - totalTime));
		tempList.add(bestSolution);
		Map<File, List<Solution>> tempMap = new LinkedHashMap<>();
		tempMap.put(outputFile, tempList);
		try {
			SpreadsheetGenerator.writeSpreadsheet(outputFile, tempMap, false);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public Solution getBestSolution() {
		if(hhType==HierarchicalHeuristicType.hmc) {
			bestSolution.setSolverType(SolverType.Heuristic_HMC);
		}
		else if(hhType==HierarchicalHeuristicType.agglomerative) {
			bestSolution.setSolverType(SolverType.Heuristic_Agg);
		}
		else if(hhType==HierarchicalHeuristicType.divisive) {
			bestSolution.setSolverType(SolverType.Heuristic_Div);
		}
		return bestSolution;
	}

	public List<Solution> getSolutions() {
		return solutions;
	}
}
