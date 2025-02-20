package cl.hierarchical.run;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.BasicConfigurator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cl.data.Instance;
import cl.data.Solution;
import cl.data.type.HierarchicalHeuristicType;
import cl.hierarchical.heuristic.HierarchicalHeuristic;
import cl.hierarchical.model.compact.MIP_Compact_AMSS_Hierarchical;
import cl.hierarchical.model.compact.MIP_Compact_MSS_Hierarchical;
import cl.util.run.SpreadsheetGenerator;
import ilog.concert.IloException;

public class BatchRunMIP {

	private static Logger log = LoggerFactory.getLogger(BatchRunMIP.class);

	public static void main(String[] args) throws IloException, IOException {
		BasicConfigurator.configure();
		File inputDir = new File("data/batchinstances");
		String date_file = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date());
		boolean MIP1 = true;
		if(MIP1) {
			date_file += "_MIP1";
		}
		else {
			date_file += "_MIP2";
		}
		File outputDir = new File("data/batchsolutions/"+date_file);
		outputDir.mkdirs();

		int numClusters = 5;
		
		// Add the best start solution
		HierarchicalHeuristicType hierarchicalHeuristicType = HierarchicalHeuristicType.hmc;

		int count = 0;
		Map<File, List<Solution>> solutions = new LinkedHashMap<>();
		for (File tempDir : inputDir.listFiles()) {
			if (!tempDir.isDirectory()) {
				throw new IllegalArgumentException("Place instances with same setting in a folder");
			}
			List<Solution> solutionList = new ArrayList<>();
			for(File f : tempDir.listFiles()) {
				log.info("Reading file {}", f);
				Instance instance = Instance.readFromTxt(f, numClusters);
				
				HierarchicalHeuristic hh = null;
				if(hierarchicalHeuristicType!=HierarchicalHeuristicType.none) {
					hh = new HierarchicalHeuristic(instance, hierarchicalHeuristicType);
					hh.solve();	
				}
				
				Solution solution;
				if(MIP1) {
					MIP_Compact_MSS_Hierarchical solver = new MIP_Compact_MSS_Hierarchical(instance); // MIP1
					if(hierarchicalHeuristicType!=HierarchicalHeuristicType.none) {
						solver.initStartSolution(hh.getBestSolution());
					}
					solver.solve();
					solution = solver.getBestSolution();
					solutionList.add(solution);
					solver.cleanUp();
				}
				else {
					MIP_Compact_AMSS_Hierarchical solver = new MIP_Compact_AMSS_Hierarchical(instance); // MIP2	
					if(hierarchicalHeuristicType!=HierarchicalHeuristicType.none) {
						solver.initStartSolution(hh.getBestSolution());
					}
					solver.solve();
					solution = solver.getBestSolution();
					solutionList.add(solution);
					solver.cleanUp();
				}
				
				File outputFile = new File(outputDir+"/java_"+date_file+"_"+count+".xlsx");
				List<Solution> tempList = new ArrayList<>();
				tempList.add(solution);
				Map<File, List<Solution>> tempMap = new LinkedHashMap<>();
				tempMap.put(outputFile, tempList);
				SpreadsheetGenerator.writeSpreadsheet(outputFile, tempMap, false);
				count++;
			}
			solutions.put(tempDir, solutionList);
		}

		File outputFile = new File(outputDir+"/java_"+date_file+".xlsx");
		File avgOutputFile = new File(outputDir+"/java_avg_"+date_file+".xlsx");
		SpreadsheetGenerator.writeSpreadsheet(outputFile, solutions, false);
		SpreadsheetGenerator.writeSpreadsheetAverage(avgOutputFile, solutions, false);
	}
}