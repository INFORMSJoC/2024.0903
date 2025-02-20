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
import cl.util.run.SpreadsheetGenerator;
import ilog.concert.IloException;

public class BatchRunHeuristic {

	private static Logger log = LoggerFactory.getLogger(BatchRunHeuristic.class);

	public static void main(String[] args) throws IloException, IOException {
		BasicConfigurator.configure();
		File inputDir = new File("data/batchinstances");
		String date_file = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date());
		File outputDir = new File("data/batchsolutions/"+date_file);
		outputDir.mkdirs();

		int numClusters = 5;

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
				
				long startHeuristicTime = System.currentTimeMillis();
				HierarchicalHeuristic hh = new HierarchicalHeuristic(instance, HierarchicalHeuristicType.agglomerative);
				hh.solve();	
				long heuristicTime = System.currentTimeMillis() - startHeuristicTime;
				
				Solution solution = hh.getBestSolution();
				List<Double> startObjectives = new ArrayList<>();
				for(Solution sol: hh.getSolutions()) {
					startObjectives.add(sol.getObjective());
				}
				solution.setStartObjectives(startObjectives);
				solution.setTotalTime(heuristicTime);
				solutionList.add(solution);
				
				File outputFile = new File(outputDir+"/java_"+date_file+"_"+count+".xlsx");
				List<Solution> tempList = new ArrayList<>();
				tempList.add(solution);
				Map<File, List<Solution>> tempMap = new LinkedHashMap<>();
				tempMap.put(outputFile, tempList);
				SpreadsheetGenerator.writeSpreadsheet(outputFile, tempMap, true);
				count++;
			}
			solutions.put(tempDir, solutionList);
		}

		File outputFile = new File(outputDir+"/java_"+date_file+".xlsx");
		File avgOutputFile = new File(outputDir+"/java_avg_"+date_file+".xlsx");
		SpreadsheetGenerator.writeSpreadsheet(outputFile, solutions, true);
		SpreadsheetGenerator.writeSpreadsheetAverage(avgOutputFile, solutions, true);
	}
}