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

import com.fasterxml.jackson.databind.ObjectMapper;

import cl.data.Instance;
import cl.data.Solution;
import cl.data.type.HierarchicalHeuristicType;
import cl.data.type.PricingProblemStrategy;
import cl.data.type.SolverType;
import cl.hierarchical.heuristic.HierarchicalHeuristic;
import cl.hierarchical.model.pdcg.PrimalDualColgen;
import cl.util.run.SpreadsheetGenerator;
import ilog.concert.IloException;

public class BatchRunPDCGM {

	private static Logger log = LoggerFactory.getLogger(BatchRunPDCGM.class);

	public static void main(String[] args) throws IloException, IOException {
		BasicConfigurator.configure();
		
		// Main parameters
		int numClusters = 5;
		HierarchicalHeuristicType hierarchicalHeuristicType = HierarchicalHeuristicType.hmc;
		SolverType solverType = SolverType.Colgen_BranchAndPrice;
		PricingProblemStrategy pricingProblemStrategy = PricingProblemStrategy.None;
		
		File inputDir = new File("data/batchinstances");
		String date_file = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date());
		date_file += "_"+solverType;
		File outputDir = new File("data/batchsolutions/"+date_file);
		outputDir.mkdirs();
		ObjectMapper om = new ObjectMapper();
		
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
				HierarchicalHeuristic hh = null;
				if(hierarchicalHeuristicType!=HierarchicalHeuristicType.none) {
					hh = new HierarchicalHeuristic(instance, hierarchicalHeuristicType);
					hh.solve();	
				}
				long heuristicTime = System.currentTimeMillis() - startHeuristicTime;
				
				PrimalDualColgen model = null; 
				model = new PrimalDualColgen(instance, solverType, pricingProblemStrategy);
				model.initHierarchicalStartSolutions(hh.getSolutions());
				model.solve();
				
				Solution bestSolution = model.getBestSolution(true);
				bestSolution.setHeuristicTime(heuristicTime);
				bestSolution.setHierarchicalHeuristicType(hierarchicalHeuristicType);
				
				solutionList.add(bestSolution);
				
				File outputFile = new File(outputDir+"/java_"+date_file+"_"+count+".xlsx");
				File outputFileJSON = new File(outputDir+"/java_"+date_file+"_"+count+".json");
				List<Solution> tempList = new ArrayList<>();
				tempList.add(bestSolution);
				Map<File, List<Solution>> tempMap = new LinkedHashMap<>();
				tempMap.put(outputFile, tempList);
				SpreadsheetGenerator.writeSpreadsheet(outputFile, tempMap, true);
				om.writeValue(outputFileJSON, bestSolution);
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