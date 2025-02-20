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
/**
 * Set numClusters equal to TrueNumClusters
 * @author rickw
 *
 */
public class BatchRunPDCGM_TrueNumClusters {

	private static Logger log = LoggerFactory.getLogger(BatchRunPDCGM_TrueNumClusters.class);

	public static void main(String[] args) throws IloException, IOException {
		// Retrieve input argument.
		String fileName = args[0];
		
		// Retrieve input argument.
		int size = Integer.valueOf(args[1]);
				
		BasicConfigurator.configure();
		File f = new File("data/Thrun"+size+"/"+fileName+".txt");
		String date_file = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date());
		date_file += "_PDCGM_TrueNumClusters_"+fileName+size;
		File outputDir = new File("data/batchsolutions/"+date_file);
		outputDir.mkdirs();
		ObjectMapper om = new ObjectMapper();

		// Main parameters
		HierarchicalHeuristicType hierarchicalHeuristicType = HierarchicalHeuristicType.hmc;
		SolverType solverType = SolverType.Colgen_BranchAndPrice;
		PricingProblemStrategy pricingProblemStrategy = PricingProblemStrategy.None;
		
		Map<File, List<Solution>> solutions = new LinkedHashMap<>();
		
		List<Solution> solutionList = new ArrayList<>();
		
		log.info("Reading file {}", f);
		Instance instance = Instance.readFromTxt(f, 1);
		if(instance.getTrueNumClusters()==1) {
			return;
		}
		instance.setNumClusters(instance.getTrueNumClusters());
		
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
		
		File outputFile = new File(outputDir+"/java_"+date_file+".xlsx");
		File outputFileJSON = new File(outputDir+"/java_"+date_file+".json");
		List<Solution> tempList = new ArrayList<>();
		tempList.add(bestSolution);
		Map<File, List<Solution>> tempMap = new LinkedHashMap<>();
		tempMap.put(outputFile, tempList);
		SpreadsheetGenerator.writeSpreadsheet(outputFile, tempMap, true);
		om.writeValue(outputFileJSON, bestSolution);
		
		solutions.put(f, solutionList);
	}
}