package cl.util.run;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import cl.data.GlobalParam;
import cl.data.Solution;

public class SpreadsheetGenerator {

	/**
	 * Create the spreadsheet
	 * @param output Output file
	 * @param solutions Input solutions
	 * @throws IOException
	 */
	public static void writeSpreadsheetAverage(File output, Map<File, List<Solution>> solutions, boolean columnGeneration) throws IOException {
		try (XSSFWorkbook wb = new XSSFWorkbook(); OutputStream os = new FileOutputStream(output)) {
			Sheet s = wb.createSheet();
			int rowIndex = 0;

			// Create Header
			Row headerRow = s.createRow(rowIndex++);

			if(columnGeneration) {
				writeCells(headerRow, "folderName", "maxLevel", "distanceType", "distanceFunction", "solverType", "hierarchicalHeuristicType", "pricingProblemStrategy", "numFeasible", "numTimeLimitReached", "numInstances", "AvgOriginalObjective", "AvgStartObjective", "AvgLP", "AvgObjective","AvgGapObjectiveLP",
						"AvgActualTotalTime (s)", "AvgHeuristicTime (s)", "AvgRootNodeTime (s)", "AvgMasterTime (s)", "AvgPricingTime (s)",
						"AvgNumIter", "AvgNumColumns", "AvgNumVNS" , "AvgNumQP" , "AvgNumBranchAndBound", "AvgNumKnapsack", "AvgNumBranchAndPrice", "AvgNumCuttingPlanes", "AvgNumEnumerate",
						"AvgNumTimesStartEqualsOptimal", "AvgNumStartSolutions");
			}
			else {
				writeCells(headerRow, "folderName", "maxLevel", "distanceType", "distanceFunction", "solverType", "numFeasible", "numTimeLimitReached", "AvgOriginalObjective", "AvgStartObjective", "AvgObjective", "StdObjective", "AvgComputationTime (s)", "StdComputationTime",
						"AvgNumIter", "AvgNumColumns" );	
			}

			for (Entry<File, List<Solution>> e : solutions.entrySet()) { //Go through all the solutions
				String folderName = e.getKey().getName();

				if(columnGeneration) {
					int numFeasible = 0;
					int numTimeLimitReached = 0;
					SummaryStatistics ssOriginalObjective = new SummaryStatistics();
					SummaryStatistics ssStartObjective = new SummaryStatistics();
					SummaryStatistics ssLP = new SummaryStatistics();
					SummaryStatistics ssObjective = new SummaryStatistics();
					SummaryStatistics ssGap = new SummaryStatistics();
					
					SummaryStatistics ssActualTotalTime = new SummaryStatistics();
					SummaryStatistics ssHeuristicTime = new SummaryStatistics();
					SummaryStatistics ssRootNodeTime = new SummaryStatistics();
					SummaryStatistics ssMasterTime = new SummaryStatistics();
					SummaryStatistics ssPricingTime = new SummaryStatistics();
					
					SummaryStatistics ssNumIter = new SummaryStatistics();
					SummaryStatistics ssNumColumns = new SummaryStatistics();
					SummaryStatistics ssNumQP = new SummaryStatistics();
					SummaryStatistics ssNumVNS = new SummaryStatistics();
					SummaryStatistics ssNumBranchAndBound = new SummaryStatistics();
					SummaryStatistics ssNumKnapsack = new SummaryStatistics();
					SummaryStatistics ssNumBranchAndPrice = new SummaryStatistics();
					SummaryStatistics ssNumCuttingPlanes = new SummaryStatistics();
					SummaryStatistics ssNumEnumerate = new SummaryStatistics();
					
					SummaryStatistics ssNumStartEqualsOptimal = new SummaryStatistics();
					SummaryStatistics ssNumStartSolutions = new SummaryStatistics();
					
					for(Solution sol: e.getValue()) {
						if(sol.isFeasible()) {
							numFeasible++;
						}
						if(sol.isTimeLimitReached()) {
							numTimeLimitReached++;
						}
						ssOriginalObjective.addValue(sol.getOriginalObjective());
						ssStartObjective.addValue(sol.getBestStartObjective());
						ssLP.addValue(sol.getBestLP());
						ssObjective.addValue(sol.getObjective());
						ssGap.addValue((sol.getObjective() - sol.getBestLP()) / sol.getBestLP() * 100);
						
						ssActualTotalTime.addValue((sol.getHeuristicTime() + sol.getTotalTime())/1000d);
						ssHeuristicTime.addValue(sol.getHeuristicTime()/1000d);
						ssRootNodeTime.addValue(sol.getRootNodeTime()/1000d);
						ssMasterTime.addValue(sol.getMasterTime()/1000d);
						ssPricingTime.addValue(sol.getPricingTime()/1000d);
						
						ssNumIter.addValue(sol.getNumIter());
						ssNumColumns.addValue(sol.getNumColumns());
						ssNumVNS.addValue(sol.getNumVNS());
						ssNumQP.addValue(sol.getNumQP());
						ssNumBranchAndBound.addValue(sol.getNumBranchAndBound());
						ssNumKnapsack.addValue(sol.getNumKnapsack());
						ssNumBranchAndPrice.addValue(sol.getNumBranchAndPrice());
						ssNumCuttingPlanes.addValue(sol.getNumCuttingPlanes());
						ssNumEnumerate.addValue(sol.getNumEnumerate());
						
						int numTimesStartEqualsOptimal = 0;
						for(double d: sol.getStartObjectives()) {
							if(sol.getObjective() - GlobalParam.EPSILON5 < d && sol.getObjective() + GlobalParam.EPSILON5 > d) {
								numTimesStartEqualsOptimal++;
							}
						}
						ssNumStartEqualsOptimal.addValue(numTimesStartEqualsOptimal);
						ssNumStartSolutions.addValue(sol.getStartObjectives().size());
					}

					Row row = s.createRow(rowIndex++);

					writeCells(row,
							folderName,
							e.getValue().get(0).getInstance().getNumClusters(),
							GlobalParam.DISTANCE_TYPE,
							GlobalParam.DISTANCE_FUNCTION,
							e.getValue().get(0).getSolverType(),
							e.getValue().get(0).getHierarchicalHeuristicType(),
							e.getValue().get(0).getPricingProblemStrategy(),
							numFeasible,
							numTimeLimitReached,
							e.getValue().size(),
							ssOriginalObjective.getMean(),
							ssStartObjective.getMean(),
							ssLP.getMean(),
							ssObjective.getMean(),
							ssGap.getMean(),
							ssActualTotalTime.getMean(),
							ssHeuristicTime.getMean(),
							ssRootNodeTime.getMean(),
							ssMasterTime.getMean(),
							ssPricingTime.getMean(),
							ssNumIter.getMean(),
							ssNumColumns.getMean(),
							ssNumVNS.getMean(),
							ssNumQP.getMean(),
							ssNumBranchAndBound.getMean(),
							ssNumKnapsack.getMean(),
							ssNumBranchAndPrice.getMean(),
							ssNumCuttingPlanes.getMean(),
							ssNumEnumerate.getMean(),
							ssNumStartEqualsOptimal.getMean(),
							ssNumStartSolutions.getMean()
							);
				}
				else {
					SummaryStatistics ssOriginalObjective = new SummaryStatistics();
					SummaryStatistics ssObjective = new SummaryStatistics();
					SummaryStatistics ssStartObjective = new SummaryStatistics();
					SummaryStatistics ssComputationTime = new SummaryStatistics();
					SummaryStatistics ssNumIter = new SummaryStatistics();
					SummaryStatistics ssNumColumns = new SummaryStatistics();
					int numFeasible = 0;
					int numTimeLimitReached = 0;
					for(Solution sol: e.getValue()) {
						if(sol.isFeasible()) {
							numFeasible++;
						}
						if(sol.isTimeLimitReached()) {
							numTimeLimitReached++;
						}
						ssOriginalObjective.addValue(sol.getOriginalObjective());
						ssStartObjective.addValue(sol.getBestStartObjective());
						ssObjective.addValue(sol.getObjective());
						ssComputationTime.addValue(sol.getTotalTime()/1000d);
						ssNumIter.addValue(sol.getNumIter());
						ssNumColumns.addValue(sol.getNumColumns());
					}

					Row row = s.createRow(rowIndex++);

					writeCells(row,
							folderName,
							e.getValue().get(0).getInstance().getNumClusters(),
							GlobalParam.DISTANCE_TYPE,
							GlobalParam.DISTANCE_FUNCTION,
							e.getValue().get(0).getSolverType(),
							numFeasible,
							numTimeLimitReached,
							ssOriginalObjective.getMean(),
							ssStartObjective.getMean(),
							ssObjective.getMean(),
							ssObjective.getStandardDeviation(),
							ssComputationTime.getMean(),
							ssComputationTime.getStandardDeviation(),
							ssNumIter.getMean(),
							ssNumColumns.getMean()
							);
				}
			}
			wb.write(os);
		}
	}

	/**
	 * Create the spreadsheet
	 * @param output Output file
	 * @param solutions Input solutions
	 * @throws IOException
	 */
	public static void writeSpreadsheet(File output, Map<File, List<Solution>> solutions, boolean columnGeneration) throws IOException {
		try (XSSFWorkbook wb = new XSSFWorkbook(); OutputStream os = new FileOutputStream(output)) {
			Sheet s = wb.createSheet();
			int rowIndex = 0;

			// Create Header
			Row headerRow = s.createRow(rowIndex++);
			if(columnGeneration) {
				writeCells(headerRow, "folderName", "fileName", "maxLevel", "distanceType", "distanceFunction", "solverType", "hierarchicalHeuristicType", "pricingProblemStrategy", "trueNumClusters", "ARI (max level)", "feasible", "timeLimitReached", "originalObjective", "startObjective", "LP", "objective", "gapObjectiveLP",
						"actualTotalTime (s)", "heuristicTime (s)", "rootNodeTime (s)", "masterTime (s)", "pricingTime (s)",
						"numIter", "numColumns", "numVNS", "numQP" , "numBranchAndBound", "numKnapsack", "numBranchAndPrice", "numCuttingPlanes", "numEnumerate",
						"numTimesStartEqualsOptimal", "numStartSolutions", "startSolutions");
			}
			else {
				writeCells(headerRow, "folderName", "fileName", "maxLevel", "distanceType", "distanceFunction", "solverType", "trueNumClusters", "ARI (max level)", "feasible", "timeLimitReached", "originalObjective", "startObjective", "objective", "computationTime (s)",
						"numIter", "numColumns");
			}

			for (Entry<File, List<Solution>> e : solutions.entrySet()) { //Go through all the solutions
				String folderName = e.getKey().getName();
				for(Solution sol: e.getValue()) {
					Row row = s.createRow(rowIndex++);

					if(columnGeneration) {
						double actualTotalTime = sol.getHeuristicTime() + sol.getTotalTime();
						int numTimesStartEqualsOptimal = 0;
						for(double d: sol.getStartObjectives()) {
							if(sol.getObjective() - GlobalParam.EPSILON5 < d && sol.getObjective() + GlobalParam.EPSILON5 > d) {
								numTimesStartEqualsOptimal++;
							}
						}

						double gap = (sol.getObjective() - sol.getBestLP()) / sol.getBestLP() * 100;
							
						writeCells(row,
								folderName,
								sol.getInstance().getInstanceName(),
								sol.getInstance().getNumClusters(),
								GlobalParam.DISTANCE_TYPE,
								GlobalParam.DISTANCE_FUNCTION,
								sol.getSolverType(),
								sol.getHierarchicalHeuristicType(),
								sol.getPricingProblemStrategy(),
								sol.getInstance().getTrueNumClusters(),
								sol.getARI(),
								sol.isFeasible(),
								sol.isTimeLimitReached(),
								sol.getOriginalObjective(),
								sol.getBestStartObjective(),
								sol.getBestLP(),
								sol.getObjective(),
								gap,
								actualTotalTime/1000d,
								sol.getHeuristicTime()/1000d,
								sol.getRootNodeTime()/1000d,
								sol.getMasterTime()/1000d,
								sol.getPricingTime()/1000d,
								sol.getNumIter(),
								sol.getNumColumns(),
								sol.getNumVNS(),
								sol.getNumQP(),
								sol.getNumBranchAndBound(),
								sol.getNumKnapsack(),
								sol.getNumBranchAndPrice(),
								sol.getNumCuttingPlanes(),
								sol.getNumEnumerate(),
								numTimesStartEqualsOptimal,
								sol.getStartObjectives().size(),
								Arrays.toString(sol.getStartObjectives().toArray())
								);
					}
					else {
						writeCells(row,
								folderName,
								sol.getInstance().getInstanceName(),
								sol.getInstance().getNumClusters(),
								GlobalParam.DISTANCE_TYPE,
								GlobalParam.DISTANCE_FUNCTION,
								sol.getSolverType(),
								sol.getInstance().getTrueNumClusters(),
								sol.getARI(),
								sol.isFeasible(),
								sol.isTimeLimitReached(),
								sol.getOriginalObjective(),
								sol.getBestStartObjective(),
								sol.getObjective(),
								sol.getTotalTime()/1000d,
								sol.getNumIter(),
								sol.getNumColumns()
								);
					}
				}
			}
			wb.write(os);
		}
	}

	/**
	 * Create a row in spreadsheet
	 * @param row Row number
	 * @param objects Objects to write in excel file
	 */
	private static void writeCells(Row row, Object... objects) {
		for (int i=0; i < objects.length; i++) {
			Cell c = row.createCell(i);
			Object o = objects[i];
			if (o == null) {
				continue;
			}
			else if (o instanceof Number) {
				Number n = (Number) o;
				c.setCellValue(n.doubleValue());
			}
			else {
				c.setCellValue(o.toString());
			}
		}
	}
}
