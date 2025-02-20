package cl.util.run;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class InstanceConverter {
	
	public static void main(String[] args) {
		File inputDir = new File("data/instance_converter/input");
		File outputDir = new File("data/instance_converter/output");
		outputDir.mkdirs();
		
		int numObservationLimit = 50;
		
		for (File inputFile : inputDir.listFiles()) {
			String fileName = inputFile.getName().replace(".txt", "");
			File outputFile = new File(outputDir+"/"+fileName+"_"+numObservationLimit+".txt");
			createSmallerInstances(inputFile, outputFile, numObservationLimit);
		}
	}

	/**
	 * If applicable: Get the same percentage of observations from each cluster
	 * @param inputFile
	 * @param outputFile
	 * @param numObservationLimit
	 */
	public static void createSmallerInstances(File inputFile, File outputFile, int numObservationLimit) {
		try {
			// Read original instance
			Scanner sc = new Scanner(inputFile);

			String header = sc.nextLine();
			String[] headers = header.split(" ");
			int numPoints = Integer.parseInt(headers[0]);
			int numDimensions = Integer.parseInt(headers[1]);
			
			if(numPoints <= numObservationLimit) {
				sc.close();
				throw new IllegalArgumentException("Original number of points "+ numPoints + " is equal or lower than the limit "+ numObservationLimit);
			}

			boolean trueClusters = false;
			int trueNumClusters = 0;
			if(headers.length==3) {
				trueClusters = true;
				trueNumClusters = Integer.parseInt(headers[2]);
			}

			Map<Integer, List<String>> observationsPerCluster = new LinkedHashMap<>();
			for(int i = 0; i < numPoints; i++) {
				String s = sc.nextLine();
				
				int clusterID = 0;
				if(trueClusters) {
					clusterID = s.charAt(s.length()-1);
				}
				
				if(!observationsPerCluster.containsKey(clusterID)) {
					observationsPerCluster.put(clusterID, new ArrayList<>());
				}
				observationsPerCluster.get(clusterID).add(s);
			}
			sc.close();

			// Write new instance
			String result = numObservationLimit + " "+ numDimensions;
			if(trueClusters) {
				result += " " + trueNumClusters;
			}
			result += "\n";
			
			int countObservations = 0;
			int countCluster = 0;
			for(Entry<Integer, List<String>> entry: observationsPerCluster.entrySet()) {
				// Find out how many observations we select per cluster
				List<String> observations = entry.getValue();
				int numToSelect = (int) Math.round(numObservationLimit / ((double) numPoints) * observations.size());
				if(countCluster==observationsPerCluster.size()-1) { // Due to rounding errors
					numToSelect = numObservationLimit - countObservations;
				}
				
				List<Integer> range = IntStream.range(0, numPoints).boxed().collect(Collectors.toList());
				Collections.shuffle(range);
				
				for(int i = 0; i < numToSelect; i++) {
					result += observations.get(i) + "\n";
				}
				
				countObservations += numToSelect;
				countCluster++;
			}
						
			PrintWriter out = new PrintWriter(outputFile);
			out.append(result);
			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}	
	}
}
