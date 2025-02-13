package cl.data;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class ObjectiveWeightFunction {
	private static List<Double> weightsProportionalToLevel;
	
	/**
	 * 
	 * Return the value of: sum_{d=h}^k w_{ghk}
	 * 
	 * Note that:
	 * k - h + 1
	 * So k>h
	 * 
	 * @param k
	 * @param h
	 * @return
	 */
	public static double getWeight(int k, int h) {
		if(GlobalParam.OBJECTIVE_WEIGHT==ObjectiveWeight.noWeight) {
			// w_{ghk} = 1
			// w_d = 1
			return k - h + 1; // = sum_{d=h}^k w_{d} = sum_{d=h}^k 1 = k - h + 1
		}
		else if(GlobalParam.OBJECTIVE_WEIGHT==ObjectiveWeight.noDoubleCounting) {
			// w_{ghk} = 1/(k - h + 1)
			// w_d = 1/(k - h + 1)
			return 1; // = sum_{d=h}^k w_{d} = sum_{d=h}^k 1 /(k - h + 1) = 1
		}
		else if(GlobalParam.OBJECTIVE_WEIGHT==ObjectiveWeight.proportionalToLevel) {
			// w_{ghk} = sum_{d=h}^k o_1 / o_d
			// w_d = o_1 / o_d
			return getProportionalToLevel(k, h); // sum_{d=h}^k w_{d}  = sum_{d=h}^k o_1 / o_d
		}
		else {
			throw new IllegalArgumentException("OBJECTIVE_WEIGHT not defined");
		}
	}
	
	public static void setWeightsProportionalToLevel(Instance instance, Map<Integer, List<Cluster>> hierarchicalClusters) {
		if(weightsProportionalToLevel!=null) {
			// Only overwrite if there were no previous weights
			return;
		}
		// Calculate cost per level
		List<Double> objectivesPerLevel = new ArrayList<>();
		for(Entry<Integer, List<Cluster>> entry: hierarchicalClusters.entrySet()) {
			double cost = 0;
			for(Cluster cluster: entry.getValue()) {
				cost += instance.computeClusterCost(cluster); 
			}
			objectivesPerLevel.add(cost);
		}
		
		// Define weights
		weightsProportionalToLevel = new ArrayList<>();
		for(int i = 0; i < objectivesPerLevel.size(); i++) {
			weightsProportionalToLevel.add(objectivesPerLevel.get(0) /  objectivesPerLevel.get(i));
			System.out.println("Weights "+ (objectivesPerLevel.get(0) /  objectivesPerLevel.get(i)));
		}
	}

	public static double getProportionalToLevel(int k, int h) {
		double sum = 0;
		for(int i = h; i <= k; i++) {
			if(weightsProportionalToLevel==null) {
				sum++;
			}
			else {
				sum += weightsProportionalToLevel.get(i);
			}
		}
		return sum;
	}
	
	public static double getProportionalToLevel(int k) {
		if(weightsProportionalToLevel==null) {
			return 1;
		}
		return weightsProportionalToLevel.get(k);
	}
}
