package cl.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.CombinatoricsUtils;

import cl.data.Cluster;
import cl.data.Instance;
/**
 * Calculate Adjusted Rand Index (ARI), based on:
 * - https://www.mathworks.com/matlabcentral/fileexchange/49908-adjusted-rand-index
 * - https://github.com/areslp/matlab/blob/master/code_cospectral/RandIndex.m
 * - https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html
 */
public class AdjustedRandIndex {
	public static List<Integer> convertClusterToPartition(Instance instance, List<Cluster> clusters) {
		Integer[] data = new Integer[instance.getNumPoints()];
		int count = 0;
		for(Cluster cluster: clusters) {
			for(int i: cluster.getPointId()) {
				data[i] = count; 
			}
			count++;
		}
		return Arrays.asList(data);
	}
	
	public static double calcARI(List<Integer> p1, List<Integer> p2) {
		List<Integer> unique1 = findUnique(p1);
		List<Integer> unique2 = findUnique(p2);
		
		double[][] contingencyTable = createContingencyTable(p1, p2, unique1.size(), unique2.size());
		RealMatrix contingencyMatrix = new Array2DRowRealMatrix(contingencyTable);
		
		int n = (int) MatrixOperations.sumAllElements(contingencyMatrix);
		RealMatrix rowTotals = MatrixOperations.sumElementsHorizontally(contingencyMatrix);
		RealMatrix columnTotals = MatrixOperations.sumElementsVertically(contingencyMatrix);
		
		double numerator = 0;
		for(int i = 0; i < unique1.size(); i++) {
			for(int j = 0; j < unique2.size(); j++) {
				numerator += nchoosek2((int) contingencyMatrix.getEntry(i, j), 2);
			}
		}
		
		double rowBin = 0;
		for(int i = 0; i < unique1.size(); i++) {
			rowBin += nchoosek2((int) rowTotals.getEntry(i, 0), 2);
		}
		double colBin = 0;
		for(int i = 0; i < unique2.size(); i++) {
			colBin += nchoosek2((int) columnTotals.getEntry(i, 0), 2);
		}
		
		numerator -= (rowBin * colBin) / nchoosek2(n, 2);
		
		double denominator = 0.5 * (rowBin + colBin) - (rowBin * colBin) / nchoosek2(n, 2);
		
		double ARI = numerator / denominator;
		return ARI;
	}
	
	 private static int nchoosek2(int a, int b) {
        if (a > 1) {
            return (int) Math.round(CombinatoricsUtils.binomialCoefficient(a, b));
        } else {
            return 0;
        }
    }
		 
	private static List<Integer> findUnique(List<Integer> partition) {
		Set<Integer> temp = new LinkedHashSet<>();
		
		for(int i: partition) {
			temp.add(i);
		}
		
		List<Integer> result = new ArrayList<>(temp);
		Collections.sort(result);
		
		// We always start the index at 0
		if(result.get(0)!=0) {
			List<Integer> result2 = new ArrayList<>();
			for(int i = 0; i < result.size(); i++) {
				result2.add(i);
			}
			return result2;
		}
		return result;
	}
	
	private static double[][] createContingencyTable(List<Integer> partition1, List<Integer> partition2, int num1, int num2) {
		double[][] result = new double[num1][num2];
		
		for(int i = 0; i < partition1.size(); i++) {
			result[partition1.get(i)][partition2.get(i)] += 1;
		}
		return result;
	}
}