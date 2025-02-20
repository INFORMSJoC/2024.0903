package cl.hierarchical.model.hp.qp;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import cl.data.GlobalParam;
import cl.data.Instance;
import cl.data.ObjectiveWeightFunction;
import cl.hierarchical.model.extended.BranchInformationExtra;
import cl.hierarchical.model.extended.DualInformation;
import cl.util.Pair;

/**
 * Variable neighbourhood solver to solve QUBO, based on the description in the following paper:
 * - An interior point algorithm for minimum sum-of-squares clustering by Olivier du Merle, Pierre Hansen, Brigitte Jaumard, and Nenad Mladenovic (1999)
 */
public class Dinkelbach_QP_VNS implements QP_Solver {

	private Instance instance;
	private int level;
	private DualInformation duals;
	private List<Integer> Gi;
	private BranchInformationExtra bix;

	private boolean[] startSolution = null;
	private double prevObj;
	private boolean[] bestX;

	private int k_max; //Equal to the length of Gi

	private double coeffMatrix[][]; // Store coefficients for more efficient computation

	public Dinkelbach_QP_VNS(Instance instance, int level, DualInformation duals, List<Integer> Gi, BranchInformationExtra bix, double prevObj) 
	{
		this.instance = instance;
		this.level = level;
		this.duals = duals;
		this.Gi = Gi;
		this.bix = bix;

		this.prevObj = prevObj;
		this.k_max = Math.min(Gi.size(), GlobalParam.VNS_MAX_NUM_NEIGHBOURHOODS);

		initCoeffMatrix();
	}

	public void solve()
	{
		// Starting point
		if(startSolution!=null)
		{
			bestX = startSolution;
		}
		else
		{
			bestX = new boolean[Gi.size()];
		}
		double obj = calcObjectiveQP(bestX).first;

		// Local search: greedy algorithm on Neighbourhood N1
		Pair<Double, boolean[]> pair = localSearch(bestX, obj);
		obj = pair.first;
		bestX = pair.second;
		
		// Only define it once
		ArrayList<Integer> indexList = new ArrayList<Integer>();
		for(int i = 0; i < Gi.size(); i++)
		{
			indexList.add(i);
		}

		// Start search procedure
		for(int k = 1; k <= k_max; k++)
		{
			// Shaking: generate random point from Neighbourhood Nk
			boolean[] xNew = bestX.clone();
			
			double objNew = obj;
			
			Collections.shuffle(indexList, GlobalParam.RANDOM);
			for(int i = 0; i < k; i++)
			{
				int ind = indexList.get(i);
				xNew[ind] = !xNew[ind];
				objNew = recalculateObjectiveQP(xNew, ind, objNew);				
			}

			// Local search: greedy algorithm on Neighbourhood N1
			pair = localSearch(xNew, objNew);
			objNew = pair.first;
		
			// Move or not: if better solution is found, set k=1. If not, set k++.
			if(objNew < obj - GlobalParam.EPSILON5)
			{
				obj = objNew;
				bestX = pair.second;
				k = 0; // Due to the for loop (k++)
			}
		}
	}

	public void setStartSolution(boolean[] startSolution)
	{
		this.startSolution = startSolution; 
	}

	public boolean[] getSolution()
	{
		return bestX;
	}

	private Pair<Double, boolean[]> localSearch(boolean[] xNew, double oldObj)
	{
		double objNew = oldObj; 	
		
		for(int i = 0; i < Gi.size(); i++)
		{
			xNew[i] = !xNew[i];
			double objNewTemp = recalculateObjectiveQP(xNew, i, objNew);
			if(objNewTemp < objNew)
			{
				// Update the best found solution
				objNew = objNewTemp;
			}
			else
			{
				// Revert changes
				xNew[i] = !xNew[i];
			}
		}
		return new Pair<>(objNew, xNew);
	}

	private void initCoeffMatrix()
	{
		double dist = 0;
		if(!bix.isBranch()) {
			Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
			dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
		}
		
		coeffMatrix = new double[Gi.size()][Gi.size()];
		for(int i = 0; i < Gi.size(); i++) {
			for(int j = i+1; j < Gi.size(); j++) {
				if(!bix.isBranch()) {
					if(Math.sqrt(dist) * instance.getRegularDistanceMatrix()[Gi.get(i)][Gi.get(j)] - duals.getDualsSqrtPoint()[level][Gi.get(i)] - duals.getDualsSqrtPoint()[level][Gi.get(j)]>0) {
						coeffMatrix[i][j] = GlobalParam.BIGM8;
					}
					else {
						coeffMatrix[i][j] = dist * instance.getSquaredDistanceMatrix()[Gi.get(i)][Gi.get(j)] - duals.getDualsPoint()[level][Gi.get(i)] - duals.getDualsPoint()[level][Gi.get(j)];
					}
				}
				else {
					if(bix.getReducedInstance().getRegularDistanceMatrix()[Gi.get(i)][Gi.get(j)]
							- bix.getBranchMappingWeights().get(Gi.get(j)) * bix.getReducedDuals().getDualsSqrtPointSingle()[Gi.get(i)] - bix.getBranchMappingWeights().get(Gi.get(i)) * bix.getReducedDuals().getDualsSqrtPointSingle()[Gi.get(j)]>0) {
						coeffMatrix[i][j] = GlobalParam.BIGM8;
					}
					else {
						coeffMatrix[i][j] = bix.getReducedInstance().getSquaredDistanceMatrix()[Gi.get(i)][Gi.get(j)]
								- bix.getBranchMappingWeights().get(Gi.get(j)) * bix.getReducedDuals().getDualsPointSingle()[Gi.get(i)] - bix.getBranchMappingWeights().get(Gi.get(i)) * bix.getReducedDuals().getDualsPointSingle()[Gi.get(j)];
					}
				}
			}
		}
	}
	
	/**
	 * We changed i, what is the new objective value of the QP?
	 * @param x
	 * @param i
	 * @param prevObj
	 * @return
	 */
	private double recalculateObjectiveQP(boolean[] x, int vertex, double oldObj) {
		if(oldObj>=GlobalParam.BIGM8) {
			return calcObjectiveQP(x).first;
		}
		
		// Numerator
		double numerator = 0;
		for(int i = 0; i < Gi.size(); i++) {
			if(i==vertex) {
				continue;
			}
			
			if(x[i]) {
				if(i < vertex) {
					if (coeffMatrix[i][vertex]>=GlobalParam.BIGM8)
					{
						return GlobalParam.BIGM8;
					}
					numerator += coeffMatrix[i][vertex];
				}
				else {
					if (coeffMatrix[vertex][i]>=GlobalParam.BIGM8)
					{
						return GlobalParam.BIGM8;
					}
					numerator += coeffMatrix[vertex][i];
				}
			}
		}
		
		// Denominator
		double obj_QP = 0;
		if(!bix.isBranch()) {
			numerator -= duals.getDualsPoint()[level][Gi.get(vertex)];
			obj_QP -= prevObj;
		}
		else {
			double w = bix.getBranchMappingWeights().get(Gi.get(vertex));
			numerator -= w * bix.getReducedDuals().getDualsPointSingle()[Gi.get(vertex)] - bix.getReducedInstance().getSquaredDistanceMatrix()[Gi.get(vertex)][Gi.get(vertex)];
			obj_QP -= w *prevObj;
		}
		
		// Update QP
		obj_QP += numerator;
		
		if(x[vertex]) {
			return oldObj + obj_QP;
		}
		return oldObj - obj_QP;
	}

	private Pair<Double, Double> calcObjectiveQP(boolean[] x)
	{
		double numerator = 0;
		for(int i = 0; i < Gi.size(); i++) {
			for(int j = i+1; j < Gi.size(); j++) {
				if(x[i] && x[j])
				{
					if (coeffMatrix[i][j]>=GlobalParam.BIGM8)
					{
						return new Pair<>(GlobalParam.BIGM8, 0d);
					}
					numerator += coeffMatrix[i][j];
				}
			}
		}
		double obj_QP = 0;
		double denominator = 0;
		for(int i = 0; i < Gi.size(); i++) {
			if(x[i])
			{
				if(!bix.isBranch()) {
					numerator -= duals.getDualsPoint()[level][Gi.get(i)];
					obj_QP -= prevObj;
					denominator ++;
				}
				else {
					double w = bix.getBranchMappingWeights().get(Gi.get(i));
					numerator -= w * bix.getReducedDuals().getDualsPointSingle()[Gi.get(i)] - bix.getReducedInstance().getSquaredDistanceMatrix()[Gi.get(i)][Gi.get(i)];
					obj_QP -= w *prevObj;
					denominator+= w;
				}
			}
		}

		double obj_HP = numerator / denominator;
		obj_QP += numerator;
		return new Pair<>(obj_QP, obj_HP);
	}

	public Pair<Double, Double> getObjValue()
	{
		return calcObjectiveQP(bestX);
	}

	public List<Integer> getClique() {
		List<Integer> result = new ArrayList<>();
		for(int i = 0; i < Gi.size(); i++) {
			if(bestX[i]) {
				result.add(Gi.get(i));
			}
		}
		return result;
	}
}
