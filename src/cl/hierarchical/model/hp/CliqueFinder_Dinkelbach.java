package cl.hierarchical.model.hp;

import java.util.List;

import cl.data.GlobalParam;
import cl.data.Instance;
import cl.hierarchical.model.extended.BranchInformationExtra;
import cl.hierarchical.model.extended.DualInformation;
import cl.hierarchical.model.hp.qp.Dinkelbach_QP_Hansen;
import cl.hierarchical.model.hp.qp.Dinkelbach_QP_VNS;
import cl.hierarchical.model.hp.qp.QP_Solver;
import cl.util.Pair;
import cl.util.Triple;

/**
 * Dinkelbach's algorithm, based on the description in the following papers:
 * - An interior point algorithm for minimum sum-of-squares clustering by Olivier du Merle, Pierre Hansen, Brigitte Jaumard, and Nenad Mladenovic (1999)
 * - An improved column generation algorithm for minimum sum-of-squares clustering by Daniel Aloise, Pierre Hansen, and Leo Liberti (2010) 
 */
public class CliqueFinder_Dinkelbach implements CliqueFinder {

	private Instance instance;
	private int level;
	private DualInformation duals;
	private List<Integer> Gi;
	private BranchInformationExtra bix;

	private Triple<List<Integer>, Double, Double> clique;

	public CliqueFinder_Dinkelbach(Instance instance, int level, DualInformation duals, List<Integer> Gi, BranchInformationExtra bix)
	{
		this.instance = instance;
		this.level = level;
		this.duals = duals;
		this.Gi = Gi;
		this.bix = bix;
	}

	@Override
	/**
	 * Keep using VNS
	 * When we have to prove optimality switch to Algorithm 2
	 */
	public void solve() {
		double obj_QP = 0;
		double obj_HP = 0;
		boolean cont = true;
		boolean boolVNS = true;
		boolean[] startSolution = null;
		double prevObj_QP = Double.NaN;
		while(cont) {
			QP_Solver qp;
			if(boolVNS)
			{
				GlobalParam.COUNTER_VNS++;
				qp = new Dinkelbach_QP_VNS(instance, level, duals, Gi, bix, obj_HP+GlobalParam.EPSILON5);
				if(startSolution!=null) {
					qp.setStartSolution(startSolution); // not all methods need a start solution
				}
			}
			else 
			{
				GlobalParam.COUNTER_QP++;
				qp = new Dinkelbach_QP_Hansen(instance, level, duals, Gi, bix, obj_HP+GlobalParam.EPSILON5);
			}
			qp.solve();
			Pair<Double, Double> objs = qp.getObjValue();
			
			obj_QP = objs.first;
			if(obj_QP>-GlobalParam.EPSILON3 || (!Double.isNaN(prevObj_QP) && Math.abs(prevObj_QP - obj_QP)<GlobalParam.EPSILON5))
			{
				if(!GlobalParam.ONLY_USE_VNS && boolVNS)
				{
					boolVNS = false;
					startSolution = null;
					continue;
				}
				obj_HP = objs.second;

				// Calculate actual objective of the hyperbolic program
				double actual_obj_HP = obj_HP - duals.getDualsMaxCluster()[level]; 
				clique = new Triple<>(qp.getClique(), actual_obj_HP, actual_obj_HP);
				break;
			}
			obj_HP = objs.second;
			prevObj_QP = obj_QP;
			if(boolVNS)
			{
				startSolution = qp.getSolution();
			}
		}
	}

	@Override
	public Triple<List<Integer>, Double, Double> getClique() {
		return clique;
	}
}