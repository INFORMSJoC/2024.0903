package cl.hierarchical.model.hp;

import java.util.ArrayList;
import java.util.List;

import cl.data.GlobalParam;
import cl.data.Instance;
import cl.hierarchical.model.extended.BranchInformation;
import cl.hierarchical.model.extended.BranchInformationExtra;
import cl.hierarchical.model.extended.DualInformation;
import cl.util.Pair;
import cl.util.Triple;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;

public class CliqueFinder_Knapsack implements CliqueFinder {
	private Instance instance;
	private int level;
	private int currentVertex;
	private DualInformation duals;
	private BranchInformation bi;
	private BranchInformationExtra bix;

	private IloCplex model;
	private IloObjective objective;

	private IloNumVar[] vars;

	public CliqueFinder_Knapsack(Instance instance, int level, int currentVertex, DualInformation duals, BranchInformation bi, BranchInformationExtra bix) throws IloException {
		this.instance = instance;
		this.level = level;
		this.currentVertex = currentVertex;
		this.duals = duals;
		this.bi = bi;
		this.bix = bix;

		this.model = new IloCplex();
		this.model.setOut(null);	
		initModel();
	}

	private void initModel() throws IloException {
		initVariables();
		initForbiddenConstraints();
		initObjective();
	}

	private void initVariables() throws IloException {
		vars = new IloNumVar[bix.getBranchMapping().size()];
		for(int i = 0; i < bix.getBranchMapping().size(); i++) {
			vars[i] = model.boolVar();
			if(bix.getBranchMapping().get(i).contains(currentVertex)) {
				// We always have to select currentVertex
				vars[i].setLB(1); 
			}
		}
	}

	private void initForbiddenConstraints() throws IloException {
		for(int i = 0; i < bix.getBranchMapping().size(); i++) {
			for(int j = i+1; j < bix.getBranchMapping().size(); j++) {
				List<Integer> list1 = bix.getBranchMapping().get(i);
				List<Integer> list2 = bix.getBranchMapping().get(j);

				// Check if list1 and list2 are forbidden
				boolean forbidden = false;
				outerloop: for(Integer m: list1) {
					for(Integer n: list2) {
						Pair<Integer, Integer> pair;
						if(m<n) {
							pair = new Pair<>(m, n);
						}
						else {
							pair = new Pair<>(n, m);
						}
						if(bi.getForbidden(level).contains(pair)) {
							forbidden = true;
							break outerloop;
						}
					}
				}

				if(forbidden) {
					IloNumExpr expr = model.sum(vars[i], vars[j]);
					model.addLe(expr, 1);
				}
			}
		}
	}
	
	public void cleanUp() throws IloException {
		model.clearModel();
		model.end();
	}

	private void initObjective() throws IloException {
		Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
		IloNumExpr expr = model.constant(0);

		for(int i = 0; i < bix.getBranchMapping().size(); i++) {
			List<Integer> list = bix.getBranchMapping().get(i);
			double dist = 0;
			double dualValue = 0;
			for(Integer j: list) {
				if(j==currentVertex) {
					continue;
				}
				dist += GlobalParam.DISTANCE_FUNCTION.getSquaredDistanceBetweenPoints(instance, currentVertex, j) * (pair.first - pair.second + 1);
				dualValue += duals.getDualsPoint()[level][j];
			}
			IloNumExpr temp = model.prod(dist - dualValue, vars[i]);
			expr = model.sum(expr, temp);
		}

		objective = model.addMinimize(expr);
	}

	public void solve() {
		try {
			model.solve();
		} catch (IloException e) {
			e.printStackTrace();
		}
	}

	public Triple<List<Integer>, Double, Double> getClique() {
		List<Integer> result = new ArrayList<>();
		double obj = 0;
		try {
			for(int i = 0; i < bix.getBranchMapping().size(); i++) {
				if(model.getValue(vars[i])>0.5) {
					result.addAll(bix.getBranchMapping().get(i));
				}
			}
			obj = model.getObjValue();
		} catch (UnknownObjectException e) {
			e.printStackTrace();
		} catch (IloException e) {
			e.printStackTrace();
		}
		return new Triple<>(result, obj, obj);
	}

}
