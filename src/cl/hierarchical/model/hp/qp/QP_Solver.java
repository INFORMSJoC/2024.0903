package cl.hierarchical.model.hp.qp;

import java.util.List;

import cl.util.Pair;

public interface QP_Solver {
	void solve();
	void setStartSolution(boolean[] startSolution);
	boolean[] getSolution();
	Pair<Double, Double> getObjValue();
	List<Integer> getClique();
}
