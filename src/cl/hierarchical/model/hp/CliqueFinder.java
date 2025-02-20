package cl.hierarchical.model.hp;

import java.util.List;

import cl.util.Triple;

public interface CliqueFinder {
	void solve();
	Triple<List<Integer>, Double, Double> getClique(); // First double is the negative RC without kappa. Second double is negative RC with kappa.
}
