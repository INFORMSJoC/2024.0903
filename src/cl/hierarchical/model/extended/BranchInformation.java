package cl.hierarchical.model.extended;

import java.util.Set;

import cl.util.Pair;

public interface BranchInformation {
	public Set<Pair<Integer, Integer>> getEnforced(int level);
	public Set<Pair<Integer, Integer>> getForbidden(int level);
	public void enforce(Pair<Integer, Pair<Integer, Integer>> pair);
	public void forbid(Pair<Integer, Pair<Integer, Integer>> pair);
	public void clear(Pair<Integer, Pair<Integer, Integer>> pair);
}
