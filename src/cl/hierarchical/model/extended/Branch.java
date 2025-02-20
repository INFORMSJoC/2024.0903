package cl.hierarchical.model.extended;

import cl.util.Pair;

public class Branch {
	private Pair<Integer, Pair<Integer, Integer>> branching;
	private boolean enforce;
	
	public Branch(Pair<Integer, Pair<Integer, Integer>> branching, boolean enforce) {
		this.branching = branching;
		this.enforce = enforce;
	}

	public Pair<Integer, Pair<Integer, Integer>> getBranching() {
		return branching;
	}

	public boolean isEnforce() {
		return enforce;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((branching == null) ? 0 : branching.hashCode());
		result = prime * result + (enforce ? 1231 : 1237);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Branch other = (Branch) obj;
		if (branching == null) {
			if (other.branching != null)
				return false;
		} else if (!branching.equals(other.branching))
			return false;
		if (enforce != other.enforce)
			return false;
		return true;
	}
	
	
}
