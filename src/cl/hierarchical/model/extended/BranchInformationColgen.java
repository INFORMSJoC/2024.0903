package cl.hierarchical.model.extended;

import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import cl.data.Instance;
import cl.util.Pair;

public class BranchInformationColgen implements BranchInformation {
	private Instance instance;
	private Map<Integer, Set<Pair<Integer, Integer>>> enforced;
	private Map<Integer, Set<Pair<Integer, Integer>>> forbidden;

	public BranchInformationColgen(Instance instance) {
		this.instance = instance;
		this.enforced = new LinkedHashMap<>();
		this.forbidden = new LinkedHashMap<>();
		for(int i = 0; i < instance.getStartEndLevelList().size(); i++) {
			enforced.put(i, new LinkedHashSet<>());
			forbidden.put(i, new LinkedHashSet<>());
		}
	}

	public void enforce(Pair<Integer, Pair<Integer, Integer>> pair) {
		for(int i = 0; i < instance.getStartEndLevelList().size(); i++) {
			Pair<Integer, Integer> level = instance.getStartEndLevelList().get(i);
			if(level.second<= pair.first && level.first >= pair.first) {
				enforced.get(i).add(pair.second);
			}
		}
	}

	public void forbid(Pair<Integer, Pair<Integer, Integer>> pair) {
		for(int i = 0; i < instance.getStartEndLevelList().size(); i++) {
			Pair<Integer, Integer> level = instance.getStartEndLevelList().get(i);
			if(level.second<= pair.first && level.first >= pair.first) {
				forbidden.get(i).add(pair.second);
			}
		}
	}

	public void clear(Pair<Integer, Pair<Integer, Integer>> pair) {
		for(int i = 0; i < instance.getStartEndLevelList().size(); i++) {
			Pair<Integer, Integer> level = instance.getStartEndLevelList().get(i);
			if(level.second<= pair.first && level.first >= pair.first) {
				enforced.get(i).remove(pair.second);
				forbidden.get(i).remove(pair.second);
			}
		}
	}

	public Set<Pair<Integer, Integer>> getEnforced(int level) {
		return enforced.get(level);
	}

	public Set<Pair<Integer, Integer>> getForbidden(int level) {
		return forbidden.get(level);
	}
}
