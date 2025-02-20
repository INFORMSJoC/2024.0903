package cl.hierarchical.model.extended;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import cl.data.Cluster;
import cl.util.Pair;

public class CuttingPlanes<E> {

	private Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<E, Double>> cMapOriginal;
	private Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, E>>>> cMapConstraint;
	private List<List<E>> cuts;
	private double maxNumCuts = Double.POSITIVE_INFINITY;

	public CuttingPlanes(Map<Pair<Pair<Integer, Integer>, Cluster>, Pair<E, Double>> cMapOriginal,
			Map<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, E>>>> cMapConstraint) {
		this.cMapOriginal = cMapOriginal;
		this.cMapConstraint = cMapConstraint;
		this.cuts = new ArrayList<>();
	}

	public void setMaxNumCuts(double maxNumCuts) {
		this.maxNumCuts = maxNumCuts;
	}

	private void addConstraint(Map<Pair<Integer, Integer>, Map<Cluster, E>> constraintMap, Entry<Pair<Pair<Integer, Integer>, Cluster>, Pair<E, Double>> entry) {
		List<E> cut = new ArrayList<>();
		for(Entry<Pair<Integer, Integer>, Map<Cluster, E>> entry2: constraintMap.entrySet()) {
			for(Entry<Cluster, E> entry3: entry2.getValue().entrySet()) {
				cut.add(entry3.getValue());
			}
		}
		cut.add(entry.getValue().first);
		cuts.add(cut);
	}

	public List<List<E>> findValidInequalities2() {
		for(Entry<Integer, Map<Integer, Map<Pair<Integer, Integer>, Map<Cluster, E>>>> entry: cMapConstraint.entrySet()) {
			int singleLevel = entry.getKey();
			for(Entry<Integer, Map<Pair<Integer, Integer>, Map<Cluster, E>>> entry2: entry.getValue().entrySet()) {
				// Look for a violation that contains all variables from this constraint
				Map<Pair<Integer, Integer>, Map<Cluster, E>> constraintMap = entry2.getValue();

				for(Entry<Pair<Pair<Integer, Integer>, Cluster>, Pair<E, Double>> entry3: cMapOriginal.entrySet()) {
					// Does this cluster violate all clusters in constraintMap?
					Cluster cluster2 = entry3.getKey().second;
					Pair<Integer, Integer> level2 = entry3.getKey().first;
					
					// Check whether cluster2 is already in constraintMap
					if(constraintMap.containsKey(level2) && constraintMap.get(level2).containsKey(cluster2)) {
						continue;
					}

					boolean allHaveViolation = true;
					for(Entry<Pair<Integer, Integer>, Map<Cluster, E>> element: constraintMap.entrySet()) {
						Pair<Integer, Integer> level = element.getKey();

						for(Entry<Cluster, E> entry4: element.getValue().entrySet()) {
							Cluster cluster = entry4.getKey();

							boolean violation = false;
							if(cluster.equals(cluster2) && (singleLevel>level2.first || level2.second>singleLevel )){
								violation = true;
							}
							else {
								if(level.first>=level2.second) {
									for(int i: cluster.getPointId()) {
										if(Cluster.checkValidInequalityType1(cluster, cluster2, i)) {
											violation = true;
											break;
										}
									}
								}
							}
							if(!violation) {
								allHaveViolation = false;
							}	
						}
					}
					if(allHaveViolation) {
						addConstraint(constraintMap, entry3);
						
						if(cuts.size()>=maxNumCuts) {
							return cuts;
						}
					}
				}
			}
		}
		return cuts;
	}	
}