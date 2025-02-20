package cl.hierarchical.model.pdcg;

import java.util.List;

import cl.util.Pair;

public class PDCGM_matrixAAT {
	private List<Pair<Integer, Double>> values;
	
	public PDCGM_matrixAAT(List<Pair<Integer, Double>> values) {
		this.values = values;
	}

	public List<Pair<Integer, Double>> getValues() {
		return values;
	}
	
	
}
