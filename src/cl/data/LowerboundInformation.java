package cl.data;

import java.util.Map;

public class LowerboundInformation {
	Map<Integer, Double> lowerboundMap, cumulativeLowerboundMap;

	public LowerboundInformation(Map<Integer, Double> lowerboundMap, Map<Integer, Double> cumulativeLowerboundMap) {
		this.lowerboundMap = lowerboundMap;
		this.cumulativeLowerboundMap = cumulativeLowerboundMap;
	}

	public Map<Integer, Double> getLowerboundMap() {
		return lowerboundMap;
	}

	public Map<Integer, Double> getCumulativeLowerboundMap() {
		return cumulativeLowerboundMap;
	}	
}