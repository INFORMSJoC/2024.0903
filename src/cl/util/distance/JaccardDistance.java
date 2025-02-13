package cl.util.distance;

public class JaccardDistance {
	/**
	 * Jaccard similarity is intersection / union
	 * Jaccard distance is 1 - jac sim
	 * 
	 * @param coordinate1
	 * @param coordinate2
	 * @return
	 */
	public static double getJaccardDistance(double[] coordinate1, double[] coordinate2) {
		double intersection = 0;
		double union = 0;
		
		for(int i = 0; i < coordinate1.length; i++) {
			if(coordinate1[i]==1 && coordinate2[i]==1) {
				intersection += 1;
			}
			if(coordinate1[i]==1 || coordinate2[i]==1) {
				union += 1;
			}
		}
		double jaccardSim = intersection / (union+1);
		double jaccardDistance = 1 - jaccardSim;
		return jaccardDistance;
	}
}
