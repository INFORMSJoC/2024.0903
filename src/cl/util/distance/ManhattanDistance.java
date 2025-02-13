package cl.util.distance;

public class ManhattanDistance {

	/**
	 * Officially the formula is |x|
	 * 
	 * @param coordinate1
	 * @param coordinate2
	 * @return
	 */
	public static double getManhattanDistance(double[] coordinate1, double[] coordinate2) {
		double sum = 0;
		for(int i = 0; i < coordinate1.length; i++) {
			sum += Math.abs(coordinate1[i]-coordinate2[i]);
		}
		return (sum); 
	}
}
