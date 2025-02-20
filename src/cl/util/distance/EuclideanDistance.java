package cl.util.distance;

public class EuclideanDistance {

	/**
	 * Officially the formula is euclideanDistance = sqrt(x^2)
	 * We always need the squared euclideanDistance = sqrt(x^2)^2
	 * @param coordinate1
	 * @param coordinate2
	 * @return
	 */
	public static double getSquaredEuclideanDistance(double[] coordinate1, double[] coordinate2) {
		double sum = 0;
		for(int i = 0; i < coordinate1.length; i++) {
			sum += Math.pow(coordinate1[i]-coordinate2[i], 2);
		}
		return (sum); 
	}
	
	public static double getRegularEuclideanDistance(double[] coordinate1, double[] coordinate2) {
		double sum = 0;
		for(int i = 0; i < coordinate1.length; i++) {
			sum += Math.pow(coordinate1[i]-coordinate2[i], 2);
		}
		return Math.sqrt(sum); 
	}
}
