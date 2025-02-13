package cl.util.distance;

public class CosineDistance {
	/**
	 * Cosine similarity is dot(a, b) / norm(a) * norm(b)
	 * Cosine distance is 1 - cos sim
	 * 
	 * @param coordinate1
	 * @param coordinate2
	 * @return
	 */
	public static double getCosineDistance(double[] coordinate1, double[] coordinate2) {
		double norm1 = 0;
		double norm2 = 0;
		double innerProduct = 0;
		for(int i = 0; i < coordinate1.length; i++) {
			norm1 += Math.pow(coordinate1[i], 2);
			norm2 += Math.pow(coordinate2[i], 2);
			innerProduct += coordinate1[i] * coordinate2[i];
		}
		norm1 = Math.sqrt(norm1);
		norm2 = Math.sqrt(norm2);
		double cosSim = innerProduct / (norm1 * norm2);
		double cosDistance = 1 - cosSim;
		return cosDistance;
	}
}
