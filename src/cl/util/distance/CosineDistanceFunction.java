package cl.util.distance;

import cl.data.Instance;
import cl.data.Medoid;

public class CosineDistanceFunction implements DistanceFunction {
	public double getDistanceOfMedoid(Instance instance, Medoid medoid)
	{
		double distance = 0;
		for(Integer i: medoid.getCluster()) {
			distance += instance.getSquaredDistanceMatrix()[i][medoid.getMedoidId()]; 
		}
		return distance;
	}
	
	public double getSquaredDistanceBetweenPoints(Instance instance, int coord1, int coord2) {
		return instance.getSquaredDistanceMatrix()[coord1][coord2];
	}
	
	public double getRegularDistanceBetweenPoints(Instance instance, int coord1, int coord2) {
		return 0;
	}
	
	public double getSquaredDistance(double[] coord1, double[] coord2) {
		return CosineDistance.getCosineDistance(coord1, coord2);
	}
	
	public double getRegularDistance(double[] coord1, double[] coord2) {
		return 0;
	}

	@Override
	public String toString() {
		return "CosineDistanceFunction";
	}
}
