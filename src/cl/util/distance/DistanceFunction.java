package cl.util.distance;

import cl.data.Instance;
import cl.data.Medoid;

public interface DistanceFunction {
	double getDistanceOfMedoid(Instance instance, Medoid medoid);
	double getSquaredDistanceBetweenPoints(Instance instance, int coord1, int coord2);
	double getRegularDistanceBetweenPoints(Instance instance, int coord1, int coord2);
	double getSquaredDistance(double[] coord1, double[] coord2);
	double getRegularDistance(double[] coord1, double[] coord2);
}
