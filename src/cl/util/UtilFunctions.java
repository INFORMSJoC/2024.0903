package cl.util;

public class UtilFunctions {
	public static double roundDecimals(double val, int decimals) {
		double factor = Math.pow(10, decimals);
		return (double) Math.round(val * factor) / factor;
	}
}
