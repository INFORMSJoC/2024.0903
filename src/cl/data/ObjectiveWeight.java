package cl.data;

public enum ObjectiveWeight {
	noWeight, // Default. Weight = 1
	noDoubleCounting, // weight = 1 / (k-h+1)
	proportionalToLevel // weight = (objective level 1 / objective level k)
}
