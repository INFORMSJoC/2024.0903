package cl.data.type;

public enum PricingProblemStrategy {
	None,
	StartFirstEndSecond, // Bottom up. Solve the PP's in this order: 44 -> 43 -> 42 -> 41 -> 33 -> 32 -> 31. (no restarts)
	StartFirstEndSecondWithRestart, // Bottom up. Solve the PP's in this order: 44 -> 43 -> 42 -> 41. Then restart and do 33 -> 32 -> 31.
	EndFirstStartSecond, // Top down approach
	EndFirstStartSecondWithRestart // Top down approach
}
