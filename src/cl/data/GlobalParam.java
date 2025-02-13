package cl.data;

import java.util.EnumSet;
import java.util.Random;

import cl.data.type.DistanceType;
import cl.data.type.PricingProblemStrategy;
import cl.data.type.ValidInequalityType;
import cl.util.distance.CosineDistanceFunction;
import cl.util.distance.DistanceFunction;
import cl.util.distance.JaccardDistanceFunction;
import cl.util.distance.ManhattanDistanceFunction;
import cl.util.distance.SquaredEuclideanDistanceFunction;

public class GlobalParam {
	// Random generator.
	public static Random RANDOM = new Random(32903920);
	
	// Time limit
	public static final boolean USE_GLOBAL_TIME_LIMIT = false;
	public static final long GLOBAL_TIME_LIMIT = 60*60*1000; // 1 hour
	// 60*60*1000; // 1 hour
	// 5*60*1000; // 5 minutes		
	
	// Write solution file when finding a better integer solution
	public static boolean WRITE_BETTER_INTEGER_SOLUTION = false;
	public static String WRITE_BETTER_INTEGER_SOLUTION_FILE_NAME;
	
	// Debug
	public static final boolean DEBUG = true;
	public static final boolean PRINT_DETAILS = true;
	
	// Constants
	public static final double EPSILON1 = 1e-1;
	public static final double EPSILON2 = 1e-2;
	public static final double EPSILON3 = 1e-3;
	public static final double EPSILON4 = 1e-4;
	public static final double EPSILON5 = 1e-5;
	public static final double EPSILON7 = 1e-7;
	public static final double EPSILON10 = 1e-10;
	public static final double BIGM8 = 1e8;
	public static final double BIGM6 = 1e6;
	public static final double BIGM3 = 1e3;
	
	// Distance function
	public static final DistanceType DISTANCE_TYPE = DistanceType.Centroid; // We either use Centroid with SquaredEuclideanDistance or Medoid with DistanceFunction
	public static final DistanceFunction DISTANCE_FUNCTION = new SquaredEuclideanDistanceFunction();
	
	// Weight in objective function
	public static ObjectiveWeight OBJECTIVE_WEIGHT = ObjectiveWeight.noWeight;
	
	// Column generation
	public static final boolean USE_VARIABLES_V = true; // Switch between 2 formulations (one with variablesV and one without)
	
	// Column management
	public static final boolean USE_COLUMN_MANAGEMENT = false;
	public static final double COLUMN_MANAGEMENT_THRESHOLD = 10; // Reduced cost
	public static final int COLUMN_MANAGEMENT_NUM_ITER = 25; // Number of inactive iterations
	
	// Valid inequalities
	public static final boolean APPLY_CUTTING_PLANES_AFTER_GENERATING_COLUMNS = false;
	public static final EnumSet<ValidInequalityType> VALID_INEQUALITIES = EnumSet.of(ValidInequalityType.Type1, ValidInequalityType.Type2);
	public static final int CUTTING_PLANES_NUM_CUTS = 2;
	public static final int BRANCH_PRICE_CUT_ENUMERATE_NUM_CUTTING_PLANES = 10;
	
	// Branch and price reduced cost
	public static final double NEGATIVE_RC_THRESHOLD = -EPSILON3;
	
	// Cluster enumeration
	public static double ENUMERATE_GRAPH_THRESHOLD; // Create a graph with a specific range of values (e.g. [-10, 0] or [0, 5])
	public static double ENUMERATE_RC_COST_THRESHOLD; // Accept clusters that are below this RC cost value
	public static final boolean APPLY_CUTTING_PLANES_AFTER_ENUMERATE = false;
	public static final boolean ENUMERATE_ADD_TO_EACH_LEVEL = true;
	
	// Cluster enumeration, how many columns to add
	public static final boolean ENUMERATE_ADD_ALL_COLUMNS = true; // True, if we add all clusters with negative RC
	public static final int ENUMERATE_PP_NUM_COLUMNS_ADD_LIMIT = 50; // If number of columns is below this number, we add all clusters
	public static final int ENUMERATE_PP_NUM_COLUMNS = 20; // Otherwise we add this number of clusters
	public static final int ENUMERATE_PP_NUM_COLUMNS_MAXIMUM_IN_GRAPH = Integer.MAX_VALUE; // Dont keep track of too many clusters in the graph
	
	// Branch and price stopping condition
	public static final double OPTIMALITY_GAP = 0.01;
	
	// Pricing problem 
	public static boolean PP_STOP_EARLY = false;
	public static final boolean PP_STOP_WHEN_REACHING_POTENTIAL_COLUMNS = true;
	public static final int PP_NUM_POTENTIAL_COLUMNS = 20;
	public static final int PP_NUM_COLUMNS = 2;
	public static final boolean PP_USE_EXCLUSION_CHECK = false; 
	public static final boolean PP_RESTRICTED_GRAPH = false; // Restrict the graph based on pricing problem strategies
	
	// QP solver
	public static boolean ONLY_USE_VNS = true;
	public static final boolean USE_ROOF_DUAL_BOUND = true;
	public static final boolean USE_PERSISTENCY_RESULT = true;
	public static final int VNS_MAX_NUM_NEIGHBOURHOODS = 50;
	
	// Heuristics
	public static final int HEURISTIC_RESTARTS = 50;
	
	// Counters
	public static int COUNTER_VNS;
	public static int COUNTER_QP;
	public static int COUNTER_BRANCH_AND_BOUND;
	public static int COUNTER_KNAPSACK;
	public static int COUNTER_CUTTING_PLANES;
	public static int COUNTER_ENUMERATE;
	
	// ACCP
	public static final int ACCP_LINE_SEARCH_ITERATION_LIMIT = 50;
	public static final double ACCP_ALPHA_BISECTION_U0 = 1e-5;
	public static final double ACCP_UMAX = 1e5;
	public static final double ACCP_ALPHA_BISECTION_GAMMA = 10;
	public static final double ACCP_ALPHA_BISECTION_DELTA = 1e-4;
}
