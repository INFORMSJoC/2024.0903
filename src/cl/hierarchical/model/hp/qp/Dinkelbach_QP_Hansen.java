package cl.hierarchical.model.hp.qp;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cl.data.GlobalParam;
import cl.data.Instance;
import cl.data.ObjectiveWeightFunction;
import cl.hierarchical.model.extended.BranchInformationExtra;
import cl.hierarchical.model.extended.DualInformation;
import cl.util.Pair;
import cl.util.Triple;

/**
 * QUBO solver based on the description in:
 * - A simple enumerative algorithm for unconstrained 0-1 quadratic programming by Pierre Hansen, Brigitte Jaumard, and Christophe Meyer
 */
public class Dinkelbach_QP_Hansen implements QP_Solver {
	private static Logger log = LoggerFactory.getLogger(QP_Solver.class);

	private Instance instance;
	private int level;
	private DualInformation duals;
	private List<Integer> Gi;
	private BranchInformationExtra bix;

	private double prevObj;

	// Original matrix
	private double[][] matrixA;

	// Posiform
	private double[][] matrixW;
	private double[] vectorC;
	private double psiZero; // Current lowerbound

	// Main algorithm
	private int[] permutationVector;
	private boolean[] bestX;

	private double coeffMatrix[][]; // Store coefficients for more efficient computation

	// Fix variables
	private List<Integer> persistencyFix = new ArrayList<>();

	public Dinkelbach_QP_Hansen(Instance instance, int level, DualInformation duals, List<Integer> Gi, BranchInformationExtra bix, double prevObj)
	{
		this.instance = instance;
		this.level = level;
		this.duals = duals;
		this.Gi = Gi;
		this.bix = bix;

		this.prevObj = prevObj;

		initCoeffMatrix();
	}

	public void solve() {
		int n = Gi.size();

		// Initialise best known minimum and minimiser
		bestX = new boolean[n]; // Algorithm always starts with false everywhere
		double bestObj = calcObjectiveQP(bestX).first;

		// Current solution and current permutation vector
		boolean[] curX = bestX.clone();
		permutationVector = new int[n];
		for(int i = 0; i < n; i++) {
			permutationVector[i] = i;
		}

		// Initialise bounds on partial derivative
		double[] partialDerivativeLB = new double[n];
		double[] partialDerivativeUB = new double[n];

		// Create a stack
		Stack<Triple<Integer, double[], Double>> stack = new Stack<>();

		int level = -1;
		boolean update = false;

		// Create posiform. 
		// Note that the quadratic part stays fixed
		// Only the linear part will be modified
		createPosiform(); 

		if(GlobalParam.USE_ROOF_DUAL_BOUND) {
			// Calculate roofDualBound
			long time = System.currentTimeMillis(); 
			roofDualBound();
			if(GlobalParam.PRINT_DETAILS) {
				log.info("Finished roofDualBound: {} in seconds", (System.currentTimeMillis() - time)/1000d);
			}

			if(GlobalParam.USE_PERSISTENCY_RESULT) {
				// Fix entries to 0, due to roofDualBound (persistency result. Section 3.1 Hansen 2000)
				fixEntriesBasedOnPersistency();
			}
		}

		boolean cont = true;
		while(cont) {
			if(psiZero>=bestObj || level==n-1) {
				if(psiZero<bestObj) {
					bestObj = psiZero; 
					bestX = curX.clone();
				}

				if(stack.isEmpty()) {
					// If stack is empty, we are done
					break;
				}
				Triple<Integer, double[], Double> triple = stack.pop();
				level = triple.first;
				vectorC = triple.second;
				psiZero = triple.third;

				curX[permutationVector[level]] = !curX[permutationVector[level]];
				update = false;
			}
			else {
				int idxFix = -1;
				int idxFix0 = -1;
				if(!persistencyFix.isEmpty()) {
					idxFix0 = persistencyFix.get(0);
					idxFix = persistencyFix.get(0);
					persistencyFix.remove(0);
					
				}
				if(idxFix0==-1) {
					int idx = -1; 
					double max1 = Double.NEGATIVE_INFINITY;
					for(int i = level+1; i < n; i++) {
						double bound = psiZero+Math.max(vectorC[2*permutationVector[i]], vectorC[2*permutationVector[i]+1]);
						if(bound>=bestObj) {
							if(bound>max1) {
								idx = i;
								idxFix = i;
								max1 = bound; 
							}
						}
					}
					if(idx!=-1) {
						if(vectorC[2*permutationVector[idx]] > vectorC[2*permutationVector[idx]+1]) {
							curX[permutationVector[idx]] = false;
						}
						else {
							curX[permutationVector[idx]] = true;
						}
						update = false;
					}
					else {
						// Recompute bounds on partial derivates of free variables
						if(update) { //If only one variable has changed, use update procedure of Pardalos and Rodgers
							if(curX[permutationVector[level]]) {
								for(int i = level+1; i < n; i++) {
									partialDerivativeLB[permutationVector[i]] += Math.max(matrixA[permutationVector[i]][permutationVector[level]], 0);
									partialDerivativeUB[permutationVector[i]] += Math.min(matrixA[permutationVector[i]][permutationVector[level]], 0); 
								}
							}
							else {
								for(int i = level+1; i < n; i++) {
									partialDerivativeLB[permutationVector[i]] -= Math.min(matrixA[permutationVector[i]][permutationVector[level]], 0);
									partialDerivativeUB[permutationVector[i]] -= Math.max(matrixA[permutationVector[i]][permutationVector[level]], 0); 
								}
							}
						}
						else { //If more than one variable has changed, recompute from scratch
							for(int i = level+1; i < n; i++) {
								int pi = permutationVector[i];
								partialDerivativeLB[pi] = matrixA[pi][pi];
								partialDerivativeUB[pi] = matrixA[pi][pi]; 
								for(int j = 1; j < level ; j++) {
									int pj = permutationVector[j];
									if(curX[pj]) {
										partialDerivativeLB[pi] += matrixA[pi][pj];
										partialDerivativeUB[pi] += matrixA[pi][pj];
									}
								}
								for(int j = level+1; j < n; j++) {
									if(i==j) {
										continue;
									}
									int pj = permutationVector[j];
									partialDerivativeLB[pi] += Math.min(matrixA[pi][pj], 0);
									partialDerivativeUB[pi] += Math.max(matrixA[pi][pj], 0);
								}
							}
						}
						update = true;

						// Check if a variable can be fixed by using partial derivatives
						int idx2 = -1;
						double largestAbs = Double.NEGATIVE_INFINITY;
						for(int i = level+1; i < n; i++) {
							if(partialDerivativeLB[permutationVector[i]]>=0 || partialDerivativeUB[permutationVector[i]]<=0) {
								if(partialDerivativeLB[permutationVector[i]]>=largestAbs)
								{
									idx2 = i;
									idxFix = i;
									largestAbs = partialDerivativeLB[permutationVector[i]];
								}
								else if(-1*partialDerivativeUB[permutationVector[i]]<=largestAbs)
								{
									idx2 = i;
									idxFix = i;
									largestAbs = -1*partialDerivativeUB[permutationVector[i]];
								}
							}
						}
						if(idx2!=-1) {
							if(partialDerivativeLB[permutationVector[idx2]]>=0) {
								curX[permutationVector[idx2]] = false;
							}
							else {
								curX[permutationVector[idx2]] = true;
							}
						}
						else {
							// Else, we have to branch
							double kmax = Double.NEGATIVE_INFINITY;
							for(int k = level+1; k < n; k++) {
								double candidate = Math.min(-partialDerivativeLB[permutationVector[k]], partialDerivativeUB[permutationVector[k]]);
								if(candidate>kmax) {
									kmax = candidate;
								}
							}

							int idx3 = -1;
							double maxkCriterion = Double.NEGATIVE_INFINITY;
							for(int j = level+1; j < n; j++) {
								double jmin = Math.min(-partialDerivativeLB[permutationVector[j]], partialDerivativeUB[permutationVector[j]]);
								if(jmin == kmax) {
									if(Math.max(-partialDerivativeLB[permutationVector[j]], partialDerivativeUB[permutationVector[j]])>maxkCriterion) {
										idx3 = j;
										idxFix = j; 
										maxkCriterion = Math.max(-partialDerivativeLB[permutationVector[j]], partialDerivativeUB[permutationVector[j]]);
									}
								}
							}

							// Alternative permutation vector
							int[] altPermutationVector = permutationVector.clone();
							int temp = altPermutationVector[idx3];
							altPermutationVector[idx3] = altPermutationVector[level+1];
							altPermutationVector[level+1] = temp;

							// Update the linear part of the posiform and the lowerbound psiZero
							curX[permutationVector[idx3]] = false;
							Pair<double[], Double> pairZero = updateStep(curX, vectorC, altPermutationVector, psiZero, level+1, n);
							double psiZeroTempZero = pairZero.second;

							curX[permutationVector[idx3]] = true;
							Pair<double[], Double> pairOne = updateStep(curX, vectorC, altPermutationVector, psiZero, level+1, n);
							double psiZeroTempOne = pairOne.second;

							if(psiZeroTempZero<psiZeroTempOne) {
								curX[permutationVector[idx3]] = false;
							}
							else {
								curX[permutationVector[idx3]] = true;
							}

							// Stack the alternative branch, if the lowerbound is below the best known objective
							boolean[] alternativeX = curX.clone();
							alternativeX[permutationVector[idx3]] = !alternativeX[permutationVector[idx3]];
							// Update the linear part of the posiform and the lowerbound psiZero
							double[] altVectorC = vectorC.clone();
							double altPsiZero = psiZero;
							Pair<double[], Double> pair = updateStep(alternativeX, altVectorC, altPermutationVector, altPsiZero, level+1, n);
							altVectorC = pair.first;
							altPsiZero = pair.second;
							if(altPsiZero<bestObj) {
								stack.push(new Triple<>(level+1, altVectorC, altPsiZero));
							}
						}
					}
				}
				level++;

				// Update permutation vector 
				int temp = permutationVector[idxFix];
				permutationVector[idxFix] = permutationVector[level];
				permutationVector[level] = temp;

				// Update the linear part of the posiform and the lowerbound psiZero
				Pair<double[], Double> pair = updateStep(curX, vectorC, permutationVector, psiZero, level, n);
				vectorC = pair.first;
				psiZero = pair.second;
			}
		}
	}

	private void roofDualBound() {
		Dinkelbach_QP_Hansen_RoofDualBound roofDualBoundSolver = new Dinkelbach_QP_Hansen_RoofDualBound(Gi, matrixW, vectorC, psiZero);

		// Find the roof dual bound
		roofDualBoundSolver.solve(); 
	
		// Update the posiform
		matrixW = roofDualBoundSolver.getMatrixW();
		vectorC = roofDualBoundSolver.getVectorC();
		psiZero = roofDualBoundSolver.getPsiZero();
	}

	private void fixEntriesBasedOnPersistency() {
		int n = Gi.size();
		boolean performFix = false;
		for(int i = 0; i < 2*n; i++) {
			if(vectorC[i]>0) {
				performFix = true;
				break;
			}
		}

		if(performFix) {
			for(int i = 0; i < 2*n; i++) {
				if(vectorC[i]>0) {
					if(i%2==0) {
						persistencyFix.add(i/2);
					}
					else {
						persistencyFix.add((i-1)/2);
					}
				}
			}
		}
	}

	/**
	 * Return vectorC and psiZero
	 * @param curX
	 * @param vectorC
	 * @param psiZero
	 * @param i
	 * @param level
	 * @param n
	 * @return
	 */
	private Pair<double[], Double> updateStep(boolean[] curX, double[] vectorC, int[] permutationVector, double psiZero, int level, int n) {
		double[] newVectorC = vectorC.clone();
		int pi = permutationVector[level];
		if(curX[pi]) {
			for(int j = level+1; j < n; j++) {
				int pj = permutationVector[j];
				int a = 2*pi;
				int b = 2*pj;
				if(pi>pj) {
					a = 2*pj;
					b = 2*pi;
				}
				newVectorC[2*pj] += matrixW[a][b];
			}
			for(int j = level+1; j < n; j++) {
				int pj = permutationVector[j];
				int a = 2*pi;
				int b = 2*pj+1;
				if(pi>pj) {
					a = 2*pj+1;
					b = 2*pi;
				}
				newVectorC[2*pj+1] += matrixW[a][b];
			}
		}
		else {
			for(int j = level+1; j < n; j++) {
				int pj = permutationVector[j];
				int a = 2*pi+1;
				int b = 2*pj;
				if(pi>pj) {
					a = 2*pj;
					b = 2*pi+1;
				}
				newVectorC[2*pj] += matrixW[a][b];
			}
			for(int j = level+1; j < n; j++) {
				int pj = permutationVector[j];
				int a = 2*pi+1;
				int b = 2*pj+1;
				if(pi>pj) {
					a = 2*pj+1;
					b = 2*pi+1;
				}
				newVectorC[2*pj+1] += matrixW[a][b];
			}
		}

		double newPsiZero = psiZero;
		for(int j = level+1; j<n; j++) {
			int pj = permutationVector[j];
			double m = Math.min(newVectorC[2*pj], newVectorC[2*pj+1]);
			newPsiZero += m;
			newVectorC[2*pj] -= m;
			newVectorC[2*pj+1] -= m;
		}

		// Not in Hansen (2000)
		// But i believe we should update psiZero using the current index as well, since we fix this variable
		if(curX[pi]) {
			newPsiZero += newVectorC[2*pi];
		}
		else {
			newPsiZero += newVectorC[2*pi+1];
		}

		return new Pair<>(newVectorC, newPsiZero);
	}

	/**
	 * indices 0, 2, 4, 6... correspond to x
	 * indices 1, 3, 5, 7... correspond to x bar
	 */
	private void createPosiform()
	{
		int n = Gi.size();
		matrixA = new double[n][n];
		matrixW = new double[2*n][2*n];
		vectorC = new double[2*n];

		psiZero = 0;

		double dist = 0;
		if(!bix.isBranch()) {
			Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
			dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
		}

		for(int i = 0; i < n; i++) {
			for(int j = i+1; j < n; j++) {
				double coeff = 0;
				if(!bix.isBranch()) {	
					coeff = dist * instance.getSquaredDistanceMatrix()[Gi.get(i)][Gi.get(j)] - duals.getDualsPoint()[level][Gi.get(i)] - duals.getDualsPoint()[level][Gi.get(j)];
					if(Math.sqrt(dist) * instance.getRegularDistanceMatrix()[Gi.get(i)][Gi.get(j)] - duals.getDualsSqrtPoint()[level][Gi.get(i)] - duals.getDualsSqrtPoint()[level][Gi.get(j)]>0) {
						coeff = GlobalParam.BIGM8;
					}
				}
				else {
					coeff = bix.getReducedInstance().getSquaredDistanceMatrix()[Gi.get(i)][Gi.get(j)]
							- bix.getBranchMappingWeights().get(Gi.get(j)) * bix.getReducedDuals().getDualsPointSingle()[Gi.get(i)]
									- bix.getBranchMappingWeights().get(Gi.get(i)) * bix.getReducedDuals().getDualsPointSingle()[Gi.get(j)];
					if(bix.getReducedInstance().getRegularDistanceMatrix()[Gi.get(i)][Gi.get(j)]
							- bix.getBranchMappingWeights().get(Gi.get(j)) * bix.getReducedDuals().getDualsSqrtPointSingle()[Gi.get(i)]
									- bix.getBranchMappingWeights().get(Gi.get(i)) * bix.getReducedDuals().getDualsSqrtPointSingle()[Gi.get(j)]>0) {
						coeff = GlobalParam.BIGM8;
					}
				}

				if(coeff>0) {
					matrixW[2*i][2*j] = coeff;
				}
				else {
					matrixW[2*i+1][2*j] = -coeff;
					vectorC[2*j] += coeff;
				}
				matrixA[i][j] = coeff;
				matrixA[j][i] = coeff;
			}
		}

		for(int i = 0; i < n; i++) {
			double coeff = 0;
			double coeffA = 0;
			if(!bix.isBranch()) {
				coeffA = - duals.getDualsPoint()[level][Gi.get(i)] - prevObj;
				coeff = vectorC[2*i] + coeffA;
			}
			else {
				double w = bix.getBranchMappingWeights().get(Gi.get(i));
				coeffA = -w * bix.getReducedDuals().getDualsPointSingle()[Gi.get(i)] - bix.getReducedInstance().getSquaredDistanceMatrix()[Gi.get(i)][Gi.get(i)] - w * prevObj;
				coeff = vectorC[2*i] + coeffA;
			}
			vectorC[2*i] = 0; // Overwrite
			if(coeff>0) {
				vectorC[2*i] = coeff;
			}
			else {
				vectorC[2*i+1] = -coeff;
				psiZero += coeff;
			}
			matrixA[i][i] = coeffA;
		}
	}

	@Override
	public void setStartSolution(boolean[] startSolution) {
		throw new IllegalArgumentException("Not implemented");
	}

	@Override
	public boolean[] getSolution() {
		throw new IllegalArgumentException("Not implemented");
	}

	private void initCoeffMatrix()
	{
		double dist = 0;
		if(!bix.isBranch()) {
			Pair<Integer, Integer> pair = instance.getStartEndLevelList().get(level);
			dist = ObjectiveWeightFunction.getWeight(pair.first, pair.second);
		}

		coeffMatrix = new double[Gi.size()][Gi.size()];
		for(int i = 0; i < Gi.size(); i++) {
			for(int j = i+1; j < Gi.size(); j++) {
				if(!bix.isBranch()) {
					coeffMatrix[i][j] = dist * instance.getSquaredDistanceMatrix()[Gi.get(i)][Gi.get(j)] - duals.getDualsPoint()[level][Gi.get(i)] - duals.getDualsPoint()[level][Gi.get(j)];
				}
				else {
					coeffMatrix[i][j] = bix.getReducedInstance().getSquaredDistanceMatrix()[Gi.get(i)][Gi.get(j)]
							- bix.getBranchMappingWeights().get(Gi.get(j)) * bix.getReducedDuals().getDualsPointSingle()[Gi.get(i)]
									- bix.getBranchMappingWeights().get(Gi.get(i)) * bix.getReducedDuals().getDualsPointSingle()[Gi.get(j)];
				}
			}
		}
	}

	private Pair<Double, Double> calcObjectiveQP(boolean[] x)
	{
		double numerator = 0;
		for(int i = 0; i < Gi.size(); i++) {
			for(int j = i+1; j < Gi.size(); j++) {
				if(x[i] && x[j])
				{
					numerator += coeffMatrix[i][j];
				}
			}
		}
		double obj_QP = 0;
		double denominator = 0;
		for(int i = 0; i < Gi.size(); i++) {
			if(x[i])
			{
				if(!bix.isBranch()) {
					numerator -= duals.getDualsPoint()[level][Gi.get(i)];
					obj_QP -= prevObj;
					denominator ++;
				}
				else {
					double w = bix.getBranchMappingWeights().get(Gi.get(i));
					numerator -= w * bix.getReducedDuals().getDualsPointSingle()[Gi.get(i)] - bix.getReducedInstance().getSquaredDistanceMatrix()[Gi.get(i)][Gi.get(i)];
					obj_QP -= w *prevObj;
					denominator+= w;
				}
			}
		}

		double obj_HP = numerator / denominator;
		obj_QP += numerator;
		return new Pair<>(obj_QP, obj_HP);
	}

	public Pair<Double, Double> getObjValue() {
		return calcObjectiveQP(bestX);
	}

	public List<Integer> getClique() {
		List<Integer> result = new ArrayList<>();
		for(int i = 0; i < Gi.size(); i++) {
			if(bestX[i]) {
				result.add(Gi.get(i));
			}
		}
		return result;
	}

}
