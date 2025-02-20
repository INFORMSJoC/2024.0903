package cl.hierarchical.model.pdcg;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cl.data.GlobalParam;
import cl.util.MatrixOperations;
import cl.util.Pair;

/**
 * Primal dual interior point method, based on the description in the following sources:
 * - Primal-Dual Interior-Point Methods for Linear Programming based on Newton's Method by Robert M. Freund (March, 2004)
 * - https://en.wikipedia.org/wiki/Mehrotra_predictor%E2%80%93corrector_method
 */
public class PrimalDualInteriorPointMethod {
	private static Logger log = LoggerFactory.getLogger(PrimalDualInteriorPointMethod.class);
	
	private RealMatrix startMatrixA, startVectorB, startVectorC;
	private PrimalDualPoint startPrimalDualPoint;
	
	// Constants
	private final double epsilon;
	private final double r = 0.99;
	private final double alpha = 1/10d; // Not used
	private final double beta = 3/40d; // Not used
	
	private double feasibleEpsilon = GlobalParam.EPSILON2;
	
	// The solution has to be primal feasible
	private final boolean primalFeasible; 
	private boolean isPrimalFeasible;
	
	//Solution
	private PrimalDualPoint solutionPrimalDualPoint;
	private int previousIndex;
	private PDCGM_matrixAAT[][] matrixAAT_final = null;
	
	public PrimalDualInteriorPointMethod(RealMatrix startMatrixA, RealMatrix startVectorB, RealMatrix startVectorC, double epsilon, boolean allFeasible) {
		this.startMatrixA = startMatrixA;
		this.startVectorB = startVectorB;
		this.startVectorC = startVectorC;
		this.epsilon = epsilon;
		this.primalFeasible = allFeasible;
	}

	public void setStartPrimalDualPoint(PrimalDualPoint startPrimalDualPoint) {
		this.startPrimalDualPoint = startPrimalDualPoint;
	}
	
	private PrimalDualPoint initPrimalDualPoint() {
		if(startPrimalDualPoint!=null) {
			return startPrimalDualPoint;
		}
		
		RealMatrix startingPointX = new Array2DRowRealMatrix(startMatrixA.getColumnDimension(), 1);
		RealMatrix startingPointPi = new Array2DRowRealMatrix(startMatrixA.getRowDimension(), 1);
		RealMatrix startingPointS = new Array2DRowRealMatrix(startMatrixA.getColumnDimension(), 1);
		MatrixOperations.setRandomElementsMatrix(startingPointX, 0.1);
		MatrixOperations.setRandomElementsMatrix(startingPointPi, 10);
		MatrixOperations.setRandomElementsMatrix(startingPointS, 0.1);
		return new PrimalDualPoint(startingPointX, startingPointPi, startingPointS);
	}

	public PDCGM_matrixAAT[][] getMatrixAAT() {
		return matrixAAT_final;
	}
	
	public void setMatrixAAT(PDCGM_matrixAAT[][] matrixAAT_start, int previousIndex) {
		this.matrixAAT_final = matrixAAT_start;
		this.previousIndex = previousIndex;
	}

	public void setFeasibleEpsilon(double feasibleEpsilon) {
		this.feasibleEpsilon = feasibleEpsilon;
	}

	public void solve() {
		// Step 1
		PrimalDualPoint curPoint = initPrimalDualPoint();
		RealMatrix matrixA = startMatrixA;
		RealMatrix vectorB = startVectorB;
		RealMatrix vectorC = startVectorC;
		
		if(matrixAAT_final==null) {
			matrixAAT_final = precomputeAAT_fromScratch(matrixA);
		}
		else {
			precomputeAAT_extend(matrixAAT_final, matrixA, previousIndex);
		}
		
		double n = curPoint.getPointX().getRowDimension();
		
		isPrimalFeasible = true;
		int iter = 0;
		while(true) {
			// Step 2
			double checkPrimalFeasible = (matrixA.multiply(curPoint.getPointX()).subtract(vectorB)).getNorm(); //primal feasible
			double checkDualFeasible = ((matrixA.transpose().multiply(curPoint.getPointPi())).add(curPoint.getPointS().subtract(vectorC))).getNorm(); // dual feasible
			double checkDualityGap = curPoint.getPointS().transpose().multiply(curPoint.getPointX()).getEntry(0, 0); // duality gap

			if(checkPrimalFeasible < feasibleEpsilon && checkDualFeasible < feasibleEpsilon && checkDualityGap < epsilon) {
				break;
			}
			if(!primalFeasible && checkDualFeasible < GlobalParam.EPSILON1 && checkDualityGap < epsilon) {
				isPrimalFeasible = false;
				break;
			}
			if(iter>GlobalParam.BIGM3) { // Set an iteration limit on PDIPM
				isPrimalFeasible = false;
				break;
			}
			
			iter++;
			
			// Step 3
			double mu = (curPoint.getPointX().transpose().multiply(curPoint.getPointS()).getEntry(0, 0) / n );
			
			// Step 4. Determine affine direction
			PrimalDualPoint curGradAff = getNewtonDirectionFaster(matrixA, vectorB, vectorC, curPoint, 0);
			
			// Determine sigma
			double alphaPrimalAff = getStepSize(curPoint.getPointX(), curGradAff.getPointX());
			double alphaDualAff = getStepSize(curPoint.getPointS(), curGradAff.getPointS());
			PrimalDualPoint curPointAff = updateCurPoint(curPoint, curGradAff, alphaPrimalAff, alphaDualAff);
			double muAff = (curPointAff.getPointX().transpose().multiply(curPointAff.getPointS()).getEntry(0, 0) / n );
			double sigma = Math.pow(muAff / mu, 3);
			
			// Determine search direction
			PrimalDualPoint curGrad = getNewtonDirectionFaster(matrixA, vectorB, vectorC, curPoint, sigma*mu, true, curGradAff);
			
			// Step 5. Determine step size
			double alphaPrimal = getStepSize(curPoint.getPointX(), curGrad.getPointX());
			double alphaDual = getStepSize(curPoint.getPointS(), curGrad.getPointS());
			
			// Step 6. Update values
			curPoint = updateCurPoint(curPoint, curGrad, alphaPrimal, alphaDual);
		}
		solutionPrimalDualPoint = curPoint;
	}
	
	public boolean isPrimalFeasible() {
		return isPrimalFeasible;
	}

	public PrimalDualPoint getSolutionPrimalDualPoint() {
		return solutionPrimalDualPoint;
	}

	private PrimalDualPoint updateCurPoint(PrimalDualPoint curPoint, PrimalDualPoint curGrad, double alphaPrimal, double alphaDual) {
		RealMatrix newX = curPoint.getPointX().add(curGrad.getPointX().scalarMultiply(alphaPrimal));
		RealMatrix newPi = curPoint.getPointPi().add(curGrad.getPointPi().scalarMultiply(alphaDual));
		RealMatrix newS = curPoint.getPointS().add(curGrad.getPointS().scalarMultiply(alphaDual));
		
		return new PrimalDualPoint(newX, newPi, newS);
	}
	
	private double getStepSize(RealMatrix curPoint, RealMatrix curDirection) {
		double result = Double.POSITIVE_INFINITY;
		for(int i = 0; i < curPoint.getRowDimension(); i++) {
			if(curDirection.getEntry(i, 0)<0) {
				double temp = -1 * curPoint.getEntry(i, 0) / curDirection.getEntry(i, 0);
				if(temp < result) {
					result = temp;
				}
			}
		}
		return Math.min(1, r*result);
	}
	
	/**
	 * Extend matrixAAT
	 * 
	 * @param matrixAAT
	 * @param matrixA
	 * @return
	 */
	private void precomputeAAT_extend(PDCGM_matrixAAT[][] matrixAAT, RealMatrix matrixA, int startIndex) {
		for(int i = 0; i < matrixA.getRowDimension(); i++) {
			for(int j = i; j < matrixA.getRowDimension(); j++) {
				if(i==j) {
					List<Pair<Integer, Double>> tempMap = matrixAAT[i][i].getValues();
					for(int m = startIndex; m < matrixA.getColumnDimension(); m++) {
						double val = Math.pow(matrixA.getEntry(i, m), 2);
						if(val!=0) {
							tempMap.add(new Pair<>(m, val));
						}
					}
				}
				else {
					List<Pair<Integer, Double>> tempMap = matrixAAT[i][j].getValues();
					for(int m = startIndex; m < matrixA.getColumnDimension(); m++) {
						double val = matrixA.getEntry(i, m) * matrixA.getEntry(j, m);
						if(val!=0) {
							tempMap.add(new Pair<>(m, val));
						}
					}
				}
			}
		}
	}
	
	/**
	 * Compute LHS from scratch
	 * 
	 * We only have to compute part of the LHS once
	 * LHS = (AS^-1 X A^T)
	 * 
	 * Note that LHS is a matrix multiplication where the middle matrix is a diagonal
	 * The end result is a positive matrix
	 * @returnm
	 */
	private PDCGM_matrixAAT[][] precomputeAAT_fromScratch(RealMatrix matrixA) {
		PDCGM_matrixAAT[][] matrixAAT = new PDCGM_matrixAAT[matrixA.getRowDimension()][matrixA.getRowDimension()];
		for(int i = 0; i < matrixA.getRowDimension(); i++) {
			for(int j = i; j < matrixA.getRowDimension(); j++) {
				if(i==j) {
					List<Pair<Integer, Double>> tempMap = new ArrayList<>();
					for(int m = 0; m < matrixA.getColumnDimension(); m++) {
						double val = Math.pow(matrixA.getEntry(i, m), 2);
						if(val!=0) {
							tempMap.add(new Pair<>(m, val));
						}
					}
					matrixAAT[i][i] = new PDCGM_matrixAAT(tempMap);
				}
				else {
					List<Pair<Integer, Double>> tempMap = new ArrayList<>();
					for(int m = 0; m < matrixA.getColumnDimension(); m++) {
						double val = matrixA.getEntry(i, m) * matrixA.getEntry(j, m);
						if(val!=0) {
							tempMap.add(new Pair<>(m, val));
						}
					}
					matrixAAT[i][j] = new PDCGM_matrixAAT(tempMap);
				}
			}
		}
		return matrixAAT;
	}
	
	private PrimalDualPoint getNewtonDirectionFaster(RealMatrix matrixA, RealMatrix vectorB, RealMatrix vectorC, PrimalDualPoint curPoint,
			double theta) {
		return getNewtonDirectionFaster(matrixA, vectorB, vectorC, curPoint, theta, false, null);
	}
	
	/**
	 * Solve Equation 11
	 * 
	 * This should be slightly faster compared to solving Equation 10
	 * 
	 * delta_s = r_2 - A^T delta_pi
	 * delta_x = -x + theta S^-1 e - S^-1 X delta_s
	 * (AS^-1 X A^T) delta_pi = r_1 - A S^-1 (r_3 -X r_2) 
	 * 
	 * centerMehrotraCorrection:
	 * r3 := r3 - delta_x_aff delta_s_aff e => element wise product
	 * 
	 * @param matrixA
	 * @param vectorB
	 * @param vectorC
	 * @param curPoint
	 * @param theta
	 * @return
	 */
	private PrimalDualPoint getNewtonDirectionFaster(RealMatrix matrixA, RealMatrix vectorB, RealMatrix vectorC, PrimalDualPoint curPoint,
			double theta, boolean centerMehrotraCorrection, PrimalDualPoint gradAff) {
		RealMatrix r1 = vectorB.subtract(matrixA.multiply(curPoint.getPointX()));
		RealMatrix r2 = (vectorC.subtract(matrixA.transpose().multiply(curPoint.getPointPi()))).subtract(curPoint.getPointS());
		RealMatrix r3 = determineR3(curPoint, theta, centerMehrotraCorrection, gradAff);
		
		// Precalc
		RealMatrix S_inv = new Array2DRowRealMatrix(curPoint.getPointS().getRowDimension(), 1);
		double[] SX_positive2 = new double[curPoint.getPointS().getRowDimension()];
		for(int i = 0; i < curPoint.getPointS().getRowDimension(); i++) {
			double val = 1/curPoint.getPointS().getEntry(i, 0);			
			S_inv.setEntry(i, 0, val);
			val = (val*curPoint.getPointX().getEntry(i, 0)); 
			SX_positive2[i] = val;
		}
			
		// Calculate LHS
		// Note that LHS is a matrix multiplication where the middle matrix is a diagonal
		// The end result is a positive matrix
		RealMatrix lhs = new Array2DRowRealMatrix(matrixA.getRowDimension(), matrixA.getRowDimension());		
		for(int i = 0; i < matrixA.getRowDimension(); i++) {
			for(int j = i; j < matrixA.getRowDimension(); j++) {
				double val = 0;
				for(Pair<Integer, Double> entry2: matrixAAT_final[i][j].getValues()) {
					val += entry2.second * SX_positive2[entry2.first];
				}
				if(i==j) {
					val += GlobalParam.EPSILON3; // Avoid numerical issues
				}
				lhs.setEntry(i, j, val);
				lhs.setEntry(j, i, val);
			}
		}
		
		// Calculate RHS
		RealMatrix r3_prime = r3.subtract(MatrixOperations.elementWiseProduct(curPoint.getPointX(), r2)); // r_3 -X r_2	
		RealMatrix rhs = new Array2DRowRealMatrix(r1.getRowDimension(), 1);
		for(int i = 0; i < r1.getRowDimension(); i++) {
			double val2 = 0;
			for(int j = 0; j < matrixA.getColumnDimension(); j++) {
				val2 += matrixA.getEntry(i, j) * S_inv.getEntry(j, 0)* r3_prime.getEntry(j, 0);
			}
			double val = r1.getEntry(i, 0) - val2;
			rhs.setEntry(i, 0, val);
		}
		
		// Use cholesky decomposition. Note that the package incorporates both solves already
		// Positive Definite matrix B = LL^T
		// Solve Ly = rhs
		// Solve L^T x = y
		DecompositionSolver solver = null;
		try {
			CholeskyDecomposition cholSolver = new CholeskyDecomposition(lhs);
			solver = cholSolver.getSolver();
		} catch (org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException ex) {
		    
		}
		if (solver==null) {
			RealMatrix directionPi = new Array2DRowRealMatrix(curPoint.getPointPi().getRowDimension(), 1);
			RealMatrix directionS = new Array2DRowRealMatrix(curPoint.getPointS().getRowDimension(), 1);
			RealMatrix directionX = new Array2DRowRealMatrix(curPoint.getPointX().getRowDimension(), 1);
			return new PrimalDualPoint(directionX, directionPi, directionS);
		}
		
		// Solve pi
		RealMatrix directionPi = solver.solve(rhs);
		
		// Solve s
		RealMatrix directionS = r2.subtract(matrixA.transpose().multiply(directionPi));
		
		// Solve x
		RealMatrix directionX = S_inv.scalarMultiply(theta).subtract(curPoint.getPointX()).subtract((MatrixOperations.elementWiseProduct(curPoint.getPointX(), directionS, S_inv)));
		
		return new PrimalDualPoint(directionX, directionPi, directionS);
	}
	
	private RealMatrix determineR3(PrimalDualPoint curPoint, double theta, boolean centerMehrotraCorrection, PrimalDualPoint gradAff) {
		RealMatrix result = new Array2DRowRealMatrix(curPoint.getPointX().getRowDimension(), 1);
		for(int i = 0; i < curPoint.getPointX().getRowDimension(); i++) {
			double val = -curPoint.getPointX().getEntry(i, 0) * curPoint.getPointS().getEntry(i, 0) + theta;
			result.setEntry(i, 0, val);
		}		
		return result;
	}
	
	private void updateMatrix(RealMatrix largeMatrix, int startRow, int startCol, RealMatrix block) {
		for(int i = 0; i < block.getRowDimension(); i++) {
			for(int j = 0; j < block.getColumnDimension(); j++) {
				largeMatrix.setEntry(startRow+i, startCol+j, block.getEntry(i, j));	
			}
		}
	}
	
	
}
