package cl.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cl.data.GlobalParam;

public class MatrixOperations {
	private static Logger log = LoggerFactory.getLogger(MatrixOperations.class);
	
	/**
	 * Multiply the elements of two row vectors
	 */
	public static RealMatrix elementWiseProduct(RealMatrix matrix1, RealMatrix matrix2) {
		RealMatrix result = new Array2DRowRealMatrix(matrix1.getRowDimension(), 1);
		for(int i = 0; i < matrix1.getRowDimension(); i++) {
			result.setEntry(i, 0, matrix1.getEntry(i, 0) * matrix2.getEntry(i, 0));
		}
		return result;
	}
	
	/**
	 * Multiply the elements of three row vectors
	 */
	public static RealMatrix elementWiseProduct(RealMatrix matrix1, RealMatrix matrix2, RealMatrix matrix3) {
		RealMatrix result = new Array2DRowRealMatrix(matrix1.getRowDimension(), 1);
		for(int i = 0; i < matrix1.getRowDimension(); i++) {
			result.setEntry(i, 0, matrix1.getEntry(i, 0) * matrix2.getEntry(i, 0) * matrix3.getEntry(i, 0));
		}
		return result;
	}
	
	/**
	 * Add a row to the end of the vector
	 * @param matrix
	 * @param vector
	 * @return
	 */
	public static RealMatrix addRow(RealMatrix matrix, double element, int nElements) {
		RealMatrix result = new Array2DRowRealMatrix(matrix.getRowDimension()+nElements, 1);
		result.setSubMatrix(matrix.getData(), 0, 0);
		for(int i = 0; i < nElements; i++) {
			result.setEntry(matrix.getRowDimension()+i, 0, element);
		}
		return result;
	}
	
	/**
	 * Add a row to the end of the vector
	 * @param matrix
	 * @param vector
	 * @return
	 */
	public static RealMatrix addRow(RealMatrix matrix, double element) {
		RealMatrix result = new Array2DRowRealMatrix(matrix.getRowDimension()+1, 1);
		result.setSubMatrix(matrix.getData(), 0, 0);
		result.setEntry(matrix.getRowDimension(), 0, element);
		return result;
	}
	
	/**
	 * Add a row vector to the end of the matrix
	 * @param matrix
	 * @param vector
	 * @return
	 */
	public static RealMatrix addRow(RealMatrix matrix, RealMatrix vector) {
		RealMatrix result = new Array2DRowRealMatrix(matrix.getRowDimension()+1, matrix.getColumnDimension());
		result.setSubMatrix(matrix.getData(), 0, 0);
		result.setRow(matrix.getRowDimension(), vector.getColumn(0));
		return result;
	}
	
	/**
	 * Add a column vector to the end of the matrix
	 * @param matrix
	 * @param vector
	 * @return
	 */
	public static RealMatrix addColumn(RealMatrix matrix, RealMatrix vector) {
		RealMatrix result = new Array2DRowRealMatrix(matrix.getRowDimension(), matrix.getColumnDimension()+1);
		result.setSubMatrix(matrix.getData(), 0, 0);
		result.setColumn(matrix.getColumnDimension(), vector.getColumn(0));
		return result;
	}
	
	
	/**
	 * From a column/row vector to diagonal matrix
	 * @param matrix
	 * @return
	 */
	public static RealMatrix diagonal(RealMatrix matrix) {
		if(matrix.getRowDimension()!=1 && matrix.getColumnDimension()!=1) {
			throw new IllegalArgumentException("Input should be a vector");
		}
		int newDimension = matrix.getRowDimension();
		boolean useRow = true;
		if(matrix.getRowDimension()<matrix.getColumnDimension()) {
			newDimension = matrix.getColumnDimension();
			useRow = false;
		}


		RealMatrix newMatrix = new Array2DRowRealMatrix(newDimension, newDimension);
		for(int i = 0; i < newDimension; i++) {
			if(useRow) {
				newMatrix.setEntry(i, i, matrix.getEntry(i, 0));
			}
			else {
				newMatrix.setEntry(i, i, matrix.getEntry(0, i));
			}
		}

		return newMatrix;
	}

	/**
	 * From a diagonal matrix to a column
	 * @param matrix
	 * @return
	 */
	public static RealMatrix selectDiagonalElements(RealMatrix matrix) {
		if(matrix.getRowDimension()!=matrix.getColumnDimension()) {
			throw new IllegalArgumentException("Input should be a square matrix");
		}

		RealMatrix newMatrix = new Array2DRowRealMatrix(matrix.getRowDimension(), 1);

		for(int i = 0; i < matrix.getRowDimension(); i++) {
			newMatrix.setEntry(i, 0, matrix.getEntry(i, i));
		}
		return newMatrix;
	}

	public static RealMatrix toThePower(RealMatrix matrix, double pow) {
		RealMatrix newMatrix = new Array2DRowRealMatrix(matrix.getRowDimension(), matrix.getColumnDimension());

		for(int i = 0; i < matrix.getRowDimension(); i++) {
			for(int j = 0; j < matrix.getColumnDimension(); j++) {
				newMatrix.setEntry(i, j, Math.pow(matrix.getEntry(i, j), pow));
			}
		}
		return newMatrix;
	}

	public static double sumAllElements(RealMatrix matrix) {
		double sum = 0;
		for(int i = 0; i < matrix.getRowDimension(); i++) {
			for(int j = 0; j < matrix.getColumnDimension(); j++) {
				sum += matrix.getEntry(i, j);
			}
		}
		return sum;
	}
	
	/**
	 * Sum elements vertically
	 * Get column sum
	 * @param matrix
	 * @return
	 */
	public static RealMatrix sumElementsVertically(RealMatrix matrix) {
		double[][] sum = new double[matrix.getColumnDimension()][1];
		for(int i = 0; i < matrix.getColumnDimension(); i++) {
			for(int j = 0; j < matrix.getRowDimension(); j++) {
				sum[i][0] += matrix.getEntry(j, i);
			}
		}
		return MatrixUtils.createRealMatrix(sum);
	}
	
	/**
	 * Sum elements horizontally
	 * Get row sum
	 * @param matrix
	 * @return
	 */
	public static RealMatrix sumElementsHorizontally(RealMatrix matrix) {
		double[][] sum = new double[matrix.getRowDimension()][1];
		for(int j = 0; j < matrix.getRowDimension(); j++) {
			for(int i = 0; i < matrix.getColumnDimension(); i++) {
				sum[j][0] += matrix.getEntry(j, i);
			}
		}
		return MatrixUtils.createRealMatrix(sum);
	}

	public static RealMatrix roundToNearestInteger(RealMatrix matrix) {
		RealMatrix newMatrix = new Array2DRowRealMatrix(matrix.getRowDimension(), matrix.getColumnDimension());
		for(int i = 0; i < matrix.getRowDimension(); i++) {
			for(int j = 0; j < matrix.getColumnDimension(); j++) {
				newMatrix.setEntry(i, j, Math.round(matrix.getEntry(i, j)));
			}
		}
		return newMatrix;
	}
	
	/**
	 * Multiplies entries of two matrices. Two matrices must have the same size
	 * @param matrix1
	 * @param matrix2
	 * @return
	 */
	public static RealMatrix multiplyElements(RealMatrix matrix1, RealMatrix matrix2) {
		if(matrix1.getRowDimension()!=matrix2.getRowDimension() || matrix1.getColumnDimension()!=matrix2.getColumnDimension()) {
			throw new IllegalArgumentException("Input matrices should have the same size");
		}
		RealMatrix newMatrix = new Array2DRowRealMatrix(matrix1.getRowDimension(), matrix1.getColumnDimension());
		for(int i = 0; i < newMatrix.getRowDimension(); i++) {
			for(int j = 0; j < newMatrix.getColumnDimension(); j++) {
				newMatrix.setEntry(i, j, matrix1.getEntry(i, j)*matrix2.getEntry(i, j));
			}
		}
		return newMatrix;
	}
	
	public static void setRandomElementsMatrix(RealMatrix matrix, double val) {
		for(int i = 0; i < matrix.getRowDimension(); i++) {
			for(int j = 0; j < matrix.getColumnDimension(); j++) {
				matrix.setEntry(i, j, val*GlobalParam.RANDOM.nextDouble());
			}
		}
	}
	
	/**
	 * Set all elements of matrix equal to input
	 * Remaining elements are assigned val
	 * @param matrix
	 * @param input
	 * @param val
	 */
	public static void setElementsMatrix(RealMatrix matrix, RealMatrix inputMatrix, double val) {
		for(int i = 0; i < matrix.getRowDimension(); i++) {
			for(int j = 0; j < matrix.getColumnDimension(); j++) {
				if(i < inputMatrix.getRowDimension() && j < inputMatrix.getColumnDimension()) {
					matrix.setEntry(i, j, inputMatrix.getEntry(i, j));
				}
				else {
					matrix.setEntry(i, j, val);
				}
			}
		}
	}

	public static void printMatrix(RealMatrix matrix) {
		if(matrix==null) {
			throw new IllegalArgumentException("Null");
		}
		for(int i = 0; i < matrix.getRowDimension(); i++) {
			for(int j = 0; j < matrix.getColumnDimension(); j++) {
				System.out.print(matrix.getEntry(i, j) + " ");
			}
			System.out.println();
		}
		System.out.println();
	}

	public static void writeMatrixToTxt(RealMatrix matrix, String name) {
		if(matrix==null) {
			throw new IllegalArgumentException("Null");
		}
		try {
			String fileName = "data/matrixF/matrixF_"+name+".txt";
			PrintWriter pw = new PrintWriter(new File(fileName));
			for(int i = 0; i < matrix.getRowDimension(); i++) {
				for(int j = 0; j < matrix.getColumnDimension(); j++) {
					pw.print(matrix.getEntry(i, j)+ " ");
				}
				pw.println();
			}
			pw.close();
			log.info("Write matrix F to {}", fileName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

}