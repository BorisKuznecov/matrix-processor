package com.epam.tat.matrixprocessor.impl;

import com.epam.tat.matrixprocessor.IMatrixProcessor;
import com.epam.tat.matrixprocessor.exception.MatrixProcessorException;

import java.math.BigDecimal;
import java.math.RoundingMode;

/**
 * Function Description:
 * Complete the functions below. All methods must work with matrices of the double type.
 *
 * Constraints:
 * 0 < m < 10
 * 0 < n < 10
 * where m - number of rows in matrix
 * where n - number of columns in matrix
 *
 * In case of incorrect input values or inability to perform a calculation, the method should throw an appropriate
 * exception.
 *
 */
public class MatrixProcessor implements IMatrixProcessor {

	public static final int MAX_DIMENSION = 10;
	public static final int MIN_DIMENSION = 0;
	public static final int PRECISION = 3;

	/**
	 * Matrix transpose is an operation on a matrix where its rows become columns with the same numbers.
	 * Ex.:
	 * |1 2|			|1 3 5|
	 * |3 4|   ====>	|2 4 6|
	 * |5 6|
	 *
	 * @param matrix - matrix for transposition
	 * @return the transposed matrix
	 */
	@Override
	public double[][] transpose(double[][] matrix) throws MatrixProcessorException {
		if (transposeExamInput(matrix)) {
			throw new MatrixProcessorException("Transpose exam input failed");
		}

		final int ROWS_TRANSPOSE_MATRIX = matrix[0].length;
		final int COLUMNS_TRANSPOSE_MATRIX = matrix.length;
		double[][] transposeMatrix = new double[ROWS_TRANSPOSE_MATRIX][COLUMNS_TRANSPOSE_MATRIX];

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				transposeMatrix[j][i] = matrix[i][j];
			}
		}

		return transposeMatrix;
	}

	private boolean transposeExamInput(double[][] matrix) {
		return matrix == null
				|| matrix.length == MIN_DIMENSION
				|| matrix.length >= MAX_DIMENSION
				|| matrix[0].length == MIN_DIMENSION
				|| matrix[0].length >= MAX_DIMENSION;
	}

	/**
	 * The method flips the matrix clockwise.
	 * Ex.:
	 * * |1 2|			|5 3 1|
	 * * |3 4|   ====>	|6 4 2|
	 * * |5 6|
	 *
	 * @param matrix - rotation matrix
	 * @return rotated matrix
	 */
	@Override
	public double[][] turnClockwise(double[][] matrix) throws MatrixProcessorException {
		if (turnClockwiseExamInput(matrix)) {
			throw new MatrixProcessorException("Turn clockwise exam input failed");
		}

		final int ROWS_CLOCKWISE_MATRIX = matrix[0].length;
		final int COLUMNS_CLOCKWISE_MATRIX = matrix.length;
		double[][] turnClockwiseMatrix = new double[ROWS_CLOCKWISE_MATRIX][COLUMNS_CLOCKWISE_MATRIX];

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				turnClockwiseMatrix[j][COLUMNS_CLOCKWISE_MATRIX-i-1] = matrix[i][j];
			}
		}

		return turnClockwiseMatrix;
	}

	private boolean turnClockwiseExamInput(double[][] matrix) {
		return matrix == null
				|| matrix.length == MIN_DIMENSION
				|| matrix.length >= MAX_DIMENSION
				|| matrix[0].length == MIN_DIMENSION
				|| matrix[0].length >= MAX_DIMENSION;
	}

	/**
	 * This method multiplies matrix firstMatrix by matrix secondMatrix
	 *
	 * See {https://en.wikipedia.org/wiki/Matrix_multiplication}
	 *
	 * @param firstMatrix  - first matrix to multiply
	 * @param secondMatrix - second matrix to multiply
	 * @return result matrix
	 */
	@Override
	public double[][] multiplyMatrices (double[][] firstMatrix, double[][] secondMatrix) throws MatrixProcessorException {
		if (multiplyMatricesExamInput(firstMatrix, secondMatrix)) {
			throw new MatrixProcessorException("Multiply matrices exam input failed");
		}

		double[][] multiplyResultMatrix = new double[firstMatrix.length][secondMatrix[0].length];

		for (int i = 0; i < firstMatrix.length; i++) {
			for (int j = 0; j < secondMatrix[0].length; j++) {
				for (int k = 0; k < secondMatrix.length; k++) {
					multiplyResultMatrix[i][j] = roundValue(multiplyResultMatrix[i][j]
							+ (firstMatrix[i][k]) * (secondMatrix [k][j]));
				}
			}
		}

		return multiplyResultMatrix;
	}

	private boolean multiplyMatricesExamInput(double[][] firstMatrix, double[][] secondMatrix) {
		return firstMatrix[0].length != secondMatrix.length
				|| firstMatrix == null
				|| secondMatrix == null // testSecondMultiplyMatrixIsNull ????
				|| firstMatrix.length == MIN_DIMENSION /*testMultiplyMatricesAreEmpty, testFirstMultiplyMatrixIsNull*/
				|| firstMatrix.length >= MAX_DIMENSION /*testMultiplyMatricesAreEmpty*/
				|| firstMatrix[0].length == MIN_DIMENSION
				|| firstMatrix[0].length >= MAX_DIMENSION
				|| secondMatrix.length == MIN_DIMENSION
				|| secondMatrix.length >= MAX_DIMENSION
				|| secondMatrix[0].length == MIN_DIMENSION
				|| secondMatrix[0].length >= MAX_DIMENSION;
	}

	/**
	 * This method returns the inverse of the matrix
	 *
	 * See {https://en.wikipedia.org/wiki/Invertible_matrix}
	 *
	 * @param matrix - input matrix
	 * @return inverse matrix for input matrix
	 */
	@Override
	public double[][] getInverseMatrix(double[][] matrix) throws MatrixProcessorException {
		if (getInverseMatrixExamInput(matrix)) {
			throw new MatrixProcessorException("Inverse matrix exam input failed");
		}

		double determinant = getMatrixDeterminant(matrix);

		if (determinant == 0) {
			throw new MatrixProcessorException("Determinant = 0");
		}

		double[][] minorMatrix = new double[matrix.length][matrix[0].length];

		for (int i = 0; i < matrix.length; i++) {
			for (int j=0; j < matrix[0].length; j++) {
				double[][] minor = getMinorMatrix(matrix, i, j);
				minorMatrix[i][j] = getMatrixDeterminant(minor);
			}
		}

		double[][] algAddMatrix = getAlgAddMatrix(minorMatrix);
		double[][] algAddTranspose = transpose(algAddMatrix);
		double[][] inverseMatrix = new double[matrix.length][matrix[0].length];

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				inverseMatrix[i][j] = roundValue(1 / determinant * algAddTranspose[i][j]);
			}
		}

		return inverseMatrix;
	}

	private boolean getInverseMatrixExamInput(double[][] matrix) {
		return matrix == null
				|| matrix[0].length != matrix.length
				|| matrix.length >= MAX_DIMENSION
				|| matrix.length == MIN_DIMENSION;//testInverseMatrixInputIsEmpty
	}

	private double[][] getAlgAddMatrix(double[][] minorMatrix) {
		double[][] algAddMatrix = new double[minorMatrix.length][minorMatrix[0].length];
		int signBeforeNumber;

		for (int i = 0; i < minorMatrix.length; i++) {
			if (i %2 == 0) {
				signBeforeNumber = 1;
			} else {
				signBeforeNumber = -1;
			}
			for (int j = 0; j < minorMatrix.length; j++) {
				algAddMatrix[i][j] = signBeforeNumber * minorMatrix[i][j];
				signBeforeNumber =  signBeforeNumber * -1;
			}
		}

		return algAddMatrix;
	}

	/**
	 * This method returns the determinant of the matrix
	 *
	 * See {https://en.wikipedia.org/wiki/Determinant}
	 *
	 * @param matrix - input matrix
	 * @return determinant of input matrix
	 */
	@Override
	public double getMatrixDeterminant(double[][] matrix) throws MatrixProcessorException {
		if (getMatrixDeterminantExamInput(matrix)) {
			throw new MatrixProcessorException("Matrix determinant input exam failed");
		}

		if (matrix.length == 1) {
			return matrix[0][0];
		}

		double determinant = 0;
		double[][] matrixMinor = new double[matrix.length-1][matrix[0].length-1];
		int signBeforeNumber = 1;
		int x;
		int y;

		for (int i = 0; i < matrix.length; i++) {
			x = 0;
			y = 0;
			for (int j = 1; j < matrix.length; j++) {
				for(int k = 0; k < matrix.length; k++) {
					if (i == k) {
						continue;
					}
					matrixMinor[x][y] = roundValue(matrix[j][k]);
					y++;
					if (y == (matrix.length - 1)) {
						y = 0;
						x++;
					}
				}
			}
			determinant += roundValue(signBeforeNumber * matrix[0][i] * getMatrixDeterminant(matrixMinor));
			signBeforeNumber *= (-1);
		}

		return roundValue(determinant);
	}

	private boolean getMatrixDeterminantExamInput(double[][] matrix) {
		return matrix == null
				|| matrix[0].length != matrix.length
				|| matrix.length >= MAX_DIMENSION
				|| matrix.length == MIN_DIMENSION;//testMatrixDeterminantInputIsEmpty
	}

	public static double[][] getMinorMatrix(double[][] matrix, int rowDeleted, int columnDeleted) {
		final int ROWS_MINOR_MATRIX = matrix.length - 1;
		final int COLUMNS_MINOR_MATRIX = matrix[0].length - 1;
		double[][] resultMinorMatrix = new double[ROWS_MINOR_MATRIX][COLUMNS_MINOR_MATRIX];
		int i;
		int j;
		int x;
		int y = 0;

		for (i = 0; i < matrix.length; i++) {
			x = 0;
			for (j = 0; j < matrix[0].length; j++) {
				if (i != rowDeleted && j != columnDeleted) {
					resultMinorMatrix[x][y] = matrix[i][j];
					x++;
				}
			}
			if (i != rowDeleted && j != columnDeleted) {
				y++;
			}
		}

		return resultMinorMatrix;
	}

	private static double roundValue(double value) {
		if (PRECISION < 0) throw new IllegalArgumentException();
		BigDecimal bigDecimal = new BigDecimal(Double.toString(value));
		bigDecimal = bigDecimal.setScale(PRECISION, RoundingMode.HALF_DOWN);
		return bigDecimal.doubleValue();
	}
}