package com.epam.tat.matrixprocessor.impl;

import com.epam.tat.matrixprocessor.IMatrixProcessor;

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
	public double[][] transpose(double[][] matrix) throws UnsupportedOperationException {
		if (transposeExamInputIsNotPassed(matrix)) {
			throw new UnsupportedOperationException("turnClockwise Failed");
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

	private boolean transposeExamInputIsNotPassed(double[][] matrix) {
		boolean examIsNotPassed = false;

		if (matrix == null) {
			examIsNotPassed = true;
		}
		if (matrix.length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix.length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix[0].length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix[0].length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}

		int numberOfEmptyValues = 0;

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				if (matrix[i][j] == 0) {
					numberOfEmptyValues++;
				}
			}
		}
		if (numberOfEmptyValues == (matrix.length + matrix[0].length)) {
			examIsNotPassed = true;
		}

		return examIsNotPassed;
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
	public double[][] turnClockwise(double[][] matrix) throws UnsupportedOperationException {
		if (turnClockwiseExamInputIsNotPassed(matrix)) {
			throw new UnsupportedOperationException("turnClockwise Failed");
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

	private boolean turnClockwiseExamInputIsNotPassed(double[][] matrix) {
		boolean examIsNotPassed = false;

		if (matrix == null) {
			examIsNotPassed = true;
		}
		if (matrix.length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix.length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix[0].length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix[0].length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}


		int numberOfEmptyValues = 0;

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				if (matrix[i][j] == 0) {
					numberOfEmptyValues++;
				}
			}
		}
		if (numberOfEmptyValues == (matrix.length + matrix[0].length)) {
			examIsNotPassed = true;
		}

		return examIsNotPassed;
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
	public double[][] multiplyMatrices (double[][] firstMatrix, double[][] secondMatrix) throws UnsupportedOperationException {
		if (multiplyMatricesExamInput(firstMatrix, secondMatrix)) {
			throw new UnsupportedOperationException("multiplyMatrices Failed");
		}

		double[][] multiplyResultMatrix = new double[firstMatrix.length][secondMatrix[0].length];

		for (int i = 0; i < multiplyResultMatrix[0].length; i++) {
			for (int j = 0; j < multiplyResultMatrix.length; j++) {
				for (int k = 0; k < firstMatrix[0].length; k++) {
					multiplyResultMatrix[i][j] = roundValue(multiplyResultMatrix[i][j] + roundValue(firstMatrix[i][k]) * roundValue(secondMatrix [k][j]));
				}
			}
		}

		return multiplyResultMatrix;
	}

	private boolean multiplyMatricesExamInput(double[][] firstMatrix, double[][] secondMatrix) {
		boolean examIsNotPassed = false;

		if (firstMatrix[0].length != secondMatrix.length) {
			examIsNotPassed = true;
		}
		if (firstMatrix.length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (firstMatrix.length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		if (firstMatrix[0].length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (firstMatrix[0].length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		if (secondMatrix.length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (secondMatrix.length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		if (secondMatrix[0].length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (secondMatrix[0].length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		return examIsNotPassed;
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
	public double[][] getInverseMatrix(double[][] matrix) throws UnsupportedOperationException {
		if (getInverseMatrixExamInput(matrix)) {
			throw new UnsupportedOperationException("getInverseMatrix Failed");
		}

		double[][] minorMatrix = new double[matrix.length][matrix[0].length];
		double determinant = getMatrixDeterminant(matrix);

        if (determinant == 0) {
            throw new UnsupportedOperationException("getInverseMatrix Failed. Determinant = 0");
        }

		for (int i = 0; i < matrix.length; i++) {
			for (int j=0; j < matrix[0].length; j++) {
				double[][] minor = getMinorMatrix(matrix, i, j);
				minorMatrix[i][j] = getMatrixDeterminant(minor);
			}
		}

		double[][] algAdd = new double[matrix.length][matrix[0].length];
		int signBeforeNumber = 1;

		for (int i = 0; i < matrix.length; i++) {
			if (i %2 == 0) {
				signBeforeNumber = 1;
			} else {
				signBeforeNumber = -1;
			}
			for (int j = 0; j < matrix.length; j++) {
				algAdd[i][j] = signBeforeNumber * minorMatrix[i][j];
                signBeforeNumber =  signBeforeNumber * -1;
			}
		}

		double[][] algAddTranspose = transpose(algAdd);
		double[][] inverseMatrix = new double[matrix.length][matrix[0].length];

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
                inverseMatrix[i][j] = roundValue(1 / determinant * algAddTranspose[i][j]);
			}
		}

		return inverseMatrix;
	}

	private boolean getInverseMatrixExamInput(double[][] matrix) {
		boolean examIsNotPassed = false;

		if (matrix == null) {
			examIsNotPassed = true;
		}
		if (matrix.length != matrix[0].length) {
			examIsNotPassed = true;
		}
		if (matrix.length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix.length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix[0].length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix[0].length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		return examIsNotPassed;
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
	public double getMatrixDeterminant(double[][] matrix) throws UnsupportedOperationException {
		if (getMatrixDeterminantExamInput(matrix)) {
			throw new UnsupportedOperationException("getMatrixDeterminant Failed");
		}

		if (matrix.length == 1) {
			return matrix[0][0];
		}

		double determinant = 0;
		double[][] matrixMinor = new double[matrix.length-1][matrix[0].length-1];
		int signBeforeNumber = 1;
        int x = 0;
        int y = 0;

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
			determinant += signBeforeNumber * roundValue(matrix[0][i]) * roundValue(getMatrixDeterminant(matrixMinor));
			signBeforeNumber *= (-1);
		}
        //Arrays.stream(inverseMatrix).map(Arrays::toString).forEach(System.out::println);

		return roundValue(determinant);
	}

	private boolean getMatrixDeterminantExamInput(double[][] matrix) {
		boolean examIsNotPassed = false;

		if (matrix[0].length != matrix.length) {
			examIsNotPassed = true;
		}
		if (matrix == null) {
			examIsNotPassed = true;
		}
		if (matrix.length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix.length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix[0].length <= MIN_DIMENSION) {
			examIsNotPassed = true;
		}
		if (matrix[0].length >= MAX_DIMENSION) {
			examIsNotPassed = true;
		}
		return examIsNotPassed;
	}

	public static double[][] getMinorMatrix(double[][] matrix, int rowDeleted, int columnDeleted) {

		final int ROWS_MINOR_MATRIX = matrix.length - 1;
		final int COLUMNS_MINOR_MATRIX = matrix[0].length - 1;
		double[][] resultMinorMatrix = new double[ROWS_MINOR_MATRIX][COLUMNS_MINOR_MATRIX];
		int i;
        int j;
        int x = 0;
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
        bigDecimal = bigDecimal.setScale(PRECISION, RoundingMode.HALF_UP);
        return bigDecimal.doubleValue();
    }
}
