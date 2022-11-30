package com.epam.tat.matrixprocessor.impl;

import com.epam.tat.matrixprocessor.IMatrixProcessor;

import java.util.Arrays;

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

		if (matrix.length <= 0) {
			throw new UnsupportedOperationException("Number of rows in matrix <= 0.\n");
		}

		if (matrix.length >= 10) {
			throw new UnsupportedOperationException("Number of rows in matrix >= 10.\n");
		}

		if (matrix[0].length <= 0) {
			throw new UnsupportedOperationException("Number of columns in matrix <= 0.\n");
		}

		if (matrix[0].length >= 10) {
			throw new UnsupportedOperationException("Number of columns in matrix >= 10.\n");
		}

		final int ROWS = matrix[0].length;
		final int COLUMNS = matrix.length;
		double[][] transposeMatrix = new double[ROWS][COLUMNS];

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				transposeMatrix[j][i] = matrix[i][j];
			}
		}

		System.out.println("transpose");
		Arrays.stream(matrix).map(Arrays::toString).forEach(System.out::println);
		System.out.println(" ");
		Arrays.stream(transposeMatrix).map(Arrays::toString).forEach(System.out::println);

	return transposeMatrix;
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

		if (matrix.length <= 0) {
			throw new UnsupportedOperationException("Number of rows in matrix <= 0.\n");
		}

		if (matrix.length >= 10) {
			throw new UnsupportedOperationException("Number of rows in matrix >= 10.\n");
		}

		if (matrix[0].length <= 0) {
			throw new UnsupportedOperationException("Number of columns in matrix <= 0.\n");
		}

		if (matrix[0].length >= 10) {
			throw new UnsupportedOperationException("Number of columns in matrix >= 10.\n");
		}

		final int ROWS = matrix[0].length;
		final int COLUMNS = matrix.length;
		double[][] turnClockwiseMatrix = new double[ROWS][COLUMNS];

		for (int i = 0; i < COLUMNS; i++) {
			for (int j = 0; j < ROWS; j++) {
				turnClockwiseMatrix[j][COLUMNS-i-1] = matrix[i][j];
			}
		}

		System.out.println("turnClockwise");
		Arrays.stream(matrix).map(Arrays::toString).forEach(System.out::println);
		System.out.println(" ");
		Arrays.stream(turnClockwiseMatrix).map(Arrays::toString).forEach(System.out::println);


		return turnClockwiseMatrix;
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

		if (firstMatrix[0].length != secondMatrix.length || firstMatrix[0].length == 0 || secondMatrix.length == 0) {
			throw new UnsupportedOperationException("Multiplying error: wrong demention. " +
					"Number of columns first matrix is " + firstMatrix[0].length + "\n. " +
					"Number of rows second matrix is " + secondMatrix.length +"\n");
		}

		if (firstMatrix.length <= 0) {
			throw new UnsupportedOperationException("Number of rows in first matrix <= 0.\n");
		}

		if (firstMatrix.length >= 10) {
			throw new UnsupportedOperationException("Number of rows in first matrix >= 10.\n");
		}

		if (firstMatrix[0].length <= 0) {
			throw new UnsupportedOperationException("Number of columns in first matrix <= 0.\n");
		}

		if (firstMatrix[0].length >= 10) {
			throw new UnsupportedOperationException("Number of columns in first matrix >= 10.\n");
		}

		if (secondMatrix.length <= 0) {
			throw new UnsupportedOperationException("Number of rows in second matrix <= 0.\n");
		}

		if (secondMatrix.length >= 10) {
			throw new UnsupportedOperationException("Number of rows in second matrix >= 10.\n");
		}

		if (secondMatrix[0].length <= 0) {
			throw new UnsupportedOperationException("Number of columns in second matrix <= 0.\n");
		}

		if (secondMatrix[0].length >= 10) {
			throw new UnsupportedOperationException("Number of columns in second matrix >= 10.\n");
		}

		double[][] multiplyResultMatrix = new double[firstMatrix.length][secondMatrix[0].length];

		for (int i = 0; i < multiplyResultMatrix[0].length; i++) {
			for (int j = 0; j < multiplyResultMatrix.length; j++) {
				for (int k = 0; k < firstMatrix[0].length; k++) {
					multiplyResultMatrix[i][j] = multiplyResultMatrix[i][j] + firstMatrix[i][k] * secondMatrix [k][j];
				}
			}
		}

		System.out.println("multiplyMatrices");
		Arrays.stream(firstMatrix).map(Arrays::toString).forEach(System.out::println);
		System.out.println(" ");
		Arrays.stream(secondMatrix).map(Arrays::toString).forEach(System.out::println);
		System.out.println(" ");
		Arrays.stream(multiplyResultMatrix).map(Arrays::toString).forEach(System.out::println);


		return multiplyResultMatrix;
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
		if (matrix.length != matrix[0].length) {
			throw new UnsupportedOperationException("Error. Wrong dimension. Rows not equals columns\n");
		}

		if (matrix.length <= 0) {
			throw new UnsupportedOperationException("Number of rows in matrix <= 0.\n");
		}

		if (matrix.length >= 10) {
			throw new UnsupportedOperationException("Number of rows in matrix >= 10.\n");
		}

		if (matrix[0].length <= 0) {
			throw new UnsupportedOperationException("Number of columns in matrix <= 0.\n");
		}

		if (matrix[0].length >= 10) {
			throw new UnsupportedOperationException("Number of columns in matrix >= 10.\n");
		}


		double[][] inverseMatrix = new double[matrix.length][matrix[0].length];

		double determinant = getMatrixDeterminant(matrix);
		double minor[][] = getMatrixMinor(matrix);

		return inverseMatrix;
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

		if (matrix.length <= 0) {
			throw new UnsupportedOperationException("Number of rows in matrix <= 0.\n");
		}
		if (matrix.length >= 10) {
			throw new UnsupportedOperationException("Number of rows in matrix >= 10.\n");
		}
		if (matrix[0].length <= 0) {
			throw new UnsupportedOperationException("Number of columns in matrix <= 0.\n");
		}
		if (matrix[0].length >= 10) {
			throw new UnsupportedOperationException("Number of columns in matrix >= 10.\n");
		}
		if (matrix.length == 1) {
			return matrix[0][0];
		}
		double determinant = 0;

		double matrixMinor[][] = new double[matrix.length-1][matrix[0].length-1];
		int signBeforeNumber = 1;

		for (int i = 0; i < matrix.length; i++) {
			int x = 0, y = 0;
			for (int j = 1; j < matrix.length; j++) {
				for(int k = 0; k < matrix.length; k++) {
					if (i == k) {
						continue;
					}
					matrixMinor[x][y] = matrix[j][k];
					y++;
					if (y == (matrix.length - 1)) {
						y = 0;
						x++;
					}
				}
			}
			determinant += signBeforeNumber * matrix[0][i] * getMatrixDeterminant(matrixMinor);
			signBeforeNumber *= (-1);
		}
		return determinant;
	}
public static double[][] getMatrixMinor(double[][] matrix) {

		int n = matrix.length - 1;
		int m = matrix[0].length - 1;
		int rowNum = matrix.length;
		int colNum = matrix[0].length;

		double[][] result = new double[ n ][ m ];

		for (int i = 0; i < matrix.length; i++) {
			boolean isRowDeleted = rowNum < i;
			int resultRowIndex = isRowDeleted ? i - 1 : i;
			for (int j = 0; j < matrix[i].length; j++) {
				boolean isColDeleted = colNum < j;
				int resultColIndex = isColDeleted ? j - 1 : j;
				if (rowNum != i && colNum != j)
					result[resultRowIndex][resultColIndex] = matrix[i][j];
			}
		}
		return result;
	}
}
