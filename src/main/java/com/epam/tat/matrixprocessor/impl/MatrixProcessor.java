package com.epam.tat.matrixprocessor.impl;

import com.epam.tat.matrixprocessor.IMatrixProcessor;

import java.sql.SQLOutput;
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
	public double[][] getInverseMatrix(double[][] matrix) {
		throw new UnsupportedOperationException("You need to implement this method");
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
	public double getMatrixDeterminant(double[][] matrix) {
		if (matrix.length != matrix.length || matrix[0].length == 0 || matrix.length == 0) {
			throw new UnsupportedOperationException("Determinant error: wrong demention. " +
					"Number of columns matrix is " + matrix[0].length + "\n. " +
					"Number of rows matrix is " + matrix.length +"\n");
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

		int[] index = new int[matrix.length];
		int numberOfValues = 0;
		double k;
		double[][] matrixA = new double[matrix.length][matrix[0].length];
		double[] vectorB = new double[matrix.length];
		double determinant = 1;

		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				matrixA[i][j] = matrix[i][j];
			}
		}

		for (int i = 0; i < matrix.length; i++) {
			index[i] = i;
		}

		for (int p = 0; p < (matrix.length - 1); p++) {
			gaussMethod(matrixA, vectorB, matrix.length, p, index, numberOfValues);
			for (int i = p+1; i < matrix.length; i++) {
				if (matrixA[p][p] != 0) {
					k = matrixA[i][p]/matrixA[p][p];
					for ( int j = p; j < matrix.length; j++) {
						matrixA[i][j] = matrixA[i][j] - matrixA[p][j] * k;
					}
				}
			}
		}


		for (int i = 0; i < matrix.length; i++) {
			determinant *= matrixA[i][i];
		}

		determinant = Math.pow(numberOfValues, -1) * determinant;

		return determinant;
	}

	private void gaussMethod(double[][] matrixA, double[] vectorB, int length, int p, int[] index, int numberOfValues) {
		double maxElement;
		double tempValue;
		int indexP;
		int maxElementIndex_i = p;
		int maxElementIndex_j = p;

		maxElement = Math.abs(matrixA[p][p]);

		for (int i = p; i < length; i++) {
			for (int j = p; j < length; j++) {
				if (Math.abs(matrixA[i][j]) > maxElement) {
					maxElement = Math.abs(matrixA[i][j]);
					maxElementIndex_i = i;
					maxElementIndex_j = j;
				}
			}
		}

		if (maxElementIndex_i != p) {
			for (int i = 0; i < length; i++) {
				tempValue = matrixA[p][i];
				matrixA[p][i] = matrixA[maxElementIndex_i][i];
				matrixA[maxElementIndex_i][i] = tempValue;
			}

			tempValue = vectorB[p];
			vectorB[p] = vectorB[maxElementIndex_i];
			vectorB[maxElementIndex_i] = tempValue;
			numberOfValues++;
		}

		if (maxElementIndex_j != p) {
			for (int i = 0; i < length; i++) {
				tempValue = matrixA[i][p];
				matrixA[i][p] = matrixA[i][maxElementIndex_j];
				matrixA[i][maxElementIndex_j] = tempValue;
			}
			indexP = index[p];
			index[p] = index[maxElementIndex_j];
			index[maxElementIndex_j] = indexP;
			numberOfValues++;
		}
	}
}
