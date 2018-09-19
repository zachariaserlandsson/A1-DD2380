package main.java.helpers;

/**
 * Helper class for operations on Arrays, such as matrix multiplication and transformation of matrices to and 
 * from strings formatted according to Kattis input/output requirements. 
 */
public class ArrayOperations {
	/**
   * Function that turns a matrix represented by a string into an actual Java
   * 2-D array.
   * @param  matString The string representation of the matrix
   * @return A Java 2-D array representation of the matrix.
   */
  public static Double[][] stringToMatrix(String matString) {
    String[] matStringSplit = matString.split(" ");
    Double[] matStringDouble = new Double[matStringSplit.length];
    for (int i = 0; i < matStringSplit.length; i++) {
      matStringDouble[i] = Double.parseDouble(matStringSplit[i]);
    }

    int numRows = Integer.parseInt(matStringSplit[0]);
    int numCols = Integer.parseInt(matStringSplit[1]);
    Double[][] outputMat = new Double[numRows][numCols];

    int offset = 2; // First two elements of input are matrix dimensions.
    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        outputMat[i][j] = matStringDouble[offset+j];
      }
      offset += numCols; // 'Jump' one row worth of values for next iteration.
    }
    return outputMat;
  }
  
  /**
   * Function that calculates the result of a matrix multiplication between two
   * matrices.
   * @param  mat1 The first of the two matrices to be multiplied.
   * @param  mat2 The second of the two matrices to be multiplied.
   * @return The resulting matrix.
   */
  public static Double[][] matrixMultiply(Double[][] mat1, Double[][] mat2) {
    int numRows = mat1.length;
    /**
     * The number of multiplications needed for every single cell of the
     * resulting matrix is equal to the # of columns in the first matrix (and
     * equally the # of columns in the second matrix).
     */
    int numMultiplications = mat1[0].length;
    int numCols = mat2[0].length;
    Double[][] outputMat = new Double[numRows][numCols];

    for (int i = 0; i < numRows; i++) {
      for (int j = 0; j < numCols; j++) {
        Double columnSum = 0.0;
        for (int k = 0; k < numMultiplications; k++) {
          columnSum += mat1[i][k] * mat2[k][j];
        }
        outputMat[i][j] = columnSum;
      }
    }
    return outputMat;
  }
  
  /**
   * Inverse function of the stringToMatrix function -> turns a 2D-array into the
   * string representation needed to get the solution accepted in Kattis.
   * @param  mat The 2D-array to turn into string representation.
   * @return String representation of the matrix according to the Kattis standards.
   */
  public static String matrixToString(Double[][] mat) {
    String outputString = Integer.toString(mat.length) + " " + Integer.toString(mat[0].length);
    for (int i = 0; i < mat.length; i++) {
      for (int j = 0; j < mat[0].length; j++) {
        outputString += " " + mat[i][j];
      }
    }
    return outputString;
  }
  
  /**
   * Function that tunrs an array represented by a string into an actual array.
   * @param  arrString The string representation of the array to be converted.
   * @return The string converted into an actual integer array.
   */
  public static int[] stringToArray(String arrString) {
    String[] arrStringSplit = arrString.split(" ");
    int arrLen = Integer.parseInt(arrStringSplit[0]);
    int[] outputArr = new int[arrLen];

    for (int i = 0; i < arrLen; i++) {
      outputArr[i] = Integer.parseInt(arrStringSplit[i+1]);
    }
    return outputArr;
  }
  
  /**
   * Transforms an array from actual Array representation to string representation.
   * Inverse function of stringToArray().
   * @param  arr Integer array to be transformed.
   * @return String representation of the integer array. 
   */
  public static String arrayToString(int[] arr) {
    String outputString = "";
    for (int i = 0; i < arr.length; i++) {
      outputString += Integer.toString(arr[i]) + " ";
    }
    return outputString;
  }
}
