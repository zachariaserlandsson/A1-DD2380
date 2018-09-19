package main.java;
import java.util.Scanner;

import main.java.helpers.ArrayOperations;

public class HMM1 {
  /**
   * Function that calculates the probability of an emission sequence using the
   * forward pass (alpha-pass) algorithm.
   * @param  A The A-matrix of an HMM, namely the transition matrix.
   * @param  B The B-matrix of an HMM, namely the emission probability matrix.
   * @param  pi The initial state distribution represented as a 1*N matrix.
   * @param  emissions The observed emission sequence in an array.
   * @return Returns the probability of the observed sequence as a scalar value.
   */
  public static Double alphaPass(Double[][] A, Double[][] B, Double[][] pi, int[] emissions) {
    int numStates = A.length;
    int numEmissions = emissions.length;
    /**
     * Matrix holding the alpha values. Since we need to keep track of both the
     * current values and the values at the time-step before, an array will not
     * suffice. The matrix has t=|emissions| columns to keep hold of the alpha
     * values at every step. To free up RAM we could have a N*2 matrix instead
     * where we just shift the column at every time step.
     */
    Double[][] alphaMat = new Double[numStates][numEmissions];

    for (int step = 0; step < numEmissions; step++) {
      for (int i = 0; i < numStates; i++) {
        /**
         * If we're at the first step we initialize the alpha matrix with the
         * rules defined, i.e. making use of the initial state distribution.
         */
        if (step == 0) {
          alphaMat[i][step] = B[i][emissions[step]] * pi[0][i];
        } else { // Else we use the normal method of using the previous time step's values
          Double transitionSum = 0.0;
          for (int j = 0; j < numStates; j++) {
            transitionSum += A[j][i] * alphaMat[j][step-1];
          }
          alphaMat[i][step] = B[i][emissions[step]] * transitionSum;
        }
      }
    }

    /**
     * What we're outputting is the scalar value of the emission probability,
     * so we sum the values of the last alpha-vector.
     */
    double probSum = 0.0;
    for (int i = 0; i < numStates; i++) {
      probSum += alphaMat[i][numEmissions-1];
    }
    return probSum;
  }

  public static void main(String[] args) {
    Scanner sc = new Scanner(System.in);
    Double[][] A = ArrayOperations.stringToMatrix(sc.nextLine());
    Double[][] B = ArrayOperations.stringToMatrix(sc.nextLine());
    Double[][] pi = ArrayOperations.stringToMatrix(sc.nextLine());
    int[] emissions = ArrayOperations.stringToArray(sc.nextLine());
    sc.close();
    
    Double res = alphaPass(A, B, pi, emissions);
    System.out.println(res);
  }
}
