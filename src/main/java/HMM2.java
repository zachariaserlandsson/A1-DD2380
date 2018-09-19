package main.java;
import java.util.Scanner;

import main.java.helpers.ArrayOperations;

public class HMM2 {
  /**
   * Function that implements the delta-pass algorithm, used for estimating
   * the sequence of states given an emission sequence.
   * @param  A The A-matrix of an HMM, namely the transition matrix.
   * @param  B The B-matrix of an HMM, namely the emission probability matrix.
   * @param  pi The initial state distribution represented as a 1*N matrix.
   * @param  emissions The observed emission sequence in an array.
   * @return Returns an integer array with the estimation of the state sequence.
   */
  public static int[] deltaPass(Double[][] A, Double[][] B, Double[][] pi, int[] emissions) {
    int numStates = A.length;
    int numEmissions = emissions.length;
    /**
     * Keep two matrices here for simplicity. One is used for holding the actual
     * delta probabilities, while the other one is used for keeping track of the
     * most likely state to have come from at any given step.
     */
    Double[][] deltaMat = new Double[numStates][numEmissions];
    int[][] deltaIndexMat = new int[numStates][numEmissions];

    for (int step = 0; step < numEmissions; step++) {
      for (int i = 0; i < numStates; i++) {
        if (step == 0) {
          deltaMat[i][step] = B[i][emissions[step]] * pi[0][i];
        } else {
          Double maxProb = 0.0; // Value to be stored in the delta matrix.
          int argMax = -1; // Value to be stored in the delta index matrix.
          for (int j = 0; j < numStates; j++) {
            Double currProb = A[j][i] * deltaMat[j][step-1] * B[i][emissions[step]];
            if (currProb >= maxProb) {
              maxProb = currProb;
              argMax = j;
            }
          }
          deltaMat[i][step] = maxProb;
          /**
           * The delta index matrix 'slacks' one step behind since the indices
           * at step k come from the probability distribution at step k-1.
           */
          deltaIndexMat[i][step-1] = argMax;
        }
      }
    }

    /**
     * Two loops for 1) calculating the max probability at the last step (=T) and
     * 2) Setting all the values in the last column of the index matrix to that
     * number (could equally well just store this index in a variable and let the
     * matrix have dimensions N*(T-1) instead).
     */
    Double maxProb = 0.0;
    Double currProb;
    int argMax = -1;
    for (int i = 0; i < numStates; i++) {
      currProb = deltaMat[i][numEmissions-1];
      if (currProb > maxProb) {
        maxProb = currProb;
        argMax = i;
      }
    }
    for (int i = 0; i < numStates; i++) {
      deltaIndexMat[i][numEmissions-1] = argMax;
    }

    /**
     * Fills the probableStates array with values taken from the delta index matrix.
     */
    int[] probableStates = new int[numEmissions];
    for (int i = numEmissions-1; i >= 0; i--){
      if (i == numEmissions - 1) {
        probableStates[i] = deltaIndexMat[0][i];
      } else {
        probableStates[i] = deltaIndexMat[probableStates[i+1]][i];
      }
    }
    return probableStates;
  }

  public static void main(String[] args) {
    Scanner sc = new Scanner(System.in);
    Double[][] A = ArrayOperations.stringToMatrix(sc.nextLine());
    Double[][] B = ArrayOperations.stringToMatrix(sc.nextLine());
    Double[][] pi = ArrayOperations.stringToMatrix(sc.nextLine());
    int[] emissions = ArrayOperations.stringToArray(sc.nextLine());
    sc.close();

    int[] probableStates = deltaPass(A, B, pi, emissions);
    String res = ArrayOperations.arrayToString(probableStates);
    System.out.println(res);
  }
}
