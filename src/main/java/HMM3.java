package main.java;
import java.util.Scanner;
import main.java.helpers.AlphaPass;
import main.java.helpers.Gamma;
import main.java.helpers.ArrayOperations;
import java.lang.Math;

public class HMM3 {
  /**
   * Function that calculates the probability of an emission sequence using the
   * forward pass (alpha-pass) algorithm, as well as scaling the alpha matrix.
   * @param  A The A-matrix of an HMM, namely the transition matrix.
   * @param  B The B-matrix of an HMM, namely the emission probability matrix.
   * @param  pi The initial state distribution represented as a 1*N matrix.
   * @param  emissions The observed emission sequence in an array.
   * @return Returns an AlphaPass object containing the Alpha matrix as well as
   * the scaling factors.
   */
  public static AlphaPass alphaPassScale(Double[][] A, Double[][] B, Double[][] pi, int[] emissions) {
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
    Double[] scalingFactors = new Double[numEmissions];
    Double c; // Variable used for scaling.
    for (int step = 0; step < numEmissions; step++) {
      c = 0.0;
      for (int i = 0; i < numStates; i++) {
        /**
         * If we're at the first step we initialize the alpha matrix with the
         * rules defined, i.e. making use of the initial state distribution.
         */
        if (step == 0) {
          alphaMat[i][step] = B[i][emissions[step]] * pi[0][i];
          c += alphaMat[i][step];
        } else { // Else we use the normal method of using the previous time step's values
          Double transitionSum = 0.0;
          for (int j = 0; j < numStates; j++) {
            transitionSum += A[j][i] * alphaMat[j][step-1];
          }
          alphaMat[i][step] = B[i][emissions[step]] * transitionSum;
          c += alphaMat[i][step];
        }
        /**
         * If last value in current column, we scale the entirety of the current
         * column of the matrix, as per the Stamp tutorial.
         */
        if (i == numStates-1) {
          c = 1/c;
          for (int scaleI = 0; scaleI < numStates; scaleI++) {
            alphaMat[scaleI][step] *= c;
          }
        }
      }
      scalingFactors[step] = c;
    }
    return new AlphaPass(alphaMat, scalingFactors);
  }

  /**
   * Function that performs the Beta-pass (or backwards-pass) algorithm for computing
   * the most likely state sequence, starting from the back and progressing to the front
   * (time-step-wise).
   * @param  A The A-matrix of an HMM, namely the transition matrix.
   * @param  B The B-matrix of an HMM, namely the emission probability matrix.
   * @param  pi The initial state distribution represented as a 1*N matrix.
   * @param  emissions The observed emission sequence in an array.
   * @return Returns a 2D array containing the values of the Beta matrix.
   */
  public static Double[][] betaPassScale(Double[][] A, Double[][] B, Double[][] pi,
                                         int[] emissions, Double[] scalingFactors) {
    int numStates = A.length;
    int numEmissions = emissions.length;
    Double[][] betaMat = new Double[numStates][numEmissions];
    /**
     * Setting the last column of the beta matrix to the scaling factors.
     */
    for (int i = 0; i < numStates; i++) {
      betaMat[i][numEmissions-1] = scalingFactors[numEmissions-1];
    }

    Double currProb;
    /**
     * Starting from the back instead of the front since this is a backwards-pass.
     */
    for (int step = numEmissions-2; step >= 0; step--) {
      for (int i = 0; i < numStates; i++) {
        currProb = 0.0;
        for (int j = 0; j < numStates; j++) {
          currProb += A[i][j] * B[j][emissions[step+1]] * betaMat[j][step+1];
        }
        currProb *= scalingFactors[step];
        betaMat[i][step] = currProb;
      }
    }
    return betaMat;
  }

  /**
   * Function for computing the two Gamma matrices, here called mono-gamma (the
   * probability of being in state i at time t given the emission sequence and the
   * current estimation of the HMM model) and di-gamma (the probability of being in
   * state i at time t and transitioning to state j at the next time step, given the
   * same parameters).
   * @param  A The A-matrix of an HMM, namely the transition matrix.
   * @param  B The B-matrix of an HMM, namely the emission probability matrix.
   * @param  pi The initial state distribution represented as a 1*N matrix.
   * @param  emissions The observed emission sequence in an array.
   * @param  alphaObj An object containing the alpha matrix and the scaling factors.
   * @param  beta A 2D-array containing the beta values.
   * @return Returns a Gamma object containing the 2D-array monogamma and the 3D-array
   * digamma (extra dimension since we're considering time as an extra dimension).
   */
  public static Gamma computeGamma(Double[][] A, Double[][] B, Double[][] pi,
                                        int[] emissions, AlphaPass alphaObj,
                                        Double[][] beta) {
    int numStates = A.length;
    int numEmissions = emissions.length;
    Double[][] alpha = alphaObj.alphaMat;

    Double[][] monoGamma = new Double[numStates][numEmissions];
    Double[][][] diGamma = new Double[numStates][numStates][numEmissions];
    Double denom;
    Double currVal;
    for (int step = 0; step < numEmissions-1; step++) {
      denom = 0.0;
      for (int i = 0; i < numStates; i++) {
        for (int j = 0; j < numStates; j++) {
          denom += alpha[i][step] * A[i][j] * B[j][emissions[step+1]] * beta[j][step+1];
        }
      }
      for (int i = 0; i < numStates; i++) {
        currVal = 0.0;
        for (int j= 0; j < numStates; j++) {
          diGamma[i][j][step] = (alpha[i][step] * A[i][j] * B[j][emissions[step+1]]
                                 * beta[j][step+1]) / denom;
          currVal += diGamma[i][j][step];
        }
        monoGamma[i][step] = currVal;
      }
    }

    denom = 0.0;
    for (int i = 0; i < numStates; i++) {
      denom += alpha[i][numEmissions-1];
    }
    for (int i = 0; i < numStates; i++) {
      monoGamma[i][numEmissions-1] = alpha[i][numEmissions-1] / denom;
    }
    return new Gamma(monoGamma, diGamma);
  }

  /**
   * The main function for training the HMM model. Calls all the pass algorithms as well
   * as the gamma, re-estimation and logarithm probability functions. Void function
   * so only changes the values of the matrices instead of returning anything.
   * @param aApprox The approximated A-matrix of an HMM, namely the transition matrix.
   * @param bApprox The approximated B-matrix of an HMM, namely the emission probability matrix.
   * @param piApprox The estimated initial state distribution represented as a 1*N matrix.
   * @param emissions The observed emission sequence in an array.
   */
  public static void trainModel(Double[][] aApprox, Double[][] bApprox, Double[][] piApprox, int[] emissions) {
    /**
     * Value set to fit the Kattis time limits, in case the log probability doesn't
     * converge quickly enough.
     */
    int maxIters = 100;
    Double oldLogProb = (-1) * Double.MAX_VALUE;

    Double logProb;
    /**
     * Runs for 1 to maxIters iterations, breaking early if the logarithm probability
     * converges.
     */
    for (int iter = 0; iter < maxIters; iter++) {
      AlphaPass alpha = alphaPassScale(aApprox, bApprox, piApprox, emissions);
      Double[][] beta = betaPassScale(aApprox, bApprox, piApprox, emissions, alpha.scalingFactors);
      Gamma gamma = computeGamma(aApprox, bApprox, piApprox, emissions, alpha, beta);
      reEstimate(aApprox, bApprox, piApprox, gamma, emissions);
      logProb = calculateLogProb(alpha.scalingFactors);
      if (logProb > oldLogProb) {
        oldLogProb = logProb;
      } else {
        break;
      }
    }

  }

  /**
   * Function for re-estimating the A, B and pi matrices using the gamma matrices.
   * @param aApprox The approximated A-matrix of an HMM, namely the transition matrix.
   * @param bApprox The approximated B-matrix of an HMM, namely the emission probability matrix.
   * @param piApprox The estimated initial state distribution represented as a 1*N matrix.
   * @param gammaObj Gamma object containing the 2D-array monogamma and the 3D-array
   * digamma (extra dimension since we're considering time as an extra dimension).
   * @param emissions The observed emission sequence in an array.
   */
  public static void reEstimate(Double[][] aApprox, Double[][] bApprox, Double[][] piApprox,
                                Gamma gammaObj, int[] emissions) {
    int numStates = aApprox.length;
    int numEmissions = emissions.length;
    Double[][] monoGamma = gammaObj.monoGamma;
    Double[][][] diGamma = gammaObj.diGamma;

    /**
     * Re-estimates the values of the initial probability distribution using the first
     * time-step of the mono-gamma array.
     */
    for (int i = 0; i < numStates; i++) {
      piApprox[0][i] = monoGamma[i][0];
    }

    Double numer;
    Double denom;
    /**
     * Re-estimates the A matrix using both the di-gamma and the mono-gamma arrays.
     */
    for (int i = 0; i < numStates; i++) {
      for (int j = 0; j < numStates; j++) {
        numer = 0.0;
        denom = 0.0;
        for (int step = 0; step < numEmissions-1; step++) {
          numer += diGamma[i][j][step];
          denom += monoGamma[i][step];
        }
        aApprox[i][j] = numer/denom;
      }
    }

    /**
    * Re-estimates the B matrix using the mono-gamma arrays.
     */
    for (int i = 0; i < numStates; i++) {
      for (int j = 0; j < bApprox[0].length; j++) {
        numer = 0.0;
        denom = 0.0;
        for (int step = 0; step < numEmissions; step++) {
          if (emissions[step] == j) {
            numer += monoGamma[i][step];
          }
          denom += monoGamma[i][step];
        }
        bApprox[i][j] = numer/denom;
      }
    }
  }

  /**
   * Function for calculating the natural logarithm probability of having observed the
   * given emission sequence given the current HMM estimation.
   * @param  scalingFactors Array holding the scaling factors (c in the stamp tutorial).
   * @return Returns the logarithm probability as a Double value.
   */
  public static Double calculateLogProb(Double[] scalingFactors) {
    Double logProb = 0.0;
    for (int i = 0; i < scalingFactors.length; i++) {
      logProb += Math.log(scalingFactors[i]); // TODO: Natural logarithm here?
    }
    return logProb * (-1);
  }

  public static void main(String[] args) {
    Scanner sc = new Scanner(System.in);
    Double[][] A = ArrayOperations.stringToMatrix(sc.nextLine());
    Double[][] B = ArrayOperations.stringToMatrix(sc.nextLine());
    Double[][] pi = ArrayOperations.stringToMatrix(sc.nextLine());
    int[] emissions = ArrayOperations.stringToArray(sc.nextLine());
    sc.close();

    trainModel(A, B, pi, emissions);
    System.out.println(ArrayOperations.matrixToString(A));
    System.out.println(ArrayOperations.matrixToString(B));
  }
}
