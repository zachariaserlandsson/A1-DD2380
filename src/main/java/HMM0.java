package main.java;
import java.util.Scanner;

import main.java.helpers.ArrayOperations;

public class HMM0 {
  public static void main(String[] args) {
    Scanner sc = new Scanner(System.in);
    Double[][] A = ArrayOperations.stringToMatrix(sc.nextLine());
    Double[][] B = ArrayOperations.stringToMatrix(sc.nextLine());
    Double[][] pi = ArrayOperations.stringToMatrix(sc.nextLine());
    sc.close();

    Double[][] transitionDistribution = ArrayOperations.matrixMultiply(pi, A);
    Double[][] emissionDistribution = ArrayOperations.matrixMultiply(transitionDistribution, B);
    String emissionString = ArrayOperations.matrixToString(emissionDistribution);
    System.out.println(emissionString);
  }
}
