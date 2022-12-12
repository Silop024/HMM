import java.util.Scanner;
import java.util.stream.Stream;

public class HMM3
{
    public static void main(String[] args)
    {
        Scanner in = new Scanner(System.in);

        // Parse data (1. Initialize ¸ lambda=(A,B,pi))
        // [row][col]
        double[][] transitionMatrix = stringToMatrix(in.nextLine()); // Transition probability matrix
        double[][] emissionMatrix = stringToMatrix(in.nextLine()); // Observation probability matrix
        double[][] initialProbMatrix = stringToMatrix(in.nextLine()); // Initial probability vector

        int[] emissionSequence = stringToVector(in.nextLine()); // Observation sequence

        in.close();

        /*
         *   2. Compute alpha, beta, digamma, gamma values
         *   3. Re-estimate ¸ lambda=(A,B,pi)
         *   4. Repeat from 2. until convergence
         */
        for (int loop = 0; loop < 2; loop++) {
            baumWelch(transitionMatrix, emissionMatrix, initialProbMatrix, emissionSequence);
        }
        printAnswer(transitionMatrix);
        printAnswer(emissionMatrix);
    }

    public static double[][] alphaPass(double[][] transitionMatrix, double[][] emissionMatrix,
                                       double[][] initialProbMatrix, int[] observationSeq)
    {
        // Initialize alpha
        int Tr = observationSeq.length;
        int Nc = transitionMatrix.length;

        double[][] alpha = new double[Tr][Nc];

        for (int i = 0; i < Nc; i++) {
            alpha[0][i] = initialProbMatrix[0][i] * emissionMatrix[i][observationSeq[0]];
        }

        // Compute rest of alpha
        for (int t = 1; t < Tr; t++) {
            for (int i = 0; i < Nc; i++) {
                double sum = 0;

                for (int j = 0; j < Nc; j++) {
                    sum += alpha[t - 1][j] * transitionMatrix[j][i];
                }
                alpha[t][i] = sum * emissionMatrix[i][observationSeq[t]];
            }
        }
        return alpha;
    }

    public static double[][] betaPass(double[][] transitionMatrix, double[][] emissionMatrix, int[] emissionSequence)
    {
        // Initialize beta
        //int Tr = emissionSequence.length;
        int Tr = transitionMatrix.length;
        int Nc = transitionMatrix.length;

        double[][] beta = new double[Tr][Nc];

        for (int i = 0; i < Nc; i++) {
            beta[Tr - 1][i] = 1;
        }

        // Compute rest of beta
        for (int t = Tr - 2; t >= 0; t--) {
            for (int i = 0; i < Nc; i++) {
                double sum = 0;

                for (int j = 0; j < Nc; j++) {
                    sum += beta[t + 1][j] * emissionMatrix[i][emissionSequence[t]] * transitionMatrix[j][i];
                }
                beta[t][i] = sum;
            }
        }

        return beta;
    }

    public static double[][][] diGamma(double[][] alpha, double[][] beta,
                                       double[][] transitionMatrix, double[][] emissionMatrix,
                                       int[] emissionSequence)
    {
        //int Tr = emissionSequence.length;
        int Tr = transitionMatrix.length;
        int Nc = transitionMatrix.length;

        double[][][] gamma = new double[Tr][Nc][Nc];

        for (int t = 0; t < Tr - 1; t++) {
            for (int i = 0; i < Nc; i++) {
                for (int j = 0; j < Nc; j++) {
                    double numerator = alpha[t][i] * transitionMatrix[i][j] * emissionMatrix[j][emissionSequence[t + 1]] * beta[t + 1][j];

                    double denominator = 0;

                    for (int k = 0; k < Nc; k++) {
                        denominator += alpha[Tr - 1][k];
                    }
                    gamma[t][i][j] = numerator / denominator;
                }
            }
        }
        return gamma;
    }

    public static double[][] gamma(double[][][] diGamma)
    {
        int Tr = diGamma.length;
        int N = diGamma[0].length;

        double[][] gamma = new double[Tr][N];

        for (int t = 0; t < Tr; t++) {
            for (int i = 0; i < N; i++) {
                double sum = 0;

                for (int j = 0; j < N; j++) {
                    sum += diGamma[t][i][j];
                }
                gamma[t][i] = sum;
            }
        }
        return gamma;
    }

    public static void baumWelch(double[][] transitionMatrix, double[][] emissionMatrix, double[][] initialProbMatrix, int[] emissionSequence)
    {
        double[][] alpha = alphaPass(transitionMatrix, emissionMatrix, initialProbMatrix, emissionSequence);
        double[][] beta = betaPass(transitionMatrix, emissionMatrix, emissionSequence);
        double[][][] diGamma = diGamma(alpha, beta, transitionMatrix, emissionMatrix, emissionSequence);
        double[][] gamma = gamma(diGamma);

        int N = transitionMatrix.length;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double sumDiGamma = 0;
                double sumGamma = 0;

                for (int t = 0; t < gamma.length; t++) {
                    sumDiGamma += diGamma[t][i][j];
                    sumGamma += gamma[t][i];
                }
                transitionMatrix[i][j] = sumDiGamma / sumGamma;
            }
        }

        int K = emissionMatrix[0].length;
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < K; k++) {
                double nominator = 0;
                double denominator = 0;

                for (int t = 0; t < gamma.length; t++) {
                    int idk = 0;
                    if (emissionSequence[t] == k) {
                        idk = 1;
                    }
                    nominator += idk * gamma[t][j];

                    denominator += gamma[t][j];
                }
                emissionMatrix[j][k] = nominator / denominator;
            }
        }

        System.arraycopy(gamma[0], 0, initialProbMatrix[0], 0, N);
    }

    public static double[][] stringToMatrix(String str)
    {
        str = str.trim(); // Remove empty spaces at beginning and end of line

        double[] array = Stream.of(str.split(" ")) // Split line into individual number strings
                .mapToDouble(Double::parseDouble) // Parse strings to doubles
                .toArray(); // Save as array

        // Nr of rows defined as first number in string, column as second
        int rows = (int) array[0];
        int cols = (int) array[1];

        double[][] matrix = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = array[(i * cols) + j + 2]; // +2 to skip row/col definitions at beginning of line
            }
        }
        return matrix;
    }

    public static int[] stringToVector(String str)
    {
        return Stream.of(str.split(" ")) // Split line into individual number strings
                .skip(1) // Skip length definition at beginning of line
                .mapToInt(Integer::parseInt) // Parse strings to integers
                .toArray(); // Return as array
    }

    public static void printAnswer(double[][] ans)
    {
        int rows = ans.length;
        int cols = ans[0].length;

        StringBuilder sb = new StringBuilder();
        sb.append(rows);
        sb.append(' ');
        sb.append(cols);
        sb.append(' ');

        for (double[] row : ans) {
            for (double col : row) {
                sb.append(col);
                sb.append(' ');
            }
        }
        System.out.println(sb);
    }
}
