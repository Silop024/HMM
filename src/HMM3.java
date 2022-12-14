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
        baumWelch(transitionMatrix, emissionMatrix, initialProbMatrix, emissionSequence);
    }

    public static void alphaPass(double[][] alpha, double[][] transitionMatrix, double[][] emissionMatrix,
                                 double[][] initialProbMatrix, int[] emissionSequence,
                                 double[] c, int Tr, int Nc)
    {
        // Compute alpha[0][i]
        c[0] = 0;
        for (int i = 0; i < Nc; i++) {
            alpha[0][i] = initialProbMatrix[0][i] * emissionMatrix[i][emissionSequence[0]];
            c[0] += alpha[0][i];
        }

        // Scale the alpha[0][i]
        c[0] = 1 / c[0];
        for (int i = 0; i < Nc; i++) {
            alpha[0][i] *= c[0];
        }

        // Compute a[t][i]
        for (int t = 1; t < Tr; t++) {
            c[t] = 0;
            for (int i = 0; i < Nc; i++) {
                alpha[t][i] = 0;
                for (int j = 0; j < Nc; j++) {
                    alpha[t][i] += alpha[t - 1][j] * transitionMatrix[j][i];
                }
                alpha[t][i] *= emissionMatrix[i][emissionSequence[t]];
                c[t] += alpha[t][i];
            }
            // Scale alpha[t][i]
            c[t] = 1 / c[t];
            for (int i = 0; i < Nc; i++) {
                alpha[t][i] *= c[t];
            }
        }
    }

    public static void betaPass(double[][] beta, double[][] transitionMatrix, double[][] emissionMatrix,
                                int[] emissionSequence, double[] c,
                                int Tr, int Nc)
    {
        // Let beta[Tr - 1][i] = 1, scaled by c[Tr - 1]
        for (int i = 0; i < Nc; i++) {
            beta[Tr - 1][i] = c[Tr - 1];
        }

        // Beta-pass
        for (int t = Tr - 2; t >= 0; t--) {
            for (int i = 0; i < Nc; i++) {
                beta[t][i] = 0;
                for (int j = 0; j < Nc; j++) {
                    beta[t][i] += transitionMatrix[i][j] * emissionMatrix[j][emissionSequence[t + 1]] * beta[t + 1][j];
                }
                beta[t][i] *= c[t];
            }
        }
    }

    public static void diGamma(double[][] alpha, double[][] beta, double[][] gamma, double[][][] diGamma,
                               double[][] transitionMatrix, double[][] emissionMatrix,
                               int[] emissionSequence, int Tr, int Nc)
    {

        // No need to normalize gamma[t][i][j] since using scaled alpha and beta
        for (int t = 0; t < Tr - 1; t++) {
            for (int i = 0; i < Nc; i++) {
                gamma[t][i] = 0;
                for (int j = 0; j < Nc; j++) {
                    diGamma[t][i][j] = alpha[t][i] * transitionMatrix[i][j] * emissionMatrix[j][emissionSequence[t + 1]] * beta[t + 1][j];

                    gamma[t][i] += diGamma[t][i][j];
                }
            }
        }
        // Special case for gamma[Tr - 1][i] (as above, no need to normalize)
        System.arraycopy(alpha[Tr - 1], 0, gamma[Tr - 1], 0, Nc);
    }

    public static void baumWelch(double[][] transitionMatrix, double[][] emissionMatrix,
                                 double[][] initialProbMatrix, int[] emissionSequence)
    {
        int iters = 0;
        int maxIters = 100;

        double oldLogProb = Double.NEGATIVE_INFINITY;

        int Tr = emissionSequence.length;
        int Nc = transitionMatrix.length;

        double[] c = new double[emissionSequence.length];
        double[][] alpha = new double[Tr][Nc];
        double[][] beta = new double[Tr][Nc];
        double[][] gamma = new double[Tr][Nc];
        double[][][] diGamma = new double[Tr][Nc][Nc];

        while (true) {
            alphaPass(alpha, transitionMatrix, emissionMatrix, initialProbMatrix, emissionSequence, c, Tr, Nc);
            betaPass(beta, transitionMatrix, emissionMatrix, emissionSequence, c, Tr, Nc);
            diGamma(alpha, beta, gamma, diGamma, transitionMatrix, emissionMatrix, emissionSequence, Tr, Nc);

            // Re-estimate probability matrix
            System.arraycopy(gamma[0], 0, initialProbMatrix[0], 0, Nc);

            double denominator;
            double numerator;
            // Re-estimate transition matrix
            for (int i = 0; i < Nc; i++) {
                denominator = 0;
                for (int t = 0; t < Tr - 1; t++) {
                    denominator += gamma[t][i];
                }
                for (int j = 0; j < Nc; j++) {
                    numerator = 0;

                    for (int t = 0; t < Tr - 1; t++) {
                        numerator += diGamma[t][i][j];
                    }
                    transitionMatrix[i][j] = numerator / denominator;
                }
            }

            // Re-estimate emission matrix
            for (int i = 0; i < Nc; i++) {
                denominator = 0;
                for (int t = 0; t < Tr; t++) {
                    denominator += gamma[t][i];
                }
                for (int j = 0; j < emissionMatrix[0].length; j++) {
                    numerator = 0;

                    for (int t = 0; t < Tr; t++) {
                        if (emissionSequence[t] == j) {
                            numerator += gamma[t][i];
                        }
                    }
                    emissionMatrix[i][j] = numerator / denominator;
                }
            }

            // Compute log[P(O|lambda]]
            double logProb = 0;
            for (int i = 0; i < Tr; i++) {
                logProb += Math.log(c[i]);
            }
            logProb = -logProb;

            iters++;

            if (iters < maxIters && logProb > oldLogProb) {
                oldLogProb = logProb;
            } else {
                printAnswer(transitionMatrix);
                printAnswer(emissionMatrix);
                System.exit(0);
            }
        }
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
