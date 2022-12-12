import java.util.Scanner;
import java.util.stream.Stream;

public class HMM1
{
    public static void main(String[] args)
    {
        Scanner in = new Scanner(System.in);

        // Parse data
        // [row][col]
        double[][] transitionMatrix = stringToMatrix(in.nextLine()); // Transition probability matrix
        double[][] emissionMatrix = stringToMatrix(in.nextLine()); // Observation probability matrix
        double[][] initialProbMatrix = stringToMatrix(in.nextLine()); // Initial probability vector

        int[] observationSeq = stringToVector(in.nextLine()); // Observation sequence

        in.close();

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

        // Compute the probability of the given sequence as a single scalar
        double answer = 0;
        for (int i = 0; i < Nc; i++) {
            answer += alpha[Tr - 1][i];
        }
        System.out.println(answer);
    }

    public static double[][] stringToMatrix(String str)
    {
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
        str = str.substring(2); // Skip length definition at beginning of line

        return Stream.of(str.split(" ")) // Split line into individual number strings
                .mapToInt(Integer::parseInt) // Parse strings to integers
                .toArray(); // Return as array
    }
}
