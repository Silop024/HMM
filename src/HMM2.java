import java.util.Scanner;
import java.util.stream.Stream;

public class HMM2
{
    public static void main(String[] args)
    {
        Scanner in = new Scanner(System.in);

        // Parse data
        // [row][col]
        double[][] transitionMatrix = stringToMatrix(in.nextLine()); // Transition probability matrix
        double[][] emissionMatrix = stringToMatrix(in.nextLine()); // Observation probability matrix
        double[][] initialProbMatrix = stringToMatrix(in.nextLine()); // Initial probability vector

        int[] emissionSequence = stringToVector(in.nextLine()); // Observation sequence

        in.close();

        // Initialize delta
        int Tr = emissionSequence.length;
        int Nc = transitionMatrix.length;

        double[][] delta = new double[Tr][Nc];
        int[][] deltaIndices = new int[Tr][Nc];

        for (int i = 0; i < Nc; i++) {
            delta[0][i] = initialProbMatrix[0][i] * emissionMatrix[i][emissionSequence[0]];
        }

        // Compute rest of delta
        for (int t = 1; t < Tr; t++) {
            for (int i = 0; i < Nc; i++) {
                double max = -1;
                int max_index = -1;

                for (int j = 0; j < Nc; j++) {
                    double tmp = delta[t - 1][j] * transitionMatrix[j][i] * emissionMatrix[i][emissionSequence[t]];

                    if (tmp > max) {
                        max = tmp;
                        max_index = j;
                    }
                }
                delta[t][i] = max;
                deltaIndices[t][i] = max_index;
            }
        }

        // Compute answer
        int[] answer = new int[Tr];

        double max = -1;
        for (int j = 0; j < Nc; j++) {
            double curr = delta[Tr - 1][j];

            if (curr > max) {
                max = curr;
                answer[Tr - 1] = j;
            }
        }

        for (int t = Tr - 2; t >= 0; t--) {
            answer[t] = deltaIndices[t + 1][answer[t + 1]];
        }

        // Print answer
        StringBuilder sb = new StringBuilder();
        for (int val : answer) {
            sb.append(val);
            sb.append(' ');
        }
        System.out.println(sb);
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
        return Stream.of(str.split(" ")) // Split line into individual number strings
                .skip(1) // Skip length definition at beginning of line
                .mapToInt(Integer::parseInt) // Parse strings to integers
                .toArray(); // Return as array
    }
}
