import java.util.Scanner;
import java.util.stream.Stream;

public class HHM0
{
    public static void main(String[] args)
    {
        Scanner in = new Scanner(System.in);

        // Parse data
        // [row][col]
        double[][] transitionMatrix = stringToMatrix(in.nextLine()); // Transition probability matrix
        double[][] emissionMatrix = stringToMatrix(in.nextLine()); // Observation probability matrix
        double[][] initialProbMatrix = stringToMatrix(in.nextLine()); // Initial probability vector

        in.close();

        // Compute x2 -> answer
        double[][] x2 = multiplyMatrices(initialProbMatrix, transitionMatrix);

        double[][] ans = multiplyMatrices(x2, emissionMatrix);

        printAnswer(ans);
    }

    public static double[][] multiplyMatrices(double[][] mA, double[][] mB)
    {
        int colsA = mA[0].length;
        int rowsA = mA.length;

        int colsB = mB[0].length;
        int rowsB = mB.length;

        if (colsA != rowsB) {
            System.out.printf("Cannot multiply a %d x %d matrix with a %d x %d matrix",
                    rowsA, colsA, rowsB, colsB);
            System.exit(1);
        }

        double[][] res = new double[rowsA][colsB];

        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < colsB; j++) {
                double sum = 0;
                for (int k = 0; k < colsA; k++) {
                    sum += mA[i][k] * mB[k][j];
                }
                res[i][j] = sum;
            }
        }
        return res;
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
