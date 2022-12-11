import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Stream;

public class Main
{
    public static void main(String[] args)
    {
        List<double[]> matricesData = parseInputFile(args[0]);

        // [row][col]
        double[][] transitionMatrix = arrayToMatrix(matricesData.get(0)); // Transition probability matrix
        double[][] emissionMatrix = arrayToMatrix(matricesData.get(1)); // Observation probability matrix
        double[][] stateProbDistMatrix = arrayToMatrix(matricesData.get(2)); // Initial probability vector

        printMatrix(transitionMatrix, "Transition matrix");
        printMatrix(emissionMatrix, "Emission matrix");
        printMatrix(stateProbDistMatrix, "Initial state probability distribution matrix");

        int nrOfStates = emissionMatrix.length;
        int nrOfEvents = emissionMatrix[0].length;

        double p01isF1 = 0;
        for (int i = 0; i < stateProbDistMatrix[0].length; i++) {
            p01isF1 += emissionMatrix[i][0] * stateProbDistMatrix[0][i];
        }
        System.out.println(p01isF1);

        double[][] x2 = multiplyMatrices(stateProbDistMatrix, transitionMatrix);
        printMatrix(x2, "X2");

        double pO2isF2 = 0;
        for (int i = 0; i < x2[0].length; i++) {
            pO2isF2 += emissionMatrix[i][1] * x2[0][i];
        }
        System.out.println(pO2isF2);

        double[][] x3 = multiplyMatrices(x2, transitionMatrix);

        printMatrix(x3, "X3");
    }

    public static List<double[]> parseInputFile(String file)
    {
        try {
            List<String> matrices = Files.readAllLines(Path.of(file));
            return matrices.stream()
                    .map(str -> Stream.of(str.split(" "))
                            .mapToDouble(Double::parseDouble)
                            .toArray())
                    .toList();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static double[][] arrayToMatrix(double[] array)
    {
        int rows = (int) array[0];
        int cols = (int) array[1];

        double[][] matrix = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = array[(i * cols) + j + 2]; // +2 to skip row/col definitions
            }
        }

        return matrix;
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

    public static void printMatrix(double[][] matrix, String name)
    {
        StringBuilder sb = new StringBuilder(name);
        sb.append('\n');
        sb.append(matrix.length); // Rows
        sb.append(" x ");
        sb.append(matrix[0].length); // Cols
        sb.append(" matrix\n");

        for (double[] row : matrix) {
            sb.append("| ");
            for (double col : row) {
                sb.append(String.format("%.2f", col));
                sb.append(' ');
            }
            sb.append("|\n");
        }
        System.out.println(sb);
    }
}
