//Copyright (C) [2023] [Neha Patel]
import java.util.Scanner;
class MatrixCalculation {
    public double[][] matrix;
    public int rows;
    public int cols;
    public MatrixCalculation(int rows, int cols) {
        this.rows = rows;
        this.cols = cols;
        matrix = new double[rows][cols];
    }
    // Method to input matrix elements
    public void inputMatrix() {
        Scanner scanner = new Scanner(System.in);
        System.out.println("Enter the elements of the matrix :");
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = scanner.nextDouble();
            }
        }
    }
     // Method to display matrix elements
    public void displayMatrix() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }
    public double[][] getMatrix() {
        return matrix;
    }
// Method for matrix addition
    public static MatrixCalculation add(MatrixCalculation A, MatrixCalculation B) {
        if (A.rows != B.rows || A.cols != B.cols) {
            System.out.println("Matrix dimensions are not compatible for addition.");
            return null;
        }
        MatrixCalculation result = new MatrixCalculation(A.rows, A.cols);
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {
                result.matrix[i][j] = A.matrix[i][j] + B.matrix[i][j];
            }
        }
        return result;
    }
// Method for matrix subtraction
    public static MatrixCalculation subtract(MatrixCalculation A, MatrixCalculation B) {
        if (A.rows != B.rows || A.cols != B.cols) {
            System.out.println("Matrix dimensions are not compatible for subtraction.");
            return null;
        }
        MatrixCalculation result = new MatrixCalculation(A.rows, A.cols);
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {
                result.matrix[i][j] = A.matrix[i][j] - B.matrix[i][j];
            }
        }
        return result;
    }
    // Method for matrix multiplication
    public static MatrixCalculation multiply(MatrixCalculation A, MatrixCalculation B) {
        if (A.cols != B.rows) {
            System.out.println("Matrix dimensions are not compatible for multiplication.");
            return null;
        }
        MatrixCalculation result = new MatrixCalculation(A.rows, B.cols);
        result = new MatrixCalculation(A.rows, B.cols);
        result = new MatrixCalculation(A.rows, B.cols);
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < B.cols; j++) {
                result.matrix[i][j] = 0;
                for (int k = 0; k < A.cols; k++) {
                    result.matrix[i][j] += A.matrix[i][k] * B.matrix[k][j];
                }
            }
        }
        return result;
    }
    // Method for matrix inverse
    public static MatrixCalculation inverse(MatrixCalculation A) {
        int n = A.rows;
        if (n != A.cols) {
            System.out.println("Matrix is not square, cannot find the inverse.");
            return null;
        }
        double[][] augmentedMatrix = new double[n][2 * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmentedMatrix[i][j] = A.matrix[i][j];
                augmentedMatrix[i][j + n] = (i == j) ? 1.0 : 0.0;
            }
        }
        for (int i = 0; i < n; i++) {
            double pivot = augmentedMatrix[i][i];
            for (int j = 0; j < 2 * n; j++) {
                augmentedMatrix[i][j] /= pivot;
            }
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = augmentedMatrix[k][i];
                    for (int j = 0; j < 2 * n; j++) {
                        augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                    }
                }
            }
        }
        MatrixCalculation result = new MatrixCalculation(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result.matrix[i][j] = augmentedMatrix[i][j + n];
            }
        }
        return result;
    }
    // Method for matrix adjoint
    public static MatrixCalculation adjoint(MatrixCalculation A) {
        int n = A.rows;
        if (n != A.cols) {
            System.out.println("Matrix is not square, cannot find the adjoint.");
            return null;
        }
        MatrixCalculation result = new MatrixCalculation(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double[][] minorMatrix = getMinorMatrix(A.matrix, i, j);
                result.matrix[i][j] = Math.pow(-1, i + j) * calculateDeterminant(minorMatrix);
            }
        }
        return result;
    }
    // Method for matrix cofactor
    public static MatrixCalculation cofactor(MatrixCalculation A) {
        int n = A.rows;
        if (n != A.cols) {
            System.out.println("Matrix is not square, cannot find the cofactor.");
            return null;
        }
        MatrixCalculation result = new MatrixCalculation(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result.matrix[i][j] = adjoint(A).matrix[i][j];
            }
        }
        return result;
    }

    public static MatrixCalculation transpose(MatrixCalculation A) {
        int rows = A.cols;
        int cols = A.rows;
        MatrixCalculation result = new MatrixCalculation(rows, cols);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.matrix[i][j] = A.matrix[j][i];
            }
        }

        return result;
    }

    public static MatrixCalculation skewSymmetric(MatrixCalculation A) {
        int n = A.rows;
        if (n != A.cols) {
            System.out.println("Matrix is not square, cannot find the skew-symmetric matrix.");
            return null;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (A.matrix[i][j] != -A.matrix[j][i]) {
                    System.out.println("Matrix is not skew-symmetric.");
                    return null;
                }
            }
        }
        return A;
    }
    // Method for matrix rank
    public static int findRank(MatrixCalculation A) {
        double[][] matrix = A.getMatrix();
        int rows = A.rows;
        int cols = A.cols;
        int rank = cols;

        for (int row = 0; row < rank; row++) {
            // Perform Gaussian elimination
            if (matrix[row][row] != 0) {
                for (int i = 0; i < rows; i++) {
                    if (i != row) {
                        double factor = matrix[i][row] / matrix[row][row];
                        for (int j = row; j < cols; j++) {
                            matrix[i][j] -= factor * matrix[row][j];
                        }
                    }
                }
            } else {
                boolean reduce = true;
                for (int i = row + 1; i < rows; i++) {
                    if (matrix[i][row] != 0) {
                        // Swap rows
                        for (int j = row; j < cols; j++) {
                            double temp = matrix[row][j];
                            matrix[row][j] = matrix[i][j];
                            matrix[i][j] = temp;
                        }
                        reduce = false;
                        break;
                    }
                }
                if (reduce) {
                    rank--;
                    for (int i = 0; i < rows; i++) {
                        matrix[i][row] = matrix[i][rank];
                    }
                }
                row--;
            }
        }

        return rank;
    }

    // Method for matrix skewasymmetric
    public static MatrixCalculation skewasymmetric(MatrixCalculation A) {
        int n = A.rows;
        if (n != A.cols) {
            System.out.println("Matrix is not square, cannot find the asymmetric matrix.");
            return null;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (A.matrix[i][j] != -A.matrix[j][i]) {
                    System.out.println("Matrix is not asymmetric.");
                    return null;
                }
            }
        }
        return A;
    }
    public static double[][] getMinorMatrix(double[][] matrix, int row, int col) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double[][] minorMatrix = new double[rows - 1][cols - 1];
        int newRow = 0;
        for (int i = 0; i < rows; i++) {
            if (i == row) continue;
            int newCol = 0;
            for (int j = 0; j < cols; j++) {
                if (j == col) continue;
                minorMatrix[newRow][newCol] = matrix[i][j];
                newCol++;
            }
            newRow++;
        }
        return minorMatrix;
    }
    public static double calculateDeterminant(double[][] matrix) {
        int n = matrix.length;
        if (n == 1) {
            return matrix[0][0];
        } else if (n == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        } else {
            double determinant = 0;
            for (int j = 0; j < n; j++) {
                determinant += matrix[0][j] * Math.pow(-1, j) * calculateDeterminant(getMinorMatrix(matrix, 0, j));
            }
            return determinant;
        }
    }
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.print("Enter the number of rows and columns of matrix A: ");
        int rowsA = scanner.nextInt();
        int colsA = scanner.nextInt();
        System.out.print("Enter the number of rows and columns of matrix B: ");
        int rowsB = scanner.nextInt();
        int colsB = scanner.nextInt();
        MatrixCalculation matrixA = new MatrixCalculation(rowsA, colsA);
        matrixA.inputMatrix();
        MatrixCalculation matrixB = new MatrixCalculation(rowsB, colsB);
        matrixB.inputMatrix();
        System.out.println("Choose an operation:");
        System.out.println("1. Matrix Addition");
        System.out.println("2. Matrix Subtraction");
        System.out.println("3. Matrix Multiplication");
        System.out.println("4. Matrix inverse");
        System.out.println("5. Matrix Adjoint");
        System.out.println("6. Matrix Cofactor");
        System.out.println("7. Matrix SkewSymmetric");
        System.out.println("8. Matrix SkewAsymmetric");
        System.out.println("9. Matrix Transpose");
        System.out.println("10. Matrix Rank");
        // User selects an operation (addition, subtraction, etc.) here...
        System.out.println("Selected choice");
        int s = scanner.nextInt();
        switch (s) {
                // Cases for various matrix operations...
            case 1:
                MatrixCalculation resultAddition = add(matrixA, matrixB);
                if (resultAddition != null) {
                    System.out.println("Matrix A + B:");
                    resultAddition.displayMatrix();
                }
                break;
            case 2:
                MatrixCalculation resultSubtraction = subtract(matrixA, matrixB);
        if (resultSubtraction != null) {
            System.out.println("Matrix A - B:");
            resultSubtraction.displayMatrix();
        }
            break;
            case 3:
                MatrixCalculation resultMultiplication = multiply(matrixA, matrixB);
        if (resultMultiplication != null) {
            System.out.println("Matrix A * B:");
            resultMultiplication.displayMatrix();
        }
        break;
            case 4:
                MatrixCalculation resultInverseA = inverse(matrixA);
        if (resultInverseA != null) {
            System.out.println("Inverse of Matrix A:");
            resultInverseA.displayMatrix();
        }
        break;
            case 5:
                MatrixCalculation resultAdjointA = adjoint(matrixA);
        if (resultAdjointA != null) {
            System.out.println("Adjoint of Matrix A:");
            resultAdjointA.displayMatrix();
        }
          break;
            case 6:
                MatrixCalculation resultCofactorA = cofactor(matrixA);
        if (resultCofactorA != null) {
            System.out.println("Cofactor of Matrix A:");
            resultCofactorA.displayMatrix();
        }
          break;
            case 7:
                MatrixCalculation resultSkewSymmetricA = skewSymmetric(matrixA);
        if (resultSkewSymmetricA != null) {
            System.out.println("Matrix A is skew-symmetric.");
        }
        break;
            case 8:
                MatrixCalculation resultAsymmetricA = skewasymmetric(matrixA);
        if (resultAsymmetricA != null) {
            System.out.println("Matrix A is asymmetric.");
        }
            break;
            case 9:
                MatrixCalculation resultTransposeA = transpose(matrixA);
                if (resultTransposeA != null) {
                    System.out.println("Transpose of Matrix A:");
                    resultTransposeA.displayMatrix();
                }
                break;
            case 10:
                int rankA = findRank(matrixA);
                System.out.println("Rank of Matrix A: " + rankA);
                break;
            default:
            System.out.println("Invalid entry");
            scanner.close();
        }
    }
}
