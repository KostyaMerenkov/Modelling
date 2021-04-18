package students.dvfu.ru;

import java.io.*;
import java.util.Locale;
import java.util.Random;

public class Main {
    private static final int size = 19;
    private static final int N = size * size;
    private static final double error = 5;
    final static Random random = new Random();

    private static int[][] createAndFillArray() {
        int[][] arr = new int[size][size];
        int[] value = {1, -1};
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                arr[i][j] = value[(int) (Math.random() * value.length)];
            }
        }
        return arr;
    }

    private static void printArray(int[][] arr) {
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                System.out.print(arr[i][j] + "\t");
            }
            System.out.println();
        }
        System.out.println();
    }

    private static int[] getNeighbours(int[][] arr, int x, int y) {
        int[] neighbour = new int[4]; //up, right, down, left
        if (x == 0) { //up
            neighbour[0] = arr[arr.length - 1][y];
        } else {
            neighbour[0] = arr[x - 1][y];
        }
        if (y == arr.length - 1) { //right
            neighbour[1] = arr[x][0];
        } else {
            neighbour[1] = arr[x][y + 1];
        }
        if (x == arr.length - 1) { //down
            neighbour[2] = arr[0][y];
        } else {
            neighbour[2] = arr[x + 1][y];
        }
        if (y == 0) { //left
            neighbour[3] = arr[x][arr.length - 1];
        } else {
            neighbour[3] = arr[x][y - 1];
        }
        return neighbour;
    }

    private static void printNeighbours(int element, int[] neighbour) {
        System.out.println(" \t" + neighbour[0]);
        System.out.println(neighbour[3] + "\t" + element + "\t" + neighbour[1]);
        System.out.println(" \t" + neighbour[2]);
    }

    private static double calcFullEnergy(int[][] arr) {
        float E_FULL = 0;
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                E_FULL += calcEnergy(arr, i, j);
            }
        }
        return E_FULL / N;
    }

    private static double calcEnergy(int[][] arr, int x, int y) {
        int[] neighbour = getNeighbours(arr, x, y);
        float E = 0;
        for (int i = 0; i < neighbour.length; i++) {
            E += (arr[x][y] * neighbour[i]);
        }
        E *= -1;
        return E;
    }

    private static double calcMagnetization(int[][] arr) {
        float M = 0;
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                M += arr[i][j];
            }
        }
        return Math.abs(M / N);
    }

    private static double calcThermal(double E_average, double E_average_squared, double T) {
        E_average = E_average / 5;
        E_average_squared = E_average_squared / 5;
        return (E_average_squared - (E_average * E_average)) / (T * T);
    }

    private static double calcError(double[] xi, double x) {
        x/=error;
        double sum = 0;
        for (int i = 0; i < error; i++) {
            sum += (xi[i] - x) * (xi[i] - x);
        }
        return Math.sqrt(1/(error-1) * sum);
    }

    public static void main(String[] args) throws IOException {
        int[][] arr = createAndFillArray();
        //printArray(arr);
        System.out.println();
//        FileWriter writer = new FileWriter(new File("Results.txt"));
//        writer.append("T\t\tM\t\tm\t\tE\t\te\t\tC");

        FileWriter writer = new FileWriter(new File("Results.txt"));
        writer.append("T\t\tE\t\te\t\tM\t\tm");

        FileWriter writer_TE = new FileWriter(new File("TE.txt"));
        writer_TE.append("T\t\t E");

        FileWriter writer_TM = new FileWriter(new File("TM.txt"));
        writer_TM.append("T\t\tM");

        FileWriter writer_TC = new FileWriter(new File("TC.txt"));
        writer_TM.append("T\t\tC");

        for (double T = 0.10001; T <= 5.00001; T += 0.1) {
            double E_average = 0;
            double E_average_squared = 0;
            double[] E_av = new double[(int)error];
            double[] M_av = new double[(int)error];
            double E_full = 0;
            double M_average = 0;
            double C;
            FileWriter massive = new FileWriter(
                    new File("data/" + String.format(Locale.ROOT, "%.5f", T) + ".txt"));
            for (int average = 0; average < error; average++) {
                for (int MK = 0; MK < 1000000; MK++) {
                    if (average == 0 && MK % 10000 == 0) {
                        for (int i = 0; i < arr.length; i++) {
                            for (int j = 0; j < arr.length; j++) {
                                massive.append(arr[j][i] + " ");
                            }
                        }
                        massive.append("\n");
                        massive.flush();
                    }
                    int r_i = random.nextInt(size), r_j = random.nextInt(size);
                    double E1 = calcEnergy(arr, r_i, r_j);
                    arr[r_i][r_j] = (arr[r_i][r_j]) * -1;
                    double E2 = calcEnergy(arr, r_i, r_j);
                    if (E2 >= E1) {
                        double randExp = Math.exp((-1 * (E2 - E1)) / T);
                        if (random.nextDouble() > randExp) {
                            arr[r_i][r_j] = (arr[r_i][r_j]) * -1;
                        }
                    }
                }
                //printArray(arr);
                E_full = calcFullEnergy(arr);
                E_av[average] = E_full;
                E_average += E_full;
                E_average_squared += (E_full * E_full);
                M_av[average] = calcMagnetization(arr);
                M_average += M_av[average];
                massive.close();
            }
            double E_error = calcError(E_av, E_average);
            double M_error = calcError(M_av, M_average);
            M_average /= error;
            C = calcThermal(E_average, E_average_squared, T);
//            writer.append("\n").append(String.format(Locale.ROOT, "%1.1f ", T)).append("\t")
//                    .append(String.format(Locale.ROOT, "%1.1f", M_average)).append("\t   ")
//                    .append((String.format(Locale.ROOT, "%1.2f", M_error))).append("\t   ")
//                    .append(String.format(Locale.ROOT, "%1.1f", E_full)).append("\t    ")
//                    .append((String.format(Locale.ROOT, "%2.2f", E_error))).append("\t")
//                    .append(String.format(Locale.ROOT, "%1.5f", C));

            writer.append("\n").append(String.format(Locale.ROOT, "%1.8f ", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.8f", E_full)).append("\t")
                    .append((String.format(Locale.ROOT, "%2.8f", E_error))).append("\t")
                    .append(String.format(Locale.ROOT, "%1.8f", M_average)).append("\t   ")
                    .append((String.format(Locale.ROOT, "%1.8f", M_error))).append("\t   ");

            writer_TE.append("\n").append(String.format(Locale.ROOT, "%1.2f", T)).append("   \t")
                    .append(String.format(Locale.ROOT, "%1.1f", E_full));

            writer_TM.append("\n").append(String.format(Locale.ROOT, "%1.1f ", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.1f", M_average));

            writer_TC.append("\n").append(String.format(Locale.ROOT, "%1.2f", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.5f", C));
        }
        writer.flush();
        writer_TC.flush();
        writer_TE.flush();
        writer_TM.flush();
    }
}
