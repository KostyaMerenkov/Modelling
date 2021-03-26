package students.dvfu.ru;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Locale;
import java.util.Random;

public class Main {
    private static final int size = 32;
    private static final int N = size * size;
    private static final double error = 5;
    final static Random random = new Random();
    //final static DecimalFormat df = new DecimalFormat("#.##");

    private static double[][] createAndFillArray() {
        double[][] arr = new double[size][size];
        //int[] value = {1, -1};
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                arr[i][j] = nextPi();
            }
        }
        return arr;
    }

    private static double nextPi() {
        double x = random.nextDouble();
        return x * (-Math.PI) + (1 - x) * Math.PI;
        //return Double.parseDouble(String.format(Locale.ROOT, "%1.2f",x * (-Math.PI) + (1 - x) * Math.PI));
    }

    private static void printArray(double[][] arr) {
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                System.out.print(arr[i][j] + "\t");
            }
            System.out.println();
        }
        System.out.println();
    }

    private static double[] getNeighbours(double[][] arr, int x, int y) {
        double[] neighbour = new double[4]; //up, right, down, left
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

    private static void printNeighbours(double element, double[] neighbour) {
        System.out.println(" \t" + neighbour[0]);
        System.out.println(neighbour[3] + "\t" + element + "\t" + neighbour[1]);
        System.out.println(" \t" + neighbour[2]);
    }

    private static double calcFullEnergy(double[][] arr) {
        double E_FULL = 0;
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                E_FULL += calcEnergy(arr, i, j);
            }
        }
        return E_FULL / N;
    }

    private static double calcEnergy(double[][] arr, int x, int y) {
        double[] neighbour = getNeighbours(arr, x, y);
        double E = 0;
        for (int i = 0; i < neighbour.length; i++) {
            E += Math.cos(arr[x][y] - neighbour[i]);
        }
        E *= -1;
        return E;
    }

    private static double calcMagnetization(double[][] arr) {
        double msin = 0, mcos = 0;
        double M = 0;
        double M_old = 0;
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                mcos += Math.cos(arr[i][j]);
                msin += Math.sin(arr[i][j]);
            }
        }

        M = (mcos*mcos + msin*msin) / (N*N);
        return M;
        //return M;
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
        double[][] arr = createAndFillArray();
        printArray(arr);
        System.out.println();
        FileWriter writer = new FileWriter(new File("Results.txt"));
        writer.append("T\t\tM\t\t\tm\t\tE\t\te\t\tC");
        FileWriter writer_TE = new FileWriter(new File("TE.txt"));
        writer_TE.append("T\t\t E");
        FileWriter writer_TM = new FileWriter(new File("TM.txt"));
        writer_TM.append("T\t\tM");
        FileWriter writer_TC = new FileWriter(new File("TC.txt"));
        writer_TM.append("T\t\tC");

        /*
        for (double T = 0.01; T <= 0.01; T++) {
            double E_average = 0;
            double E_average_squared = 0;
            double[] E_av = new double[(int)error];
            double[] M_av = new double[(int)error];
            double E_full = 0;
            double M_average = 0;
            double C;
            for (int average = 0; average < error; average++) {
                for (int MK = 0; MK < 1000; MK++) {
                    int ri = random.nextInt(size), rj = random.nextInt(size);
                    double E1 = calcEnergy(arr, ri, rj);
                    double old_value = arr[ri][rj];
                    arr[ri][rj] = nextPi();
                    double E2 = calcEnergy(arr, ri, rj);
                    if (E2 >= E1) {
                        if (nextPi() > Math.exp((-1 * (E2 - E1)) / T)) {
                            arr[ri][rj] = old_value;
                        }
                    }
                }
                printArray(arr);
                E_full = calcFullEnergy(arr);
                E_av[average] = E_full;
                E_average += E_full;
                E_average_squared += (E_full * E_full);
                M_av[average] = calcMagnetization(arr);
                M_average += M_av[average];
            }
            double E_error = calcError(E_av, E_average);
            double M_error = calcError(M_av, M_average);
            M_average /= error;
            C = calcThermal(E_average, E_average_squared, T);
            writer.append("\n").append(String.format(Locale.ROOT, "%1.2f", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.1f", M_average)).append("\t   ")
                    .append((String.format(Locale.ROOT, "%1.2f", M_error))).append("\t   ")
                    .append(String.format(Locale.ROOT, "%1.1f", E_full)).append("\t    ")
                    .append((String.format(Locale.ROOT, "%2.2f", E_error))).append("\t")
                    .append(String.format(Locale.ROOT, "%1.5f", C));

            writer_TE.append("\n").append(String.format(Locale.ROOT, "%1.2f", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.1f", E_full)).append("\t\t");

            writer_TM.append("\n").append(String.format(Locale.ROOT, "%1.2f", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.1f", M_average));

            writer_TC.append("\n").append(String.format(Locale.ROOT, "%1.2f", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.5f", C));
        }
*/

        for (double T = 0.1, j = 0; T <= 4; T += 0.1, j++) {
            FileWriter writer_T = new FileWriter(new File("S/S-" + getBinary7(Integer.toBinaryString((int)j)) + ".txt"));
            massiveToParaview(writer_T, arr);
            double E_average = 0;
            double E_average_squared = 0;
            double[] E_av = new double[(int)error];
            double[] M_av = new double[(int)error];
            double E_full = 0;
            double M_average = 0;
            double C = 0;
            for (int average = 0; average < error; average++) {
                for (int MK = 0; MK < 100000; MK++) {
                    int r_i = random.nextInt(size), r_j = random.nextInt(size);
                    double E1 = calcEnergy(arr, r_i, r_j);
                    double old_value = arr[r_i][r_j];
                    arr[r_i][r_j] = nextPi();
                    double E2 = calcEnergy(arr, r_i, r_j);
                    if (E2 >= E1) {
                        double randExp = Math.exp((-1 * (E2 - E1)) / T);
                        double rand = random.nextDouble();
                        if (rand > randExp) {
                            arr[r_i][r_j] = old_value;
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
            }
            double E_error = calcError(E_av, E_average);
            double M_error = calcError(M_av, M_average);
            M_average /= error;
            C = calcThermal(E_average, E_average_squared, T);
            writer.append("\n").append(String.format(Locale.ROOT, "%1.1f ", T)).append("\t\t")
                    .append(String.format(Locale.ROOT, "%1.4f", M_average)).append("\t   ")
                    .append((String.format(Locale.ROOT, "%1.2f", M_error))).append("\t   ")
                    .append(String.format(Locale.ROOT, "%1.1f", E_full)).append("\t    ")
                    .append((String.format(Locale.ROOT, "%2.2f", E_error))).append("\t")
                    .append(String.format(Locale.ROOT, "%1.5f", C));

            writer_TE.append("\n").append(String.format(Locale.ROOT, "%1.2f", T)).append("   \t")
                    .append(String.format(Locale.ROOT, "%1.1f", E_full));

            writer_TM.append("\n").append(String.format(Locale.ROOT, "%1.1f ", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.1f", M_average));

            writer_TC.append("\n").append(String.format(Locale.ROOT, "%1.2f", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.5f", C));
        }
        printArray(arr);
        writerFlushAndClose(writer, writer_TC, writer_TE, writer_TM);


    }

    private static void massiveToParaview(FileWriter writer_t, double[][] arr) throws IOException {
        //writer_t.append("x\ty\tz\tX\tY\tZ\n");
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                writer_t.append(String.valueOf(i)).append("\t").append(String.valueOf(j)).append("\t")
                        .append("0").append("\t")
                        .append(String.format(Locale.ROOT, "%1.2f", Math.cos(arr[i][j]))).append("\t")
                        .append(String.format(Locale.ROOT, "%1.2f", Math.sin(arr[i][j]))).append("\t")
                        .append("0\n");
            }
        }
        writer_t.flush();
        writer_t.close();
    }

    private static String getBinary7(String toBinaryString) {
        StringBuilder toBinaryStringBuilder = new StringBuilder(toBinaryString);
        while (toBinaryStringBuilder.length() < 7)
        {
            toBinaryStringBuilder.insert(0, "0");
        }
        toBinaryString = toBinaryStringBuilder.toString();
        return toBinaryString;
    }

    private static void writerFlushAndClose(FileWriter writer, FileWriter writer_tc, FileWriter writer_te, FileWriter writer_tm) throws IOException {
        writer.flush();
        writer_tc.flush();
        writer_te.flush();
        writer_tm.flush();
        writer.close();
        writer_tc.close();
        writer_te.close();
        writer_tm.close();
    }
}
