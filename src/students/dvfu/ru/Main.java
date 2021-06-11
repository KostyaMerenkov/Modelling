package students.dvfu.ru;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Locale;
import java.util.Random;

public class Main {
    private static final int size = 10;
    private static final int N = size * size;
    private static final double error = 5;
    final static Random random = new Random();
    //final static DecimalFormat df = new DecimalFormat("#.##");

    private static double[][] Sx;
    private static double[][] Sy;
    private static double[][] Sz;


    private static void createAndFillArray() {
        Sx = new double[size][size];
        Sy = new double[size][size];
        Sz = new double[size][size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                double psi = - Math.acos(2 * random.nextDouble() - 1);
                double fi = Math.PI*2* random.nextDouble();
                Sx[i][j] = Math.sin(fi)*Math.cos(psi);
                Sy[i][j] = Math.sin(fi)*Math.sin(psi);
                Sz[i][j] = Math.cos(fi);
            }
        }
    }

//    private static double nextPi() {
//        double x = random.nextDouble();
//        return x * (-Math.PI) + (1 - x) * Math.PI;
//        //return Double.parseDouble(String.format(Locale.ROOT, "%1.2f",x * (-Math.PI) + (1 - x) * Math.PI));
//    }

    private static void printArray(double[][] arr) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
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

    private static double calcFullEnergy() {
        double E_FULL = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                E_FULL += calcEnergy(i, j);
            }
        }
        return E_FULL / N;
    }

    private static double calcEnergy(int x, int y) {
        double[] neighbourX = getNeighbours(Sx, x, y);
        double[] neighbourY = getNeighbours(Sy, x, y);
        double[] neighbourZ = getNeighbours(Sz, x, y);
        double E = 0;
        for (int i = 0; i < 4; i++) {
            E += (Sx[x][y]*neighbourX[i]+Sy[x][y]*neighbourY[i]+Sy[x][y]*neighbourZ[i]);
        }
        E *= -1;
        return E;
    }

    private static double calcMagnetization() {
        double M = 0;
        double Mx = 0, My = 0, Mz = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                Mx += Sx[i][j];
                My += Sy[i][j];
                Mz += Sz[i][j];
            }
        }
        M = (Mx*Mx + My*My + Mz*Mz) / (N);
        return M;
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
        createAndFillArray();
        //printArray(arr);
        System.out.println();
        FileWriter writer = new FileWriter(new File("Results.txt"));
        writer.append("T\t\t M\t\t\tm\t\t\tE\t\t\te\t\t\tC");
        FileWriter writer_TE = new FileWriter(new File("TE.txt"));
        writer_TE.append("T\t\t E");
        FileWriter writer_TM = new FileWriter(new File("TM.txt"));
        writer_TM.append("T\t\tM");
        FileWriter writer_TC = new FileWriter(new File("TC.txt"));
        writer_TC.append("T\t\tC");

        for (double T = 0.1, j = 0; T <= 4; T += 0.1, j++) {
            FileWriter writer_T = new FileWriter(new File("S/S-" + getBinary7(Integer.toBinaryString((int)j)) + ".txt"));
            massiveToParaview(writer_T);
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
                    double E1 = calcEnergy(r_i, r_j);
                    double old_x = Sx[r_i][r_j];
                    double old_y = Sy[r_i][r_j];
                    double old_z = Sz[r_i][r_j];
                    double psi = - Math.acos(2 * random.nextDouble() - 1);
                    double fi = Math.PI*2 * random.nextDouble();
                    Sx[r_i][r_j] = Math.sin(fi)*Math.cos(psi);
                    Sy[r_i][r_j] = Math.sin(fi)*Math.sin(psi);
                    Sz[r_i][r_j] = Math.cos(fi);
                    double E2 = calcEnergy(r_i, r_j);
                    if (E2 >= E1) {
                        double randExp = Math.exp((-1 * (E2 - E1)) / T);
                        double rand = random.nextDouble();
                        if (rand > randExp) {
                            Sx[r_i][r_j] = old_x;
                            Sy[r_i][r_j] = old_y;
                            Sz[r_i][r_j] = old_z;
                        }
                    }
                }
                //printArray(arr);
                E_full = calcFullEnergy();
                E_av[average] = E_full;
                E_average += E_full;
                E_average_squared += (E_full * E_full);
                M_av[average] = calcMagnetization();
                M_average += M_av[average];
            }
            double E_error = calcError(E_av, E_average);
            double M_error = calcError(M_av, M_average);
            M_average /= error;
            C = calcThermal(E_average, E_average_squared, T);
            writer.append("\n").append(String.format(Locale.ROOT, "%1.1f ", T)).append("\t")
                    .append(String.format(Locale.ROOT, "%1.4f", M_average)).append("\t\t")
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
        //printArray(arr);
        writerFlushAndClose(writer, writer_TC, writer_TE, writer_TM);


    }

    private static void massiveToParaview(FileWriter writer_t) throws IOException {
        //writer_t.append("x \t y \t z \t X \t Y \t Z \n");
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                writer_t.append(String.valueOf(i)).append("\t").append(String.valueOf(j)).append("\t")
                        .append("1").append("\t")
                        .append(String.format(Locale.ROOT, "%1.2f", Sx[i][j])).append("\t")
                        .append(String.format(Locale.ROOT, "%1.2f", Sy[i][j])).append("\t")
                        .append(String.format(Locale.ROOT, "%1.2f", Sz[i][j]));
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
