package students.dvfu.ru;

import java.io.*;
import java.util.Random;

public class Main {
    private static final int size = 30;
    private static final int N = size*size;
    final static Random random = new Random();

    private static int[][] createAndFillArray() {
        int[][] arr = new int[size][size];
        int[] value = {1,-1};
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                arr[i][j] = value[(int)(Math.random()*value.length)];
            }
        }
        return arr;
    }

    private static void printArray(int[][] arr) {
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                System.out.print(arr[i][j]+"\t");
            }
            System.out.println();
        }
        System.out.println();
    }

    private static int[] getNeighbours(int[][] arr, int x, int y) {
        int[] neighbour = new int[4]; //up, right, down, left
        if (x == 0) { //up
            neighbour[0] = arr[arr.length-1][y];
        } else {
            neighbour[0] = arr[x-1][y];
        }
        if (y == arr.length-1) { //right
            neighbour[1] = arr[x][0];
        } else {
            neighbour[1] = arr[x][y+1];
        }
        if (x == arr.length-1) { //down
            neighbour[2] = arr[0][y];
        } else {
            neighbour[2] = arr[x+1][y];
        }
        if (y == 0) { //left
            neighbour[3] = arr[x][arr.length-1];
        } else {
            neighbour[3] = arr[x][y-1];
        }
        return neighbour;
    }

    private static void printNeighbours(int element, int[] neighbour) {
        System.out.println(" \t" + neighbour[0]);
        System.out.println(neighbour[3] + "\t" + element + "\t"+ neighbour[1]);
        System.out.println(" \t" + neighbour[2]);
    }

    private static double calcFullEnergy(int[][] arr) {
        float E_FULL = 0;
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr.length; j++) {
                E_FULL += calcEnergy(arr,i,j);
            }
        }
        return E_FULL/N;
    }

    private static double calcEnergy(int[][] arr, int x, int y) {
        int[] neighbour = getNeighbours(arr, x, y);
        float E = 0;
        for (int i=0; i < neighbour.length; i++) {
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
        return Math.abs(M/N);
    }

    private static double calcThermal(double E, double T) {
        return (Math.abs(Math.pow(E,2))-Math.pow(Math.abs(E),2))/Math.pow(T,2);
    }

    public static void main(String[] args) throws IOException {
        int[][] arr = createAndFillArray();
        printArray(arr);
        System.out.println();
        FileWriter writer = new FileWriter(new File("Results.txt"));
        writer.append("T\t\tM\t\tE\t\tC");

        FileWriter writer_TE = new FileWriter(new File("TE.txt"));
        writer_TE.append("T\t\t E");

        FileWriter writer_TM = new FileWriter(new File("TM.txt"));
        writer_TM.append("T\t\tM");

        for (double T=0.01; T<=0.01; T++){
            for (int MK = 0; MK < 100000; MK++) {
                int ri = random.nextInt(size), rj = random.nextInt(size);
                double E1 = calcEnergy(arr, ri, rj);
                arr[ri][rj] = arr[ri][rj] * -1;
                double E2 = calcEnergy(arr, ri, rj);
                if (E2 >= E1) {
                    if (random.nextDouble() > Math.exp((-1*(E2-E1))/T)) {
                        arr[ri][rj] = arr[ri][rj] * -1;
                    }
                }
            }
            printArray(arr);
            double E_FULL = calcFullEnergy(arr);
            double M = calcMagnetization(arr);
            double C = calcThermal(E_FULL, T);
            writer.append("\n").append(String.format("%1.2f", T)).append("\t").append(String.format("%1.1f", M))
                    .append("\t   ").append(String.format("%1.1f", E_FULL)).append("\t    ")
                    .append(String.format("%1.1f", C));
            writer_TE.append("\n").append(String.format("%1.2f", T)).append("\t")
                    .append(String.format("%1.1f", E_FULL)).append("\t\t");
            writer_TM.append("\n").append(String.format("%1.2f", T)).append("\t").append(String.format("%1.1f", M));
        }

        for (double T = 0.1; T <= 4; T+=0.1) {

            for (int MK = 0; MK < 1000000; MK++) {
                int r_i = random.nextInt(size), r_j = random.nextInt(size);
                double E1 = calcEnergy(arr, r_i, r_j);
                arr[r_i][r_j] = (arr[r_i][r_j]) * -1;
                double E2 = calcEnergy(arr, r_i, r_j);
                if (E2>=E1) {
                    double randExp = Math.exp((-1*(E2-E1))/T);
                    if(random.nextDouble() > randExp) {
                        arr[r_i][r_j] = (arr[r_i][r_j]) * -1;
                    }
                }
            }
            printArray(arr);
            double E_FULL = calcFullEnergy(arr);
            double M = calcMagnetization(arr);
            double C = calcThermal(E_FULL, T);
            writer.append("\n").append(String.format("%1.1f ", T)).append("\t").append(String.format("%1.1f", M))
                    .append("\t   ").append(String.format("%1.1f", E_FULL)).append("\t    ")
                    .append(String.format("%1.1f", C)).flush();

            writer_TE.append("\n").append(String.format("%1.2f", T)).append("   \t")
                    .append(String.format("%1.1f", E_FULL)).flush();

            writer_TM.append("\n").append(String.format("%1.1f ", T)).append("\t").append(String.format("%1.1f", M))
                    .flush();
        }
    }
}
