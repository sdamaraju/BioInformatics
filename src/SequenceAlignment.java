import java.io.*;
import java.util.Objects;

public class SequenceAlignment {
    String outputMid = "";
    int calculatedScore = 0;
    enum DIRECTION {
        DIAG,
        TOP,
        LEFT,
        NA
    };

    public static void main(String[] args) throws IOException {
        File file1 = new File(args[0]);
        File file2 = new File(args[1]);
        BufferedReader br = new BufferedReader(new FileReader(file1));
        String line1 = br.readLine();
        br = new BufferedReader(new FileReader(file2));
        String line2 = br.readLine();
        String str1 = line1;
        String str2 = line2;

        //String str1 = "ATCGTAC";
        //String str2 = "ATGTTAT";

        SequenceAlignment align = new SequenceAlignment();
        CustomObj[][] dpMatrix = align.createInitialDPMatrix(str1, str2);
        dpMatrix = align.populateDPMatrix(dpMatrix, str1, str2);
        //System.out.println("DP populated Grid looks like below :");
        //align.printDPMatrix(dpMatrix);
        align.traceBack(dpMatrix, str1, str2);
        System.out.println("\n\nFinal output is : ");
        align.printFinalOutput(str1, align.outputMid, str2);
    }

     class CustomObj {
        int value;
        String direction;

        CustomObj(int value, String direction) {
            this.value = value;
            this.direction = direction;
        }

        @Override
        public String toString() {
            return "{" + value +
                    "," + direction + "}";
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            CustomObj customObj = (CustomObj) o;
            return value == customObj.value;
        }

        @Override
        public int hashCode() {
            return Objects.hash(value);
        }
    }

    private CustomObj[][] createInitialDPMatrix(String str1, String str2) {
        CustomObj[][] dpMatrix = new CustomObj[str1.length() + 1][str2.length() + 1];
        int j = 0;
        for (int i = 0; i < str1.length() + 1; i++) dpMatrix[i][j] = new CustomObj(i, DIRECTION.NA.toString());
        j = 0;
        for (int i = 0; i < str2.length() + 1; i++) dpMatrix[j][i] = new CustomObj(i, DIRECTION.NA.toString());
        return dpMatrix;
    }

    private CustomObj[][] populateDPMatrix(CustomObj[][] dpmatrix, String str1, String str2) {
        for (int i = 1; i < dpmatrix.length; i++) {
            for (int j = 1; j < dpmatrix[0].length; j++) {
                int diag = str1.charAt(i - 1) == str2.charAt(j - 1) ? dpmatrix[i - 1][j - 1].value : dpmatrix[i - 1][j - 1].value + 1;
                int left = dpmatrix[i][j - 1].value + 1;
                int top = dpmatrix[i - 1][j].value + 1;
                int minValue = -1;
                if (diag <= top && diag <= left) {
                    minValue = diag;
                    dpmatrix[i][j] = new CustomObj(minValue, DIRECTION.DIAG.toString());
                } else if (top <= left) {
                    minValue = top;
                    dpmatrix[i][j] = new CustomObj(minValue, DIRECTION.TOP.toString());
                } else {
                    minValue = left;
                    dpmatrix[i][j] = new CustomObj(minValue, DIRECTION.LEFT.toString());
                }

            }
        }
        return dpmatrix;
    }

    private void printDPMatrix(CustomObj[][] dpmatrix) {
        for (int i = 0; i < dpmatrix.length; i++) {
            System.out.print("\n");
            for (int j = 0; j < dpmatrix[0].length; j++) {
                System.out.print(dpmatrix[i][j] + " ");
            }
        }
    }

    private void traceBack(CustomObj[][] dpmatrix, String str1, String str2) {
        int i = dpmatrix.length;
        int j = dpmatrix[0].length;
        getMidOutput(dpmatrix, i - 1, j - 1, str1, str2);
    }

    private void getMidOutput(CustomObj[][] dpmatrix, int i, int j, String str1, String str2) {
        if (i == 0 || j == 0) return;
        if (dpmatrix[i][j].direction == DIRECTION.DIAG.toString()) {
            getMidOutput(dpmatrix, i - 1, j - 1, str1, str2);
            if (str1.charAt(i - 1) == str2.charAt(j - 1)) {
                this.outputMid = this.outputMid + "|";
                calculatedScore =  calculatedScore + 1;
            }
            else {
                this.outputMid = this.outputMid + "*";
                calculatedScore = calculatedScore -2;
            }
        } else if (dpmatrix[i][j].direction == DIRECTION.LEFT.toString()) {
            getMidOutput(dpmatrix, i, j - 1, str1, str2);
            calculatedScore = calculatedScore - 3;
            this.outputMid = this.outputMid + "<";
        } else {
            getMidOutput(dpmatrix, i - 1, j, str1, str2);
            calculatedScore = calculatedScore - 3;
            this.outputMid = this.outputMid + "^";
        }
    }

    private void printFinalOutput(String str1, String mid, String str2) {
        for (int i = 0; i < mid.length(); i++) {
            if (mid.charAt(i) == '<') {
                str1 = str1.substring(0, i) + "_" + str1.substring(i);
            } else if (mid.charAt(i) == '^') {
                str2 = str2.substring(0, i) + "_" + str2.substring(i);
            }
        }
        mid = mid.replaceAll("<", " ");
        mid = mid.replaceAll("\\^", " ");
        System.out.println(str1);
        System.out.println(mid);
        System.out.println(str2);
        System.out.println("Calculated Score is " + calculatedScore);
    }

}
