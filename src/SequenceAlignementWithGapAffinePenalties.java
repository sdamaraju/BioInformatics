import java.io.*;
import java.util.Objects;

public class SequenceAlignementWithGapAffinePenalties {

    String outputMid = "";
    int calculatedScore = 0;
    int gapStart = -5;
    int gapExtend = -1;
    int matched = 1;
    int unMatched = -2;
    boolean previousGap = false;

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

        SequenceAlignementWithGapAffinePenalties align = new SequenceAlignementWithGapAffinePenalties();
        Object initialMatrices[] = align.createInitialDPMatrices(str1, str2);
        /*for (int i = 0; i < initialMatrices.length; i++) {
            align.printDPMatrix((CustomObj[][]) initialMatrices[i]);
            System.out.println("");
        }*/
        initialMatrices = align.populateDPMatrices(initialMatrices, str1, str2);
        /*for (int i = 0; i < initialMatrices.length; i++) {
            align.printDPMatrix((CustomObj[][]) initialMatrices[i]);
            System.out.println("");
        }*/
        align.traceBack((CustomObj[][]) initialMatrices[0], str1, str2);
        System.out.println("\n\nFinal output is : ");
        align.printFinalOutput(str1, align.outputMid, str2);
    }

    private Object[] populateDPMatrices(Object dpMatrices[], String str1, String str2) {
        CustomObj[][] mMatrix = (CustomObj[][]) dpMatrices[0];
        CustomObj[][] xMatrix = (CustomObj[][]) dpMatrices[1];
        CustomObj[][] yMatrix = (CustomObj[][]) dpMatrices[2];
        for (int i = 1; i < mMatrix.length; i++) {
            for (int j = 1; j < mMatrix[0].length; j++) {
                //Matrix M
                int match = str1.charAt(i - 1) == str2.charAt(j - 1) ? matched : unMatched;
                int diag = mMatrix[i - 1][j - 1].value != Integer.MIN_VALUE ? match + mMatrix[i - 1][j - 1].value : Integer.MIN_VALUE;
                int left = xMatrix[i - 1][j - 1].value != Integer.MIN_VALUE ? match + xMatrix[i - 1][j - 1].value : Integer.MIN_VALUE;
                int top = yMatrix[i - 1][j - 1].value != Integer.MIN_VALUE ? match + yMatrix[i - 1][j - 1].value : Integer.MIN_VALUE;
                /*System.out.println("M");
                System.out.println(diag);
                System.out.println(left);
                System.out.println(top);
                */int maxValue = -1;
                if (diag >= top && diag >= left) {
                    maxValue = diag;
                    mMatrix[i][j] = new CustomObj(maxValue, DIRECTION.DIAG.toString());
                } else if (top >= left) {
                    maxValue = top;
                    mMatrix[i][j] = new CustomObj(maxValue, DIRECTION.TOP.toString());
                } else {
                    maxValue = left;
                    mMatrix[i][j] = new CustomObj(maxValue, DIRECTION.LEFT.toString());
                }
                // Matrix X
                diag = mMatrix[i][j - 1].value != Integer.MIN_VALUE ? gapStart + gapExtend + mMatrix[i][j - 1].value : Integer.MIN_VALUE;
                left = xMatrix[i][j - 1].value != Integer.MIN_VALUE ? gapExtend + xMatrix[i][j - 1].value : Integer.MIN_VALUE;
                top = yMatrix[i][j - 1].value != Integer.MIN_VALUE ? gapStart + gapExtend + yMatrix[i][j - 1].value : Integer.MIN_VALUE;
                /*System.out.println("X");
                System.out.println(diag);
                System.out.println(left);
                System.out.println(top);
                */maxValue = -1;
                if (diag >= top && diag >= left) {
                    maxValue = diag;
                    xMatrix[i][j] = new CustomObj(maxValue, DIRECTION.DIAG.toString());
                } else if (top >= left) {
                    maxValue = top;
                    xMatrix[i][j] = new CustomObj(maxValue, DIRECTION.TOP.toString());
                } else {
                    maxValue = left;
                    xMatrix[i][j] = new CustomObj(maxValue, DIRECTION.LEFT.toString());
                }
                // Matrix Y
                diag = mMatrix[i - 1][j].value != Integer.MIN_VALUE ? gapStart + gapExtend + mMatrix[i-1][j].value : Integer.MIN_VALUE;
                left = xMatrix[i - 1][j].value != Integer.MIN_VALUE ? gapStart + gapExtend + xMatrix[i-1][j].value : Integer.MIN_VALUE;
                top = yMatrix[i - 1][j].value != Integer.MIN_VALUE ? gapExtend + yMatrix[i-1][j].value : Integer.MIN_VALUE;
                /*System.out.println(diag);
                System.out.println(left);
                System.out.println(top);
*/
                maxValue = -1;
                if (diag >= top && diag >= left) {
                    maxValue = diag;
                    yMatrix[i][j] = new CustomObj(maxValue, DIRECTION.DIAG.toString());
                } else if (top >= left) {
                    maxValue = top;
                    yMatrix[i][j] = new CustomObj(maxValue, DIRECTION.TOP.toString());
                } else {
                    maxValue = left;
                    yMatrix[i][j] = new CustomObj(maxValue, DIRECTION.LEFT.toString());
                }
            }
        }
        return dpMatrices;
    }


    private void printDPMatrix(CustomObj[][] dpmatrix) {
        for (int i = 0; i < dpmatrix.length; i++) {
            System.out.print("\n");
            for (int j = 0; j < dpmatrix[0].length; j++) {
                System.out.print(dpmatrix[i][j] + " ");
            }
        }
    }

    private Object[] createInitialDPMatrices(String str1, String str2) {
        //M
        CustomObj[][] mMatrix = new CustomObj[str1.length() + 1][str2.length() + 1];
        int j = 0;
        for (int i = 0; i < str1.length() + 1; i++)
            mMatrix[i][j] = new CustomObj(Integer.MIN_VALUE, DIRECTION.NA.toString());
        j = 0;
        for (int i = 0; i < str2.length() + 1; i++)
            mMatrix[j][i] = new CustomObj(Integer.MIN_VALUE, DIRECTION.NA.toString());
        // X and Y
        CustomObj[][] xMatrix = new CustomObj[str1.length() + 1][str2.length() + 1];
        j = 0;
        for (int i = 0; i < str2.length() + 1; i++)
            xMatrix[j][i] = new CustomObj(gapStart + (i * gapExtend), DIRECTION.NA.toString());
        j = 0;
        for (int i = 0; i < str1.length() + 1; i++)
            xMatrix[i][j] = new CustomObj(Integer.MIN_VALUE, DIRECTION.NA.toString());

        CustomObj[][] yMatrix = new CustomObj[str1.length() + 1][str2.length() + 1];
        j = 0;
        for (int i = 0; i < str1.length() + 1; i++)
            yMatrix[i][j] = new CustomObj(gapStart + (i * gapExtend), DIRECTION.NA.toString());
        j = 0;
        for (int i = 0; i < str2.length() + 1; i++)
            yMatrix[j][i] = new CustomObj(Integer.MIN_VALUE, DIRECTION.NA.toString());
        return new Object[]{mMatrix, xMatrix, yMatrix};
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
            SequenceAlignment.CustomObj customObj = (SequenceAlignment.CustomObj) o;
            return value == customObj.value;
        }

        @Override
        public int hashCode() {
            return Objects.hash(value);
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
            previousGap = false;
            calculatedScore = str1.charAt(i - 1) == str2.charAt(j - 1) ? calculatedScore + matched : calculatedScore + unMatched;
            if (str1.charAt(i - 1) == str2.charAt(j - 1)) this.outputMid = this.outputMid + "|";
            else this.outputMid = this.outputMid + "*";
        } else if (dpmatrix[i][j].direction == DIRECTION.LEFT.toString()) {
            getMidOutput(dpmatrix, i, j - 1, str1, str2);
            if(!previousGap){
                previousGap=true;
                calculatedScore = calculatedScore -5;
            }
            calculatedScore = calculatedScore - 1;
            this.outputMid = this.outputMid + "<";
        } else {
            getMidOutput(dpmatrix, i - 1, j, str1, str2);
            if(!previousGap){
                previousGap=true;
                calculatedScore = calculatedScore -5;
            }
            calculatedScore = calculatedScore - 1;
            this.outputMid = this.outputMid + "^";
        }
    }

    private void printFinalOutput(String str1, String mid, String str2) {
        for (int i = 0; i < mid.length(); i++) {
            if (mid.charAt(i) == '<') {
                str1 = str1.substring(0, i) + "_" + str1.substring(i);
            } else if (mid.charAt(i) == '^') str2 = str2.substring(0, i) + "_" + str2.substring(i);
        }
        mid = mid.replaceAll("<", " ");
        mid = mid.replaceAll("\\^", " ");
        System.out.println(str1);
        System.out.println(mid);
        System.out.println(str2);
        System.out.println("Calculated Score is " + calculatedScore);
    }

}
