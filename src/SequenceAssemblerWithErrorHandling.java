import java.io.*;
import java.util.ArrayList;

public class SequenceAssemblerWithErrorHandling {
    static int thresholdForOverlapMatch = 51; // given from problem.
    static int errorThreshold = 1; // given from problem (1 percent). // yet to be implemented

    public static void main(String[] args) throws IOException {
        File file = new File("/Users/sdamaraju/IdeaProjects/BioInformatics/files/HW2_reads.fasta");
        File opfile = new File("/Users/sdamaraju/Desktop/730HW2.fasta");
        FileOutputStream fos = new FileOutputStream(opfile);
        if (!opfile.exists()) {
            opfile.createNewFile();
        }
        try {
            ArrayList<String> reads = new ArrayList();
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = br.readLine();
            while (line != null) {
                if (!line.startsWith(">")) {
                    reads.add(line);
                }
                line = br.readLine();
            }
            /*reads = new ArrayList();
            reads.add("abcdefghk");
            reads.add("defghijk");
            reads.add("ghijklmn");
            reads.add("mnopqrst");
            reads.add("jklmnopq");*/
            // Step 1 : Identify the first String
            int indexOfFirstString = identifyFirstPartOfTheFullSequence(reads);
            // Step 2 : Store first String
            String superString = reads.get(indexOfFirstString);
            // Step 3 : Remove the  first string from reads.
            reads.remove(indexOfFirstString);
            // Step 4: Compare with all other strings to find the best merging match. => maximum overlap.
            superString = identifyTheSuperString(superString, reads);
            // Output
            System.out.println("Merged Super String is :");
            fos.write(">\n".getBytes());
            fos.write(superString.getBytes());
            System.out.println(superString);
            System.out.println("Length of superString is " + superString.length());
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (fos != null) {
                    fos.close();
                }
            } catch (IOException ioe) {
                System.out.println("Error in closing the Stream");
            }
        }
    }

    private static int identifyFirstPartOfTheFullSequence(ArrayList<String> reads) {
        // if its the first string, it should not be part of any other string
        for (int i = 0; i < reads.size(); i++) {
            String str1 = reads.get(i);
            int count = 0;
            for (int j = 0; j < reads.size(); j++) {
                if (i == j) continue;
                String str2 = reads.get(j);
                if (isValidOverlapPrefixSuffixSequence(str2, str1) > thresholdForOverlapMatch) {
                    break;
                } else {
                    count++;
                }
                if (count == reads.size() - 1) {
                    return i;
                }
            }
        }
        return -1;
    }

    private static String identifyTheSuperString(String superString, ArrayList<String> reads) {
        int totalStringsToMergeIntoFirstString = reads.size();
        for (int i = 0; i < totalStringsToMergeIntoFirstString; i++) {
            int minIndexForMaxOverlap = Integer.MAX_VALUE;
            int indexOfReadsWithMaxOverlap = -1;
            for (int j = 0; j < reads.size(); j++) {
                String sequenceToMatch = reads.get(j);
                int overlapLength = isValidOverlapPrefixSuffixSequence(superString, sequenceToMatch);
                if (overlapLength >= thresholdForOverlapMatch) {
                    int indexOfOverlap = superString.length() - overlapLength;
                    if (indexOfOverlap < minIndexForMaxOverlap) {
                        minIndexForMaxOverlap = indexOfOverlap;
                        indexOfReadsWithMaxOverlap = j;
                    }
                }
            }
            // We identify the index of reads with maximum overlap
            // now remove that identified string from reads.
            if (indexOfReadsWithMaxOverlap != -1) {
                String stringToMerge = reads.get(indexOfReadsWithMaxOverlap);
                int indexOfMerge = superString.indexOf(stringToMerge.substring(0, thresholdForOverlapMatch));
                superString = superString.substring(0, indexOfMerge).concat(stringToMerge);
                reads.remove(indexOfReadsWithMaxOverlap);
            }
        }
        return superString;
    }

    private static int isValidOverlapPrefixSuffixSequence(String str1, String str2) {
        if (str1 == null || str2 == null || str1.length() == 0 || str2.length() == 0) return -1;
        char firstCharOfTargetSequence = str2.charAt(0);
        int nextPossibleStartIndex = -1;
        int k = 0;
        int overlapMatch = 0;
        int errorCount = 0;
        for (int i = 0; i < str1.length(); i++) {
            if (str2.length() <= k) return -1;
            if (str1.charAt(i) == str2.charAt(k) || errorCount < errorThreshold) {
                if (str1.charAt(i) != str2.charAt(k)) {
                    errorCount++;
                }
                overlapMatch++;
                k++;
                if (str1.charAt(i) == firstCharOfTargetSequence && (k - 1 != 0) && nextPossibleStartIndex == -1) {
                    nextPossibleStartIndex = i;
                }
            } else {
                overlapMatch = 0;
                errorCount = 0;
                k = 0;
                if (str1.charAt(i) == firstCharOfTargetSequence && (k - 1 != 0) && nextPossibleStartIndex == -1) {
                    nextPossibleStartIndex = i;
                }
                if (nextPossibleStartIndex != -1) {
                    i = nextPossibleStartIndex - 1;
                    nextPossibleStartIndex = -1;
                }
            }

        }
        if (overlapMatch >= thresholdForOverlapMatch) {
            return overlapMatch;
        } else {
            return -1;
        }
    }
}
