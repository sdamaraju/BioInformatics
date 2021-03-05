import java.io.*;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.HashSet;

public class SequenceGenerator {

    static HashMap<String, String> dnaSeq = new HashMap();
    static HashMap<String, String> codonTable = new HashMap();
    static HashSet<String> nonExistingChromozome = new HashSet();

    public static void main(String[] args) {
        LocalDateTime curr1 = LocalDateTime.now();
        buildCodonTable();
        File file = new File("/Users/sdamaraju/IdeaProjects/BioInformatics/files/730Proj1.csv");
        FileOutputStream fos = null;
        File opfile;
        int count = 1;
        try {
            //Specify the file path here
            opfile = new File("/Users/sdamaraju/Desktop/730HW1.fa");
            fos = new FileOutputStream(opfile);
            if (!opfile.exists()) {
                opfile.createNewFile();
            }
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = br.readLine();

            while (line != null) {
                boolean cdsStartUsed = false;
                boolean cdsEndUsed = false;
                if (count != 1) { //ignore the first header line in file
                    line = line.replace("\"", "");
                    line = line.replace(",,", ",");
                    String[] lineValues = line.split(",");
                    int exonCount = Integer.parseInt(lineValues[8]);
                    int firstNonNegativeExonFrameCounter = getFirstNonNegativeExonFrameCounter(getExonFrames(lineValues));
                    if (firstNonNegativeExonFrameCounter == -1/* || lineValues[2].equals("chrM") || lineValues[2].equals("chrX") || lineValues[2].equals("chrY")*/) {
                        line = br.readLine();
                        //count++;
                        continue;
                    }
                    String chromozome = lineValues[2];
                    boolean isValid = populateDNASequence(chromozome);
                    if (!isValid) {
                        line = br.readLine();
                        continue;
                    }
                    byte[] bytesArray1 = (getHeader(lineValues) + "\n").getBytes();
                    fos.write(bytesArray1);
                    String strand = lineValues[3];
                    StringBuffer protienSequence = new StringBuffer();
                    if (strand.equals("+")) {
                        int exonCounter = 1;
                        while (exonCounter <= exonCount) {
                            int exonFrameValue = getExonFrameValue(lineValues, exonCounter);
                            if (isExonValidForSequencing(exonFrameValue)) {
                                int startIndex = cdsStartUsed ? getExonStartIndexForPositiveStrand(exonCounter, lineValues, exonFrameValue) : Integer.parseInt(lineValues[6]);
                                cdsStartUsed = true;
                                int endIndex = getExonEndIndexForPositiveStrand(exonCounter, exonCount, lineValues);
                                String exonSequenceSegment = getExonSequenceSegmentFromSequence(startIndex, endIndex, dnaSeq.get(chromozome));
                                Object[] checkStop_sequence = translateToProteinSequence(exonSequenceSegment);
                                protienSequence.append(checkStop_sequence[1]);
                                if ((Boolean) checkStop_sequence[0]) {
                                    break;
                                }
                            }
                            exonCounter++;
                        }
                    } else {
                        int exonCounter = exonCount;
                        while (exonCounter > 0) {
                            int exonFrameValue = getExonFrameValue(lineValues, exonCounter);
                            if (isExonValidForSequencing(exonFrameValue)) {
                                int startIndex = getExonStartIndexForNegativeStrand(exonCounter, lineValues);
                                int endIndex = cdsEndUsed ? getExonEndIndexForNegativeStrand(exonCounter, exonCount, lineValues, exonFrameValue) : Integer.parseInt(lineValues[7]) - 1;
                                cdsEndUsed = true;
                                String exonSequenceSegment = getExonSequenceSegmentFromSequence(startIndex, endIndex, dnaSeq.get(chromozome));
                                String reverseComplementExonSequence = getReverseComplement(exonSequenceSegment);
                                Object[] checkStop_sequence = translateToProteinSequence(reverseComplementExonSequence);
                                protienSequence.append(checkStop_sequence[1]);
                                if ((Boolean) checkStop_sequence[0]) {
                                    break;
                                }
                            }
                            exonCounter--;
                        }
                    }
                    byte[] bytesArray2 = (protienSequence.toString() + "\n").getBytes();
                    fos.write(bytesArray2);
                }
                line = br.readLine();
                count++;
            }
            fos.flush();
        } catch (Exception ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (fos != null) {
                    fos.close();
                }
            } catch (IOException ioe) {
                System.out.println("Error in closing the Stream");
            }
        }
        System.out.println("There are in a total " + count + " records");
        LocalDateTime curr2 = LocalDateTime.now();
        Duration duration = Duration.between(curr1, curr2);
        System.out.printf("Duration = %s seconds.%n", duration.getSeconds());
    }

    private static String getReverseComplement(String exonSequenceSegment) {
        StringBuffer exonSequence = new StringBuffer(exonSequenceSegment);
        StringBuffer reverseSequence = exonSequence.reverse();
        for (int i = 0; i < reverseSequence.length(); i++) {
            //complement
            if (reverseSequence.charAt(i) == 'A' || reverseSequence.charAt(i) == 'a') reverseSequence.setCharAt(i, 'T');
            else if (reverseSequence.charAt(i) == 'G' || reverseSequence.charAt(i) == 'g')
                reverseSequence.setCharAt(i, 'C');
            else if (reverseSequence.charAt(i) == 'C' || reverseSequence.charAt(i) == 'c')
                reverseSequence.setCharAt(i, 'G');
            else if (reverseSequence.charAt(i) == 'T' || reverseSequence.charAt(i) == 't')
                reverseSequence.setCharAt(i, 'A');
            else if (reverseSequence.charAt(i) == 'N' || reverseSequence.charAt(i) == 'n')
                reverseSequence.setCharAt(i, 'A');
        }
        return reverseSequence.toString();
    }

    private static boolean populateDNASequence(String chromozome) throws IOException {
        if (dnaSeq.get(chromozome) == null) {
            dnaSeq = new HashMap();
            if (nonExistingChromozome.contains(chromozome)) {
                return false;
            }
            String dnaSequence = getDataFromDNASequence(chromozome);
            if (dnaSequence.equals("")) {
                nonExistingChromozome.add(chromozome);
                return false;
            } else dnaSeq.put(chromozome, dnaSequence);
        }
        return true;
    }

    private static boolean isExonValidForSequencing(int exonFrame) {
        if (exonFrame == -1) return false;
        else return true;
    }

    private static int getExonFrameValue(String lineValues[], int exonCounter) {
        int[] exonFrames = getExonFrames(lineValues);
        return exonFrames[exonCounter - 1]; // exonCounter-1 because we pass exoncounter that starts from 1.
    }

    private static int[] getExonFrames(String lineValues[]) {
        int exonCount = Integer.parseInt(lineValues[8]);
        int minIndex = lineValues.length - exonCount; //last n values are exonframes (where n is exonCount)
        int[] exonFrames = new int[exonCount];
        int k = 0;
        for (int i = minIndex; i < lineValues.length; i++) {
            exonFrames[k++] = Integer.parseInt(lineValues[i]);
        }
        return exonFrames;
    }

    private static Object[] translateToProteinSequence(String exonSequenceSegment) {
        StringBuffer triplet = new StringBuffer();
        StringBuffer proteinSequence = new StringBuffer();
        Object[] checkStop_Sequence = {false, ""};
        for (int i = 1; i <= exonSequenceSegment.length(); i++) {
            triplet.append(exonSequenceSegment.charAt(i - 1));
            if (i % 3 == 0) {
                int indexOfN = triplet.indexOf("N");
                if (indexOfN != -1) triplet.setCharAt(indexOfN, 'A'); // If we get an N. convert to A.
                String decodedValue = codonTable.get(triplet.toString().toUpperCase());
                if (decodedValue == "*") {
                    checkStop_Sequence[0] = true;
                    break;
                } else {
                    checkStop_Sequence[0] = false;
                    triplet = new StringBuffer();
                    proteinSequence.append(decodedValue);
                }
            }
        }
        checkStop_Sequence[1] = proteinSequence.toString();
        return checkStop_Sequence;
    }

    private static String getHeader(String[] annotatedData) {
        int exonCount = Integer.parseInt(annotatedData[8]);
        int name2Index = (2 * exonCount) + 8 + 2; // (+8 to reach the exonCount index. 2 * exonCount gets all the exonstarts and ends +2 index is the name 2 field.)
        return ">" + annotatedData[1] + ":" + annotatedData[name2Index];
    }

    private static String getDataFromDNASequence(String chromozome) throws IOException {
        File file = new File("/Users/sdamaraju/IdeaProjects/BioInformatics/files/hg38.fa");
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line = br.readLine();
        StringBuffer seq = new StringBuffer();
        boolean stop = false;
        while (line != null && !stop) {
            if (line.length() > chromozome.length() && line.substring(1, chromozome.length() + 1).equals(chromozome) && line.length() == (1 + chromozome.length())) {
                line = br.readLine();
                while (line != null && !line.startsWith(">")) {
                    seq.append(line);
                    line = br.readLine();
                }
                stop = true;
            } else line = br.readLine();
        }
        return seq.toString();
    }

    private static String getExonSequenceSegmentFromSequence(int start, int end, String sequence) {
        if (sequence.length() == end) return sequence.substring(start, end);
        return sequence.substring(start, end + 1);
    }

    private static int getFirstNonNegativeExonFrameCounter(int[] exonFrames) {
        for (int i = 0; i < exonFrames.length; i++) {
            if (exonFrames[i] != -1) return i + 1; //(our exonCounter starts from 1)
        }
        return -1;
    }

    private static int getExonStartIndexForPositiveStrand(int exonCounter, String[] lineValues, int exonFrameValue) {
        // to reach the Exon start index: we need to go 8 + exonCounter, reason below..
        // 8 is the index of exonCount and exonStart indices are available after exonCount.
        // exoncounter determines which end index to pick up.
        int startIndex = Integer.parseInt(lineValues[8 + exonCounter]) - exonFrameValue;
        return startIndex;
    }

    private static int getExonEndIndexForPositiveStrand(int exonCounter, int exonCount, String[] lineValues) {
        // to reach the Exon end index: we need to go 8 + exonCount + exonCounter, reason below..
        // 8 is the index of exonCount.
        // exonCount start indices will be there before reaching the exon ends.
        // exoncounter determines which end index to pick up.
        int endIndex = Integer.parseInt(lineValues[8 + exonCount + exonCounter]);
        return endIndex;
    }

    private static int getExonStartIndexForNegativeStrand(int exonCounter, String[] lineValues) {
        int startIndex = Integer.parseInt(lineValues[8 + exonCounter]);
        return startIndex;
    }

    private static int getExonEndIndexForNegativeStrand(int exonCounter, int exonCount, String[] lineValues, int exonFrameValue) {
        int endIndex = Integer.parseInt(lineValues[8 + exonCount + exonCounter]) + exonFrameValue;
        return endIndex - 1;
    }

    private static void buildCodonTable() {
        codonTable.put("TTT", "F");
        codonTable.put("TTC", "F");
        codonTable.put("TTA", "L");
        codonTable.put("TTG", "L");
        codonTable.put("TAT", "Y");
        codonTable.put("TAC", "Y");
        codonTable.put("TAA", "*");
        codonTable.put("TAG", "*");
        codonTable.put("CTT", "L");
        codonTable.put("CTC", "L");
        codonTable.put("CTA", "L");
        codonTable.put("CTG", "L");
        codonTable.put("CAT", "H");
        codonTable.put("CAC", "H");
        codonTable.put("CAA", "Q");
        codonTable.put("CAG", "Q");
        codonTable.put("ATT", "I");
        codonTable.put("ATC", "I");
        codonTable.put("ATA", "I");
        codonTable.put("ATG", "M");
        codonTable.put("AAT", "N");
        codonTable.put("AAC", "N");
        codonTable.put("AAA", "K");
        codonTable.put("AAG", "K");
        codonTable.put("GTT", "V");
        codonTable.put("GTC", "V");
        codonTable.put("GTA", "V");
        codonTable.put("GTG", "V");
        codonTable.put("GAT", "D");
        codonTable.put("GAC", "D");
        codonTable.put("GAA", "E");
        codonTable.put("GAG", "E");
        codonTable.put("TCT", "S");
        codonTable.put("TCC", "S");
        codonTable.put("TCA", "S");
        codonTable.put("TCG", "S");
        codonTable.put("TGT", "C");
        codonTable.put("TGC", "C");
        codonTable.put("TGA", "*");
        codonTable.put("TGG", "W");
        codonTable.put("CCT", "P");
        codonTable.put("CCC", "P");
        codonTable.put("CCA", "P");
        codonTable.put("CCG", "P");
        codonTable.put("CGT", "R");
        codonTable.put("CGC", "R");
        codonTable.put("CGA", "R");
        codonTable.put("CGG", "R");
        codonTable.put("ACT", "T");
        codonTable.put("ACC", "T");
        codonTable.put("ACA", "T");
        codonTable.put("ACG", "T");
        codonTable.put("AGT", "S");
        codonTable.put("AGC", "S");
        codonTable.put("AGA", "R");
        codonTable.put("AGG", "R");
        codonTable.put("GCT", "A");
        codonTable.put("GCC", "A");
        codonTable.put("GCA", "A");
        codonTable.put("GCG", "A");
        codonTable.put("GGT", "G");
        codonTable.put("GGC", "G");
        codonTable.put("GGA", "G");
        codonTable.put("GGG", "G");
    }
}
