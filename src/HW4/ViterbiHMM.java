package HW4;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class ViterbiHMM {
    static String sequence = "";
    static List match = new ArrayList();
    static List insert = new ArrayList();
    static List state = new ArrayList();
    static List runToUse = new ArrayList();
    static List uniqueChars = new ArrayList();

    public static class CustomObj {
        Double match;
        Double insert;
        Double delete;

        @Override
        public String toString() {
            return "\n match : " + match + " \n insert : " + insert + " \n delete : " + delete;
        }
    }

    public static void main(String[] args) throws IOException {
        File file1 = new File(args[0]);
        File file2 = new File(args[1]);
        Boolean showEachDetail = Boolean.parseBoolean(args[2]);
        BufferedReader br1 = new BufferedReader(new FileReader(file1));
        BufferedReader br2 = new BufferedReader(new FileReader(file2));
        parseHMMFile(br1);
        parseFastaFile(br2);
        List sequenceAsList = new ArrayList();
        for (char c : sequence.toCharArray()) {
            sequenceAsList.add(c);
        }
        CustomObj obj = new CustomObj();
        obj.match = 1.0;
        obj.insert = 0.0;
        obj.delete = 0.0;
        List transitions = evaluateStates(state);
        List emissions = evaluateEmissions();
        applyViterbi(sequenceAsList, obj, transitions, emissions, showEachDetail);

    }

    public static void applyViterbi(List sequence, CustomObj start, List transitions, List emissions, boolean showEachDetail) {
        Double score = 0.0;
        List allStates = new ArrayList();
        List<Character> defaultOP = new ArrayList();
        for (int i = 0; i < sequence.size(); i++)
            defaultOP.add('*');
        Double sMatch = start.match;
        Double sInsert = start.insert;
        Double sDelete = start.delete;
        System.out.println("\nAll scores at individual sequence bit for comparison.. and corresponding state manipulations..");
        for (int i = 0; i < sequence.size() - 1; i++) {
            Character in = (Character) sequence.get(i);
            Double emissionMatch = (Double) ((HashMap) ((List) emissions.get(i)).get(0)).get(in.toString());
            Double emissionInsert = (Double) ((HashMap) ((List) emissions.get(i)).get(1)).get(in.toString());
            int index = 0;
            for (int k = 0; k < uniqueChars.size(); k++) {
                if (((String) uniqueChars.get(k)).charAt(0) == (in)) {
                    index = k;
                    break;
                }
            }
            Double transform = transform((String) runToUse.get(index));
            if (emissionMatch > transform) defaultOP.set(i, '+'); // emission probability is greater than background frequency here...
            Double firstMatch = Math.log(emissionMatch / transform);
            Double firstInsert = Math.log(emissionInsert / transform);
            CustomObj match_ = (CustomObj) ((List) transitions.get(i)).get(0);
            CustomObj insert_ = (CustomObj) ((List) transitions.get(i)).get(1);
            CustomObj delete_ = (CustomObj) ((List) transitions.get(i)).get(2);
            Double transitionMM = match_.match;
            Double transitionIM = insert_.match;
            Double transitionDM = delete_.match;
            Double transitionMI = match_.insert;
            Double transitionII = insert_.insert;
            Double transitionMD = match_.delete;
            Double transitionDD = delete_.delete;
            Double prevStartMatch = sMatch;
            Double prevStartInsert = sInsert;
            Double prevStartDelete = sDelete;

            if (transitionMM != 0.0 && transitionIM != 0.0 && transitionDM != 0.0) {
                sMatch = firstMatch + Math.max(prevStartDelete + Math.log(transitionDM), Math.max(prevStartMatch + Math.log(transitionMM), prevStartInsert + Math.log(transitionIM)));
            } else {
                sMatch = prevStartMatch;
            }

            if (transitionMI != 0.0 && transitionII != 0.0) {
                sInsert = firstInsert + Math.max(prevStartMatch + Math.log(transitionMI), prevStartInsert + Math.log(transitionII));
            } else {
                sInsert = prevStartInsert;
            }

            if (transitionMD != 0 && transitionDD != 0) {
                sDelete = Math.max(prevStartMatch + Math.log(transitionMD), prevStartDelete + Math.log(transitionDD));
            } else {
                sDelete = prevStartInsert;
            }
            if(showEachDetail) {
                System.out.println("\n {delete score " + sDelete);
                System.out.print(" match score " + sMatch + "\n");
                System.out.print(" insert score " + sInsert + "}\n");
            }

            score = Math.max(sDelete,Math.max(sMatch,sInsert));
            System.out.println(" max is : " + score);

            if(sMatch.equals(score)) allStates.add("m");
            else if(sInsert.equals(score)) allStates.add("i");
            else allStates.add("d");

            System.out.println("Inserting "+allStates.get(allStates.size()-1)+" into states.");
        }
        System.out.println();
        System.out.println("Likelihood score calculated from above is : "+score);
        System.out.println("List of all States : " +allStates);
        System.out.println("Sequence Identified from highest emission probability : " );
        System.out.println(defaultOP);
    }


    public static List evaluateEmissions() {
        List emissions = new ArrayList();
        for (int i = 0; i < match.size(); i++) {
            HashMap match_ = new HashMap();
            HashMap insert_ = new HashMap();
            for (int j = 0; j < uniqueChars.size(); j++) {
                match_.put(uniqueChars.get(j), transform((String) ((ArrayList) (match.get(i))).get(j)));
                insert_.put(uniqueChars.get(j), transform((String) ((ArrayList) (insert.get(i))).get(j)));
            }
            List emission = new ArrayList();
            emission.add(match_);
            emission.add(insert_);
            emissions.add(emission);
        }
        return emissions;
    }

    public static List evaluateStates(List states) {
        List allTransitions = new ArrayList();
        for (int i = 0; i < states.size(); i++) {
            List state = (ArrayList) states.get(0);
            CustomObj match = new CustomObj();
            match.match = transform((String) state.get(0));
            match.insert = transform((String) state.get(1));
            match.delete = transform((String) state.get(2));
            CustomObj insert = new CustomObj();
            insert.match = transform((String) state.get(3));
            insert.insert = transform((String) state.get(4));
            CustomObj delete = new CustomObj();
            delete.match = transform((String) state.get(5));
            delete.delete = transform((String) state.get(6));

            List op = new ArrayList();
            op.add(match);
            op.add(insert);
            op.add(delete);
            allTransitions.add(op);
        }
        return allTransitions;
    }

    static Double transform(String value) {
        if (value == "*") return 0.0;
        return Math.exp(-(Double.parseDouble(value)));
    }


    public static void parseFastaFile(BufferedReader br) throws IOException {
        String line = br.readLine();
        while (line != null) {
            line = line.stripLeading();
            if (line.startsWith(">")) {
                line = br.readLine();
                continue;
            } else {
                sequence = sequence + line;
            }
            line = br.readLine();
        }
        System.out.println("Sequence from Fasta format as follows:");
        System.out.println(sequence);
    }

    public static void parseHMMFile(BufferedReader br) throws IOException {
        String line = br.readLine();
        List table = new ArrayList();
        List transitions = new ArrayList();
        int start = 0;
        List<String> lines = new ArrayList();
        while (line != null) {
            lines.add(line);
            line = br.readLine();
        }
        for (int i = 0; i < lines.size(); i++) {
            start = 0;
            if (lines.get(i).startsWith("HMM ")) {
                start = i;
                break;
            }
        }
        for (int i = start; i < lines.size(); i++) {
            table.add(lines.get(i));
        }
        String temp = (String) table.get(0);
        uniqueChars = new ArrayList<>(Arrays.asList(((String) table.get(0)).split("\\s+")));
        uniqueChars.remove(0);
        transitions = new ArrayList<>(Arrays.asList(((String) table.get(1)).split("\\s+")));
        transitions.remove(0);
        for (int i = 0; i < 2; i++)
            table.remove(0);
        runToUse = new ArrayList(Arrays.asList(((String) table.get(0)).split("\\s+")));
        runToUse.remove(0);
        runToUse.remove(0);
        for (int i = 0; i < 3; i++)
            table.remove(0);
        table.remove(table.size() - 1);

        for (int j = 0; j < table.size() / 3; j++) {
            List l1 = new ArrayList(Arrays.asList(((String) table.get(j * 3)).split("\\s+")));
            l1.remove(0);
            l1.remove(0);
            l1.remove(l1.size() - 1);
            l1.remove(l1.size() - 1);
            l1.remove(l1.size() - 1);
            l1.remove(l1.size() - 1);
            match.add(l1);
            List l2 = new ArrayList(Arrays.asList((((String) table.get(j * 3 + 1)).split("\\s+"))));
            l2.remove(0);
            insert.add(l2);
            List l3 = new ArrayList(Arrays.asList((((String) table.get(j * 3 + 2)).split("\\s+"))));
            l3.remove(0);
            state.add(l3);
        }
    }


}
