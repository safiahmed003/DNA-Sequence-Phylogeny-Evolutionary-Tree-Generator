
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.Set;

/** Console entry point: DNA Sequence Phylogeny & Evolutionary Tree Generator. */
public class DNAPhylogeny {
    private static final Scanner in = new Scanner(System.in);
    private static final String LINE = "========================================================";
    private static final String DASH = "--------------------------------------------------------";

    public static void main(String[] args) {
        // Optional: java -cp out phylogeny.Main data/sample_primates.txt  -> run directly from a file
        if (args.length == 1) {
            try {
                run(DatasetLoader.load(args[0]));
            } catch (IOException e) {
                System.out.println("Error: " + e.getMessage());
            }
            return;
        }

        try {
            while (true) {
                System.out.println("\n" + LINE);
                System.out.println("   DNA PHYLOGENY & EVOLUTIONARY TREE GENERATOR");
                System.out.println("   Edit Distance (DP) + Kruskal's MST + DSU");
                System.out.println(LINE);
                System.out.println("  1. Run built-in sample (4 primates)");
                System.out.println("  2. Enter species manually");
                System.out.println("  3. Load species from a file");
                System.out.println("  0. Exit");
                System.out.print("Choose an option: ");

                String choice = in.nextLine().trim();
                switch (choice) {
                    case "1": run(sampleData()); break;
                    case "2": run(manualInput()); break;
                    case "3": loadFromFile(); break;
                    case "0": System.out.println("Goodbye!"); return;
                    default:  System.out.println("Invalid option. Please enter 0, 1, 2 or 3.");
                }
            }
        } catch (NoSuchElementException e) {
            System.out.println("\nInput ended. Goodbye!");
        }
    }

    private static List<Species> sampleData() {
        List<Species> l = new ArrayList<>();
        l.add(new Species("Human",      "AGCTTAGCTA"));
        l.add(new Species("Chimpanzee", "AGCTTAGTTA"));
        l.add(new Species("Gorilla",    "AACTTAGCTA"));
        l.add(new Species("Orangutan",  "TGCTTAGCTG"));
        return l;
    }

    private static List<Species> manualInput() {
        int n = 0;
        while (n < 2) {
            System.out.print("Enter number of species (>= 2): ");
            try {
                n = Integer.parseInt(in.nextLine().trim());
                if (n < 2) System.out.println("  Need at least 2 species.");
            } catch (NumberFormatException e) {
                System.out.println("  Please enter a valid integer.");
            }
        }
        List<Species> list = new ArrayList<>();
        Set<String> names = new HashSet<>();
        System.out.println("\nEnter species names and their DNA sequences (A, C, G, T only):");
        for (int i = 1; i <= n; i++) {
            String name;
            while (true) {
                System.out.print("Species " + i + " name: ");
                name = in.nextLine().trim();
                if (name.isEmpty()) System.out.println("  Name cannot be empty.");
                else if (names.contains(name.toLowerCase())) System.out.println("  Name already used.");
                else break;
            }
            names.add(name.toLowerCase());
            String seq;
            while (true) {
                System.out.print("Species " + i + " DNA sequence (e.g., AGCT): ");
                seq = DnaValidator.normalize(in.nextLine());
                if (DnaValidator.isValidDna(seq)) break;
                System.out.println("  Invalid! Use only A, C, G, T (no spaces, not empty).");
            }
            list.add(new Species(name, seq));
        }
        return list;
    }

    private static void loadFromFile() {
        System.out.print("Enter file path (e.g., data/sample_primates.txt): ");
        String path = in.nextLine().trim();
        try {
            run(DatasetLoader.load(path));
        } catch (IOException e) {
            System.out.println("Error: " + e.getMessage());
        }
    }

    /** Builds the tree, prints the full report, and saves it to output/result.txt. */
    private static void run(List<Species> species) {
        PhylogenyBuilder b = new PhylogenyBuilder(species);
        long t0 = System.nanoTime();
        b.build();
        long micros = (System.nanoTime() - t0) / 1000;

        StringBuilder sb = new StringBuilder();
        sb.append("\n").append(LINE).append("\n  INPUT SPECIES\n").append(LINE).append("\n");
        for (int i = 0; i < species.size(); i++)
            sb.append(String.format("%2d. %-18s %s%n", i + 1, species.get(i).getName(), species.get(i).getSequence()));

        sb.append("\n").append(LINE).append("\n  PAIRWISE GENETIC DISTANCE MATRIX (Edit Distance)\n").append(LINE).append("\n");
        sb.append(TreePrinter.distanceMatrix(b));

        sb.append("\n").append(LINE).append("\n  KRUSKAL'S ALGORITHM - STEP BY STEP\n").append(LINE).append("\n");
        int step = 1;
        for (String s : b.getKruskalLog()) sb.append(String.format("Step %2d: %s%n", step++, s));

        sb.append("\n").append(LINE).append("\n  EVOLUTIONARY TREE (MST OUTPUT)\n").append(LINE).append("\n");
        sb.append(TreePrinter.mstEdges(b));
        sb.append(DASH).append("\n");
        sb.append("Total Evolutionary Distance Cost: ").append(b.getTotalCost()).append("\n");

        sb.append("\n").append(LINE).append("\n  TREE VIEW (rooted at ").append(species.get(0).getName())
          .append("; number = distance from parent)\n").append(LINE).append("\n");
        sb.append(TreePrinter.asciiTree(b));

        Edge closest = b.getAllEdges().get(0);
        Edge farthest = b.getAllEdges().get(b.getAllEdges().size() - 1);
        sb.append("\n").append(LINE).append("\n  SUMMARY\n").append(LINE).append("\n");
        sb.append("Closest pair  : ").append(species.get(closest.u).getName()).append(" & ")
          .append(species.get(closest.v).getName()).append(" (distance ").append(closest.weight).append(")\n");
        sb.append("Farthest pair : ").append(species.get(farthest.u).getName()).append(" & ")
          .append(species.get(farthest.v).getName()).append(" (distance ").append(farthest.weight).append(")\n");
        sb.append("Species: ").append(species.size()).append(" | Pairs compared: ").append(b.getAllEdges().size())
          .append(" | MST edges: ").append(b.getMstEdges().size()).append("\n");
        sb.append("Computation time: ").append(micros).append(" microseconds\n");

        System.out.print(sb);
        saveReport(sb.toString());
    }

    private static void saveReport(String text) {
        try {
            File dir = new File("output");
            if (!dir.exists()) dir.mkdirs();
            try (PrintWriter pw = new PrintWriter(new FileWriter(new File(dir, "result.txt")))) {
                pw.print(text);
            }
            System.out.println("\n[Saved] Report written to output/result.txt");
        } catch (IOException e) {
            System.out.println("\n[Warning] Could not save report: " + e.getMessage());
        }
    }

    // ======================================================================
    //  HELPER CLASSES (kept in one file for easy sharing)
    // ======================================================================

    /** Represents one organism: its name and DNA sequence. */
    static class Species {
        private final String name;
        private final String sequence;

        public Species(String name, String sequence) {
            this.name = name;
            this.sequence = sequence;
        }

        public String getName() { return name; }
        public String getSequence() { return sequence; }
        public int length() { return sequence.length(); }

        @Override
        public String toString() { return name + " (" + sequence.length() + " bp)"; }
    }

    /** Weighted undirected edge between species u and v (weight = genetic distance). */
    static class Edge implements Comparable<Edge> {
        public final int u, v, weight;

        public Edge(int u, int v, int weight) {
            this.u = u;
            this.v = v;
            this.weight = weight;
        }

        /** Sort by weight; ties broken by (u, v) so the output is always deterministic. */
        @Override
        public int compareTo(Edge o) {
            if (weight != o.weight) return Integer.compare(weight, o.weight);
            if (u != o.u) return Integer.compare(u, o.u);
            return Integer.compare(v, o.v);
        }
    }

    /** Disjoint Set Union (Union-Find) with path compression and union by rank. */
    static class DisjointSet {
        private final int[] parent;
        private final int[] rank;

        public DisjointSet(int n) {
            parent = new int[n];
            rank = new int[n];
            for (int i = 0; i < n; i++) parent[i] = i;
        }

        public int find(int x) {
            if (parent[x] != x) parent[x] = find(parent[x]); // path compression
            return parent[x];
        }

        /** @return true if x and y were in different sets and have now been merged. */
        public boolean union(int x, int y) {
            int rx = find(x), ry = find(y);
            if (rx == ry) return false;
            if (rank[rx] < rank[ry]) parent[rx] = ry;
            else if (rank[rx] > rank[ry]) parent[ry] = rx;
            else { parent[ry] = rx; rank[rx]++; }
            return true;
        }
    }

    /** Levenshtein (edit) distance using bottom-up Dynamic Programming. */
    static final class EditDistance {
        private EditDistance() {}

        /**
         * Minimum number of insertions, deletions and substitutions to turn a into b.
         * Time: O(m*n). Space: O(m*n).
         */
        public static int compute(String a, String b) {
            int m = a.length(), n = b.length();
            int[][] dp = new int[m + 1][n + 1];

            for (int i = 0; i <= m; i++) dp[i][0] = i;
            for (int j = 0; j <= n; j++) dp[0][j] = j;

            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    if (a.charAt(i - 1) == b.charAt(j - 1)) {
                        dp[i][j] = dp[i - 1][j - 1];
                    } else {
                        dp[i][j] = 1 + Math.min(dp[i - 1][j - 1],            // substitution
                                       Math.min(dp[i - 1][j], dp[i][j - 1])); // deletion / insertion
                    }
                }
            }
            return dp[m][n];
        }

        /** Similarity in percent: 100 * (1 - distance / longerLength). */
        public static double similarityPercent(String a, String b, int distance) {
            int max = Math.max(a.length(), b.length());
            if (max == 0) return 100.0;
            return 100.0 * (1.0 - (double) distance / max);
        }
    }

    /** Input validation helpers. */
    static final class DnaValidator {
        private DnaValidator() {}

        /** Valid DNA = non-empty string containing only A, C, G, T (case-insensitive). */
        public static boolean isValidDna(String s) {
            if (s == null || s.isEmpty()) return false;
            for (int i = 0; i < s.length(); i++) {
                char c = Character.toUpperCase(s.charAt(i));
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
            }
            return true;
        }

        public static String normalize(String s) {
            return s.trim().toUpperCase();
        }
    }

    /**
     * Builds the evolutionary tree:
     *   1) pairwise edit-distance matrix
     *   2) complete weighted graph (edge list)
     *   3) Kruskal's MST using DSU
     */
    static class PhylogenyBuilder {
        private final List<Species> species;
        private final int n;
        private int[][] distanceMatrix;
        private List<Edge> allEdges;
        private List<Edge> mstEdges;
        private List<String> kruskalLog;   // step-by-step trace (accepted / rejected)
        private int totalCost;

        public PhylogenyBuilder(List<Species> species) {
            this.species = species;
            this.n = species.size();
        }

        public void build() {
            computeDistances();
            runKruskal();
        }

        private void computeDistances() {
            distanceMatrix = new int[n][n];
            allEdges = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    int d = EditDistance.compute(species.get(i).getSequence(), species.get(j).getSequence());
                    distanceMatrix[i][j] = d;
                    distanceMatrix[j][i] = d;
                    allEdges.add(new Edge(i, j, d));
                }
            }
            Collections.sort(allEdges);
        }

        private void runKruskal() {
            DisjointSet dsu = new DisjointSet(n);
            mstEdges = new ArrayList<>();
            kruskalLog = new ArrayList<>();
            totalCost = 0;

            for (Edge e : allEdges) {
                String label = species.get(e.u).getName() + " - " + species.get(e.v).getName()
                             + " (d=" + e.weight + ")";
                if (dsu.union(e.u, e.v)) {
                    mstEdges.add(e);
                    totalCost += e.weight;
                    kruskalLog.add("ACCEPT  " + label);
                    if (mstEdges.size() == n - 1) break;
                } else {
                    kruskalLog.add("REJECT  " + label + "  [would form a cycle]");
                }
            }
        }

        public List<Species> getSpecies() { return species; }
        public int[][] getDistanceMatrix() { return distanceMatrix; }
        public List<Edge> getAllEdges() { return allEdges; }
        public List<Edge> getMstEdges() { return mstEdges; }
        public List<String> getKruskalLog() { return kruskalLog; }
        public int getTotalCost() { return totalCost; }
    }

    /** Formats matrices, MST edges and an ASCII tree (pure ASCII so it works on any console). */
    static final class TreePrinter {
        private TreePrinter() {}

        public static String distanceMatrix(PhylogenyBuilder b) {
            List<Species> sp = b.getSpecies();
            int n = sp.size();
            int w = 6;
            for (Species s : sp) w = Math.max(w, s.getName().length() + 2);

            StringBuilder sb = new StringBuilder();
            sb.append(pad("", w));
            for (Species s : sp) sb.append(pad(s.getName(), w));
            sb.append('\n');
            for (int i = 0; i < n; i++) {
                sb.append(pad(sp.get(i).getName(), w));
                for (int j = 0; j < n; j++) sb.append(pad(String.valueOf(b.getDistanceMatrix()[i][j]), w));
                sb.append('\n');
            }
            return sb.toString();
        }

        public static String mstEdges(PhylogenyBuilder b) {
            StringBuilder sb = new StringBuilder();
            List<Species> sp = b.getSpecies();
            int k = 1;
            for (Edge e : b.getMstEdges()) {
                Species a = sp.get(e.u), c = sp.get(e.v);
                double sim = EditDistance.similarityPercent(a.getSequence(), c.getSequence(), e.weight);
                sb.append(String.format("%2d) %s <---> %s | Genetic Distance: %d | Similarity: %.1f%%%n",
                        k++, a.getName(), c.getName(), e.weight, sim));
            }
            return sb.toString();
        }

        /** ASCII tree rooted at species 0. */
        public static String asciiTree(PhylogenyBuilder b) {
            int n = b.getSpecies().size();
            List<List<Edge>> adj = new ArrayList<>();
            for (int i = 0; i < n; i++) adj.add(new ArrayList<>());
            for (Edge e : b.getMstEdges()) {
                adj.get(e.u).add(e);
                adj.get(e.v).add(e);
            }
            StringBuilder sb = new StringBuilder();
            sb.append(b.getSpecies().get(0).getName()).append("  [root]\n");
            draw(0, -1, "", adj, b, sb);
            return sb.toString();
        }

        private static void draw(int node, int parent, String prefix,
                                 List<List<Edge>> adj, PhylogenyBuilder b, StringBuilder sb) {
            List<Edge> kids = new ArrayList<>();
            for (Edge e : adj.get(node)) {
                int other = (e.u == node) ? e.v : e.u;
                if (other != parent) kids.add(e);
            }
            for (int i = 0; i < kids.size(); i++) {
                Edge e = kids.get(i);
                int child = (e.u == node) ? e.v : e.u;
                boolean last = (i == kids.size() - 1);
                sb.append(prefix).append(last ? "\\-- " : "+-- ")
                  .append(b.getSpecies().get(child).getName())
                  .append(" (").append(e.weight).append(")\n");
                draw(child, node, prefix + (last ? "    " : "|   "), adj, b, sb);
            }
        }

        private static String pad(String s, int w) {
            StringBuilder sb = new StringBuilder(s);
            while (sb.length() < w) sb.append(' ');
            return sb.toString();
        }
    }

    /**
     * Loads species from a text file. Each non-empty line:  Name SEQUENCE
     * Lines starting with '#' are comments. The LAST token on a line is the DNA
     * sequence, everything before it is the name (so names may contain spaces).
     */
    static final class DatasetLoader {
        private DatasetLoader() {}

        public static List<Species> load(String path) throws IOException {
            List<Species> list = new ArrayList<>();
            Set<String> seen = new HashSet<>();
            try (BufferedReader br = new BufferedReader(new FileReader(path))) {
                String line;
                int lineNo = 0;
                while ((line = br.readLine()) != null) {
                    lineNo++;
                    line = line.trim();
                    if (line.isEmpty() || line.startsWith("#")) continue;
                    int idx = Math.max(line.lastIndexOf(' '), line.lastIndexOf('\t'));
                    if (idx < 0) throw new IOException("Line " + lineNo + ": expected 'Name SEQUENCE'");
                    String name = line.substring(0, idx).trim();
                    String seq = DnaValidator.normalize(line.substring(idx + 1));
                    if (name.isEmpty()) throw new IOException("Line " + lineNo + ": missing species name");
                    if (!DnaValidator.isValidDna(seq))
                        throw new IOException("Line " + lineNo + ": invalid DNA (only A,C,G,T allowed): " + seq);
                    if (!seen.add(name.toLowerCase()))
                        throw new IOException("Line " + lineNo + ": duplicate species name: " + name);
                    list.add(new Species(name, seq));
                }
            }
            if (list.size() < 2) throw new IOException("Need at least 2 species in the file.");
            return list;
        }
    }
}
