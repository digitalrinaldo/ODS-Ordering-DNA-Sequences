import java.util.*;

/**
 * ODS — Ordering DNA Sequences by Minimizing Manhattan (L1) Distance via Simulated Annealing
 * Author: Rinaldo Di Giorgio (DiGiorgio Labs — https://digiorgiolabs.com)
 * License: Public Domain (The Unlicense) / MIT
 */
public class OrderingDnaSequences {

    public static int manhattanDistance(int[] a, int[] b) {
        int sum = 0;
        for (int i = 0; i < a.length; i++) {
            sum += Math.abs(a[i] - b[i]);
        }
        return sum;
    }

    public static int totalOrderingCost(List<int[]> sequences, int[] order) {
        int cost = 0;
        for (int i = 0; i < order.length - 1; i++) {
            cost += manhattanDistance(sequences.get(order[i]), sequences.get(order[i + 1]));
        }
        return cost;
    }

    public static int[] encodeDnaKmerCounts(String dna, int k) {
        int dim = 1 << (2 * k);
        int[] counts = new int[dim];
        for (int i = 0; i <= dna.length() - k; i++) {
            int idx = 0;
            for (int j = 0; j < k; j++) {
                char c = Character.toUpperCase(dna.charAt(i + j));
                int v = (c == 'A') ? 0 : (c == 'C') ? 1 : (c == 'G') ? 2 : 3;
                idx = (idx << 2) | v;
            }
            counts[idx]++;
        }
        return counts;
    }

    public static int[] simulatedAnnealingOrder(List<int[]> sequences, double initialTemp, double coolingRate, int iterations, long seed) {
        int n = sequences.size();
        int[] currentOrder = new int[n];
        for (int i = 0; i < n; i++) currentOrder[i] = i;
        int[] bestOrder = currentOrder.clone();

        int currentCost = totalOrderingCost(sequences, currentOrder);
        int bestCost = currentCost;

        Random rng = new Random(seed);
        double temp = initialTemp;

        for (int step = 0; step < iterations && temp > 1e-6; step++) {
            int i = rng.nextInt(n);
            int j = rng.nextInt(n);
            if (i == j) continue;

            int tmp = currentOrder[i];
            currentOrder[i] = currentOrder[j];
            currentOrder[j] = tmp;

            int candidateCost = totalOrderingCost(sequences, currentOrder);
            int delta = candidateCost - currentCost;

            if (delta < 0 || Math.exp(-delta / temp) > rng.nextDouble()) {
                currentCost = candidateCost;
                if (currentCost < bestCost) {
                    bestCost = currentCost;
                    bestOrder = currentOrder.clone();
                }
            } else {
                currentOrder[j] = currentOrder[i];
                currentOrder[i] = tmp;
            }
            temp *= coolingRate;
        }
        return bestOrder;
    }

    public static void main(String[] args) {
        List<String> rawSequences = Arrays.asList(
            "ACGTACGTACGTACGT",
            "TTTTGGGGCCCCAAAA",
            "ACGTACGTTTTTGGGG",
            "TTTTGGGGAAAACCCC",
            "ACGTACGTACGTGGGG"
        );
        List<int[]> vectors = new ArrayList<>();
        for (String seq : rawSequences) {
            vectors.add(encodeDnaKmerCounts(seq, 2));
        }

        int[] initial = {0, 1, 2, 3, 4};
        System.out.println("Initial Manhattan (L1) path cost: " + totalOrderingCost(vectors, initial));

        int[] optimal = simulatedAnnealingOrder(vectors, 100.0, 0.995, 10000, 42L);
        System.out.println("Optimized Manhattan (L1) path cost: " + totalOrderingCost(vectors, optimal));
        for (int idx : optimal) {
            System.out.println("  [" + idx + "] " + rawSequences.get(idx));
        }
    }
}
