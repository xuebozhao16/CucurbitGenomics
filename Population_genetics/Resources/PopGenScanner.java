/**
 * PopGenScanner
 * A population genetics tool for SNP & SV extraction and statistics.
 * Author: Xuebo Zhao (rewritten version)
 *
 * Functions:
 * 1. Read multi-sample VCF (SNP + SV)
 * 2. Compute allele count (AC), allele frequency (AF)
 * 3. Compute missing rate
 * 4. Compute nucleotide diversity (π)
 * 5. Compute per-population allele frequency (if population map provided)
 *
 * Usage:
 *   java PopGenScanner input.vcf.gz population.map output_folder
 */

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class PopGenScanner {

    // sample → population
    Map<String, String> popMap = new HashMap<>();

    // population → sample list
    Map<String, List<Integer>> popIndex = new HashMap<>();

    List<String> samples = new ArrayList<>();

    String vcfFile;
    String popFile;
    String outDir;

    public PopGenScanner(String vcf, String pop, String out) {
        this.vcfFile = vcf;
        this.popFile = pop;
        this.outDir = out;
    }

    /** Load sample→population map */
    void loadPopulationMap() throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(popFile));
        String line;
        while ((line = br.readLine()) != null) {
            if (line.trim().isEmpty()) continue;
            String[] t = line.split("\t");
            popMap.put(t[0], t[1]);   // sample → population
        }
        br.close();
    }

    /** Parse the VCF header to get the list of samples */
    void readSamplesFromVCF(BufferedReader br) throws Exception {
        String line;
        while ((line = br.readLine()) != null) {
            if (line.startsWith("#CHROM")) {
                String[] t = line.split("\t");
                for (int i = 9; i < t.length; i++) samples.add(t[i]);
                break;
            }
        }
        buildPopulationIndex();
    }

    /** Build mapping: population → list of genotype indices */
    void buildPopulationIndex() {
        for (int i = 0; i < samples.size(); i++) {
            String s = samples.get(i);
            String p = popMap.get(s);
            popIndex.computeIfAbsent(p, k -> new ArrayList<>()).add(i);
        }
    }

    /** Convert GT string like "0/1" → {0,1} */
    int[] parseGT(String gt) {
        if (gt.equals("./.")) return null;
        String[] t = gt.split("[/|]");
        return new int[]{ Integer.parseInt(t[0]), Integer.parseInt(t[1]) };
    }

    /** Compute π = 2p(1-p) */
    double nucleotideDiversity(double af) {
        return 2 * af * (1 - af);
    }

    /** Main processing */
    void run() throws Exception {
        loadPopulationMap();

        BufferedReader br = new BufferedReader(new InputStreamReader(
            new GZIPInputStream(new FileInputStream(vcfFile))));

        readSamplesFromVCF(br);

        // Output files
        BufferedWriter bwVar = new BufferedWriter(new FileWriter(outDir + "/variant_stats.tsv"));
        bwVar.write("CHR\tPOS\tREF\tALT\tAC\tAF\tMissingRate\tPi\n");

        BufferedWriter bwPop = new BufferedWriter(new FileWriter(outDir + "/population_AF.tsv"));
        bwPop.write("CHR\tPOS\tALT");
        for (String p : popIndex.keySet()) bwPop.write("\tAF_" + p);
        bwPop.write("\n");

        String line;
        while ((line = br.readLine()) != null) {
            if (line.startsWith("#")) continue;

            String[] t = line.split("\t");
            String chr = t[0];
            int pos = Integer.parseInt(t[1]);
            String ref = t[3];
            String alt = t[4];

            int AC = 0, AN = 0, miss = 0;

            // collect GTs
            List<int[]> genotypes = new ArrayList<>();
            for (int i = 9; i < t.length; i++) {
                String gt = t[i].split(":")[0];
                int[] g = parseGT(gt);
                if (g == null) {
                    miss++;
                    continue;
                }
                genotypes.add(g);
                AC += g[0] + g[1];
                AN += 2;
            }

            double AF = (AN == 0 ? 0 : (double) AC / AN);
            double missingRate = (double) miss / samples.size();
            double pi = nucleotideDiversity(AF);

            // output variant-level statistics
            bwVar.write(chr + "\t" + pos + "\t" + ref + "\t" + alt + "\t" +
                AC + "\t" + AF + "\t" + missingRate + "\t" + pi + "\n");

            // population AF
            bwPop.write(chr + "\t" + pos + "\t" + alt);

            for (String p : popIndex.keySet()) {
                List<Integer> idx = popIndex.get(p);

                int ACp = 0, ANp = 0;
                for (int id : idx) {
                    String gt = t[id + 9].split(":")[0];
                    int[] g = parseGT(gt);
                    if (g == null) continue;
                    ACp += g[0] + g[1];
                    ANp += 2;
                }
                double AFp = (ANp == 0 ? 0 : (double) ACp / ANp);
                bwPop.write("\t" + AFp);
            }
            bwPop.write("\n");
        }

        br.close();
        bwVar.close();
        bwPop.close();

        System.out.println("PopGenScanner finished.");
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 3) {
            System.out.println("Usage: java PopGenScanner input.vcf.gz population.map output_dir");
            System.exit(1);
        }
        PopGenScanner s = new PopGenScanner(args[0], args[1], args[2]);
        s.run();
    }
}
