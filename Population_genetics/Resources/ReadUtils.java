/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package format;

import java.util.HashMap;

/**
 *
 * @author Fei Lu
 */
public class ReadUtils {
    /**Illumina TruSeq universal adaptor*/
    public static String truSeqUAdaptor = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    /**Illumina TruSeq index adaptor part 1, before barcode*/
    public static String truSeqIAdaptor1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
    /**Illumina TruSeq index adaptor part 2, after barcode*/
    public static String truSeqIAdaptor2 = "ATCTCGTATGCCGTCTTCTGCTTG";
    
    public static final HashMap<String, String> baseCompleMap = new HashMap();
    static {
        baseCompleMap.put("A", "T");
        baseCompleMap.put("T", "A");
        baseCompleMap.put("G", "C");
        baseCompleMap.put("C", "G");
        baseCompleMap.put("N", "N");
        baseCompleMap.put("a", "t");
        baseCompleMap.put("t", "a");
        baseCompleMap.put("g", "c");
        baseCompleMap.put("c", "g");
        baseCompleMap.put("n", "n");
    }
}
