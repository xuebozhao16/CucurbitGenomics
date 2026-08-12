/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;




/**
 * Base encoder for A(00), C(01), G(10), T(11). Bases are compressed into 2 bits and loaded into long, int, and short primitives.
 * The primitives are padded with ployA (0000...) at the end, if sequence length is less than 32, 16, or 8.
 * @author fl262
 */
public class BaseEncoder {
    public static final int longChunkSize = 32;
    public static final int intChunkSize = 16;
    public static final int shortChunkSize = 8;
    public static final char[] bases = {'A', 'C', 'G', 'T'};
    
    /**
     * Build a byte converter to convert AscII byte following the BaseCoder rules
     * A(00000000), C(00000001), G(00000010), T(0000000011), others(00000100)
     * @return 
     */
    public static HashByteByteMap getAscIIByteMap () {
        HashByteByteMap ascIIByteMap = HashByteByteMaps.newMutableMap();
        for (int i = 0; i < 128; i++) {
            ascIIByteMap.put((byte)i, (byte)4);
        }
        ascIIByteMap.put((byte)65, (byte)0);
        ascIIByteMap.put((byte)97, (byte)0);
        ascIIByteMap.put((byte)67, (byte)1);
        ascIIByteMap.put((byte)99, (byte)1);
        ascIIByteMap.put((byte)71, (byte)2);
        ascIIByteMap.put((byte)103, (byte)2);
        ascIIByteMap.put((byte)84, (byte)3);
        ascIIByteMap.put((byte)116, (byte)3);
        return ascIIByteMap;
    }
    
    /**
     * Return a reverse complementary byte array
     * @param a The byte array must be in 0-3 coding for A, C, G, T, other bases are 4
     * @return 
     */
    public static byte[] getReverseComplementary (byte[] a) {
        byte[] b = new byte[a.length];
        for (int i = 0; i < a.length; i++) {
            int t;
            if (a[i] > 3) {
                t = a[i];
            } 
            else {
                t = (~a[i]<<30)>>>30;
            }
            b[a.length-i-1] = (byte)t;
        }
        return b;
    }
    
    /**
     * Return a long from a byte array. It is padded with polyA at the end if the length of this array is less than 32
     * If the length of this array is greater than 32, return -1 (polyT)
     * @param b The byte array must be in 0-3 coding for A, C, G, T
     * @return 
     */
    public static long getLongSeqFromByteArray(byte[] b) {
        long v = 0;
        if (b.length > longChunkSize) {
            return -1;
        }
        for (int i = 0; i < b.length; i++) {
            v = (v << 2) + b[i];
        }
        v = (v << (2*(longChunkSize-b.length)));
        return v;
    }
    
    /**
     * Return an int from a byte array. It is padded with polyA at the end if the length of this array is less than 16
     * If the length of this array is greater than 16, return -1 (polyT)
     * @param b The byte array must be in 0-3 coding for A, C, G, T
     * @return 
     */
    public static int getIntSeqFromByteArray(byte[] b) {
        int v = 0;
        if (b.length > intChunkSize) {
            return -1;
        }
        for (int i = 0; i < b.length; i++) {
            v = (v << 2) + b[i];
        }
        v = (v << (2*(intChunkSize-b.length)));
        return v;
    }
    
    /**
     * Return an int from a byte array. It is padded with polyA at the end if the length of this array is less than 8
     * If the length of this array is greater than 8, return -1 (polyT)
     * @param b The byte array must be in 0-3 coding for A, C, G, T
     * @return 
     */
    public static short getShortSeqFromByteArray(byte[] b) {
        int v = 0;
        if (b.length > shortChunkSize) {
            return -1;
        }
        for (int i = 0; i < b.length; i++) {
            v = (v << 2) + b[i];
        }
        v = (v << (2*(shortChunkSize-b.length)));
        return (short)v;
    }
    
    /**
     * Return a reverse complementary long. It is padded with polyA at the end if the length of this array is less than 32
     * If the length of this array is greater than 32 or less than 0, return -1 (polyT)
     * @param seq
     * @param kmerLength
     * @return 
     */
    public static long getLongReverseComplement(long seq, int kmerLength) {
        if (kmerLength<0 || kmerLength>longChunkSize) return -1;
        long rev = 0;
        long mask = 3;
        seq = ~seq;
        seq = seq >> ((longChunkSize-kmerLength)*2);
        for (int i = 0; i < kmerLength; i++) {
            rev = (rev << 2) + (seq & mask);
            seq = seq >> 2;
        }
        rev = rev << ((longChunkSize-kmerLength)*2);
        return rev;
    }
    
    /**
     * Return a reverse complementary int. It is padded with polyA at the end if the length of this array is less than 16
     * If the length of this array is greater than 16 or less than 0, return -1 (polyT)
     * @param seq
     * @param kmerLength
     * @return 
     */
    public static int getIntReverseComplement(int seq, int kmerLength) {
        if (kmerLength<0 || kmerLength>intChunkSize) return -1;
        int rev = 0;
        int mask = 3;
        seq = ~seq;
        seq = seq >> ((intChunkSize-kmerLength)*2);
        for (int i = 0; i < kmerLength; i++) {
            rev = (rev << 2) + (seq & mask);
            seq = seq >> 2;
        }
        rev = rev << ((intChunkSize-kmerLength)*2);
        return rev;
    }
    
    /**
     * Return a reverse complementary int. It is padded with polyA at the end if the length of this array is less than 8
     * If the length of this array is greater than 8 or less than 0, return -1 (polyT)
     * @param seq
     * @param kmerLength
     * @return 
     */
    public static short getShortReverseComplement(int seq, int kmerLength) {
        if (kmerLength<0 || kmerLength>shortChunkSize) return -1;
        int rev = 0;
        int mask = 3;
        seq = ~seq;
        seq = seq >> ((shortChunkSize-kmerLength)*2);
        for (int i = 0; i < kmerLength; i++) {
            rev = (rev << 2) + (seq & mask);
            seq = seq >> 2;
        }
        rev = rev << ((shortChunkSize-kmerLength)*2);
        return (short)rev;
    }
    
    /**
     * Return DNA sequence from long
     * @param val
     * @return 
     */
    public static String getSequenceFromLong(long val) {
    	StringBuilder seq = new StringBuilder(longChunkSize + 1);
    	long mask = 3L << 62;
    	for (int i = 0; i < longChunkSize; i++) {
    		byte base = (byte) (((val & mask) >>> 62));
    		seq.append(bases[base]);          
    		val = val << 2;
    	}
    	return seq.toString();
    }
    
    /**
     * Return DNA sequence from long
     * @param val
     * @return 
     */
    public static String getSequenceFromInt(int val) {
    	StringBuilder seq = new StringBuilder(intChunkSize + 1);
    	int mask = 3 << 30;
    	for (int i = 0; i < intChunkSize; i++) {
    		byte base = (byte) (((val & mask) >>> 30));
    		seq.append(bases[base]);          
    		val = val << 2;
    	}
    	return seq.toString();
    }
    
    /**
     * Return DNA sequence from long
     * @param val
     * @return 
     */
    public static String getSequenceFromShort(int val) {
    	StringBuilder seq = new StringBuilder(shortChunkSize + 1);
    	int mask = 3 << 14;
    	for (int i = 0; i < shortChunkSize; i++) {
    		byte base = (byte) (((val & mask) >>> 14));
    		seq.append(bases[base]);          
    		val = val << 2;
    	}
    	return seq.toString();
    }
}
