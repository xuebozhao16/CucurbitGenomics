/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package utils;

import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.util.Arrays;

/**
 *
 * @author fl262
 */
public class StringArrayUtils {
    
    public StringArrayUtils() {}
    /**
     * Sort string array based on its converted int value
     * String array should be convertable
     * @return 
     */
    public String[] sortStringArrayByIntValue (String[] a) {
        System.out.println("Start sorting String array based on values");
        int[] value = new int[a.length];
        StringAndIntValue[] sv = new StringAndIntValue[a.length];
        for (int i = 0; i < value.length; i++) {
            sv[i] = new StringAndIntValue(a[i], Integer.parseInt(a[i]));
        }
        Arrays.sort(sv);
        for (int i = 0; i < a.length; i++) {
            a[i] = sv[i].s;
        }
        System.out.println("Finished sort");
        return a;
    }
    
    private class StringAndIntValue implements Comparable<StringAndIntValue> {
        String s;
        int value;
        StringAndIntValue (String s, int value) {
            this.s = s;
            this.value = value;
        }
        @Override
        public int compareTo(StringAndIntValue o) {
            return value - o.value;
        }
    }
    
    
}
