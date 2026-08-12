/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format.dna;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;

/**
 * Utilities related to DNA sequence
 * @author Fei Lu
 */
public class DNAUtils {
    /**The byte value of 4 DNA bases, A, C, G, T*/
    public static final byte[] baseByte = {65, 67, 71, 84};
    
    private static final byte[] compleBaseByte = {84, 71, 67, 65};
    
    /**
     * The byte value hash map pointing to complementary bases
     */
    public static final HashByteByteMap baseCompleByteMap = 
            HashByteByteMaps.getDefaultFactory().newImmutableMap(baseByte, compleBaseByte);
    
}
