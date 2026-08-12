/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format.dna;

/**
 * The {@code SequenceByte} use one byte to store a DNA base. 
 * <p>
 * It supports standard IUPAC DNA coding (https://www.bioinformatics.org/sms/iupac.html). Non-"ATGCN-." are ignored for now.
 * @author Fei Lu
 */
public class SequenceByte implements SequenceInterface {
    
    byte[] seqByte = null;
    
    /**
     * Constructs a {@code SequenceByte} from {@code String}. The lower case bases are converted to upper case.
     * @param seq 
     *        The name of a supported DNA sequence
     */
    public SequenceByte (String seq) {
        seqByte = seq.toUpperCase().getBytes();
    }
    
    @Override
    public int getSequenceLength() {
        return seqByte.length;
    }

    @Override
    public double getProportionA() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 65) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getProportionT() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 84) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getProportionG() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 71) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getProportionC() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 67) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getGCContent() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 67 || this.seqByte[i] == 84) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public char getBase(int positionIndex) {
        return (char)this.seqByte[positionIndex];
    }

    @Override
    public String getSequence() {
        return new String(this.seqByte);
    }

    @Override
    public String getSequence(int startIndex, int endIndex) {
        return new String(this.seqByte, startIndex, endIndex-startIndex);
    }

    @Override
    public String getReverseComplementarySeq() {
        return this.getReverseComplementarySeq(0, this.getSequenceLength());
    }

    @Override
    public String getReverseComplementarySeq(int startIndex, int endIndex) {
        byte[] reverseByte = new byte[endIndex - startIndex];
        for (int i = 0; i < reverseByte.length; i++) {
            reverseByte[i] = DNAUtils.baseCompleByteMap.get(seqByte[endIndex-i-1]);
        }
        return new String(reverseByte);
    }
    
    /**
     * Return if the sequence has gaps, "-" or "."
     * @return 
     */
    public boolean isThereGap () {
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 45 || this.seqByte[i] == 46) return true;
        }
        return false;
    }
    
    /**
     * Return if the sequence has N
     * @return 
     */
    public boolean isThereN () {
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 78) return true;
        }
        return false;
    }
}
