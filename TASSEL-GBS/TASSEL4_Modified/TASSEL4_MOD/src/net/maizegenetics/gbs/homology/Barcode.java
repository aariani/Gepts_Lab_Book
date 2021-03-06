package net.maizegenetics.gbs.homology;

import java.util.Arrays;

import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Container class for storing information on GBS barcodes.
 * 
 * @author Ed Buckler
 */
public class Barcode implements Comparable<Barcode> {
    /**Barcode sequence */
    String barcodeS;
    /**Overhang sequence from the restriction enzyme */
    String[] overhangS;
    /**Taxon (sample) name */
    String taxaName;
    /**Flowcell name */
    String flowcell;
     /**Flowcell lane name */
    String lane;
    /**Barcode  encoded in 2-bit long*/
    long[] barOverLong;
    /**Length of barcode plus overhang*/
    int barOverLength;
     /**Length of barcode*/
    int barLength;

    /**
     * Constructor creating a barcode
     * @param barcodeS barcode sequence
     * @param overhangSunsort overhang sequence array unsorted (array size is 1 for
     * non-degerate restriction enzymes, degenerate enzymes >1)
     * @param taxa name of taxon (sample) 
     * @param flowcell name of the flowcell
     * @param lane name of the lane
     */
    public Barcode(String barcodeS, String[] overhangSunsort, String taxa, String flowcell, String lane) {
        this.barcodeS = barcodeS;
        Arrays.sort(overhangSunsort);
        this.overhangS = overhangSunsort;
        this.flowcell = flowcell;
        this.lane = lane;
        this.taxaName = taxa;
        barOverLong = new long[overhangS.length];
        for (int i = 0; i < overhangS.length; i++) {
            barOverLong[i] = BaseEncoder.getLongFromSeq(barcodeS + overhangS[i]);
        }
        barOverLength = barcodeS.length() + overhangS[0].length();
        barLength = barcodeS.length();
    }

    /**
     * Return the minimum sequence divergence between a query sequence and 
     * barcode with its overhang 
     * @param queryLong query sequence encoded in 2-bit long
     * @param maxDivCheck maximum divergence to search upto
     * @return minimum divergence between barcode and query
     */
    public int compareSequence(long queryLong, int maxDivCheck) {
        int div = barOverLength;
        for (long targetLong : barOverLong) {
            int c = BaseEncoder.seqDifferencesForSubset(targetLong, queryLong, barOverLength, maxDivCheck);
            if (c < div) {
                div = c;
            }
        }
        return div;
    }

    @Override
    public int compareTo(Barcode anotherBarcode) {
        if (this.barOverLong[0] < anotherBarcode.barOverLong[0]) {
            return -1;
        }
        if (this.barOverLong[0] > anotherBarcode.barOverLong[0]) {
            return 1;
        }
        return 0;
    }

    public String getTaxaName() {
        return taxaName;
    }
}
