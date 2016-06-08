/*
 *  AlleleDepth
 */
package net.maizegenetics.pal.alignment.depth;

/**
 *
 * @author terry
 */
public interface AlleleDepth {

    /**
     * Returns depth count for each diploid allele at the given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return two counts
     */
    public byte[] getDepthForAlleles(int taxon, int site);
}
