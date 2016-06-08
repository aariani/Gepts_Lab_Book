/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Fei Lu
 */
public class TagAgainstAnchorOriginal {
    BitAlignment anchor;
    double[] anchorMaf;
    int[] chromosomeNumber;
    int[] chrStartIndex;
    int[] chrEndIndex;
    TagsByTaxaByteHDF5TagGroups tbt;
    int[] tbtRedirect;
    double pThresh = 0.000001;
    int minCount = 20;
    
    int minTagAlleleIntersection = 4;
    
    Task[] jobs;
    static int threadNumPerCore = 32; //32 can be faster, but the admin will complain about too many threads
    int threadNum;
    int chunkSize;
    int chunkNum;
    int chunkStartIndex;
    int chunkEndIndex;
    
    /**
     * Constructor to run genetic mapping
     * @param hapMapHDF5 anchor map
     * @param tbtHDF5 TagsByTaxa of HDF5 format
     * @param outfileS output file of genetic mapping
     * @param pThresh P-value threshold, default: 1e-6
     * @param minCount minimum count when tag appear in taxa (In TBT file), default = 20, too low number lacks statistical power
     * @param coreNum default:-1, which means using all cores in a node. When the coreNum is set less than total core number, which means using coreNum cores, each core runs 1 thread
     * @param chunkSize number of tags in a chunk. This determines the time usage in a node
     */
    public TagAgainstAnchorOriginal (String hapMapHDF5, String tbtHDF5, String outfileS, double pThresh, int minCount, int coreNum, int chunkSize) {
        this.pThresh = pThresh;
        this.minCount = minCount;
        this.loadAnchorMap(hapMapHDF5);
        this.loadTBT(tbtHDF5);
        this.redirect();
        this.chunkSize = chunkSize;
        this.chunkNum = this.getChunkNum(chunkSize);
        this.calculateThreadNum(coreNum);
        this.chunkStartIndex = 0;
        this.chunkEndIndex = chunkNum;
        this.MTMapping(outfileS);
    }
    
    /**
     * Constructor to run genetic mapping on specific chunk of TBT. This is for "qsub" in clusters.
     * @param hapMapHDF5
     * @param tbtHDF5
     * @param outfileS
     * @param pThresh
     * @param minCount
     * @param coreNum
     * @param tagNumPerThreadInChunk
     * @param chunkStartIndex start index of chunk
     * @param chunkEndIndex end index of chunk. Note: chunk of end index is not included.
     */
    public TagAgainstAnchorOriginal (String hapMapHDF5, String tbtHDF5, String outfileS, double pThresh, int minCount, int coreNum, int chunkSize, int chunkStartIndex, int chunkEndIndex) {
        this.pThresh = pThresh;
        this.minCount = minCount;
        this.loadAnchorMap(hapMapHDF5);
        this.loadTBT(tbtHDF5);
        this.redirect();
        this.chunkSize = chunkSize;
        this.chunkNum = this.getChunkNum(chunkSize);
        this.calculateThreadNum(coreNum);
        this.chunkStartIndex = chunkStartIndex;
        this.chunkEndIndex = chunkEndIndex;
        this.MTMapping(outfileS);
    }
    
    /**
     * Return the number of chunks at current chunkSize
     * @param chunkSize
     * @return 
     */
    private int getChunkNum (int chunkSize) {
        int tagNum = tbt.getTagCount();
        int chunkNum = 0;
        int left = tagNum % chunkSize;
        if (left == 0) chunkNum = tagNum / chunkSize;
        else chunkNum = tagNum / chunkSize + 1;
        return chunkNum;
    }
    
    /**
     * Pre-calculate number of chunks when qsub genetic mapping 
     * @param tbtHDF5
     * @param coreNum
     * @param tagNumPerThreadInChunk
     * @return 
     */
    public static int getChunkNum (String tbtHDF5, int chunkSize) {
        TagsByTaxaByteHDF5TagGroups tbt = new TagsByTaxaByteHDF5TagGroups (tbtHDF5);
        int tagNum = tbt.getTagCount();
        int chunkNum = 0;
        int left = tagNum % chunkSize;
        if (left == 0) chunkNum = tagNum / chunkSize;
        else chunkNum = tagNum / chunkSize + 1;
        System.out.println("TBT has " + tagNum + " tags");
        System.out.println("TBT will be devided into " + chunkNum + " chunks, " + chunkSize + " tags each");
        if (left != 0) System.out.println("The last chunk has " + left + " tags");
        System.out.println("The index of chunk are used submit parallel computation to different node");
        System.out.println("Tags in each chunk will be multi-threaded in one node");
        return chunkNum;
    }
    
    /**
     * MT genetic mapping
     * @param outfileS 
     */
    public void MTMapping (String outfileS) {
        if (chunkStartIndex > (chunkNum-1) || chunkEndIndex > chunkNum || chunkStartIndex == chunkEndIndex) {
            System.out.println("Error in setting chunk index. Please reset");
            System.exit(0);
        }
        int tagNum = tbt.getTagCount();
        System.out.println("TBT has " + chunkNum + " chunks, mapping will start at chunk " + chunkStartIndex + ", end at chunk " + chunkEndIndex +"\n");
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
            bw.write("TestTag	TestTagNum	BlastChr	BlastPos	refDiv	LDChr	LDSite	LDPos	BinomP	SigTests	TagTaxaCnt	ChrSig	LRatioB:2	LRatioB:M	SiteOnBestChrThanNextChr	MinSigPos	MaxSigPos");
            bw.newLine();
            for (int i = this.chunkStartIndex; i < this.chunkEndIndex; i++) {
                int chunkStartTagIndex = i * chunkSize;
                int chunkEndTagIndex = chunkStartTagIndex + chunkSize;
                if (chunkEndTagIndex > tagNum) chunkEndTagIndex = tagNum;
                System.out.println("Start mapping tag chunk " + i + "(Index), tag index from " + chunkStartTagIndex + " to " + chunkEndTagIndex);
                int[] threadStartTagIndex;
                int[] threadEndTagIndex;
                int[] threadSize;
                int actualChunkSize = chunkEndTagIndex-chunkStartTagIndex;
                int left = actualChunkSize % threadNum;
                int baseSize = actualChunkSize / threadNum;
                if (baseSize == 0) {
                    threadStartTagIndex = new int[left];
                    threadEndTagIndex = new int[left];
                    for (int j = 0; j < threadStartTagIndex.length; j++) {
                        threadStartTagIndex[j] = chunkStartTagIndex + j;
                        threadEndTagIndex[j] = threadStartTagIndex[j]+1;
                    }
                }
                else {
                    threadSize = new int[threadNum];
                    threadStartTagIndex = new int[threadNum];
                    threadEndTagIndex = new int[threadNum];
                    for (int j = 0; j < left; j++) {
                        threadSize[j] = baseSize + 1;
                    }
                    for (int j = left; j < threadNum; j++) {
                        threadSize[j] = baseSize;
                    }              
                    threadStartTagIndex[0] = chunkStartTagIndex;
                    threadEndTagIndex[0] = threadStartTagIndex[0] + threadSize[0];
                    for (int j = 1; j < threadNum; j++) {
                        threadStartTagIndex[j] = threadEndTagIndex[j-1];
                        threadEndTagIndex[j] = threadStartTagIndex[j]+threadSize[j];
                    }
                }
                int actualThreadNum = threadStartTagIndex.length;                
                jobs = new Task[actualThreadNum];
                Thread[] mts = new Thread[actualThreadNum];
                long lastTimePoint = this.getCurrentTimeNano();
                for (int j = 0; j < actualThreadNum; j++) {
                    jobs[j] = new Task (threadStartTagIndex[j], threadEndTagIndex[j]);
                } 
                System.out.println("Loading this chunk to multiple threads took " + this.getTimeSpanSecond(lastTimePoint) + " seconds");
                System.out.println("Multiple threading mapping in progress...");
                lastTimePoint = this.getCurrentTimeNano();       
                for (int j = 0; j < actualThreadNum; j++) {
                    mts[j] = new Thread(jobs[j]);
                    mts[j].start();
                }
                for (int j = 0; j < actualThreadNum; j++) {
                    try {
                        mts[j].join();
                    }
                    catch (Exception e) {
                        System.out.println(e.toString());
                    }
                }
                System.out.println("Each LD compirison took " + (double)this.getTimeSpanNano(lastTimePoint)/actualChunkSize/anchor.getSiteCount() + " nano seconds");
                System.out.println("Multiple threading mapping took " + this.getTimeSpanSecond(lastTimePoint) + " seconds");
                for (int j = 0; j < jobs.length; j++) {
                    String[] result = jobs[j].getResult();
                    if (result != null)
                    for (int k = 0; k < result.length; k++) {
                        bw.write(result[k]);
                    }
                }
                bw.flush();
                System.out.println("Mapping result from chunk " + i + " was written\n");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Calculate the thread number used in this node
     * @param coreNum 
     */
    private void calculateThreadNum (int coreNum) {
        int numOfProcessors = Runtime.getRuntime().availableProcessors();
        if (coreNum < 0) {
            threadNum = numOfProcessors * threadNumPerCore;
        }
        else {
            if (coreNum == 0) {
                System.out.println("Core number = 0, This runs at least on 1 thread. Quit.");
                System.exit(0);
            }
            threadNumPerCore = 1;
            threadNum = coreNum * threadNumPerCore;
            System.out.println("TBT will be mapped by " + threadNum + " tasks");
            System.out.println("Each core runs 1 tasks, or 1 threads");
        }
        int acutualUseCoreNum = numOfProcessors;
        if (acutualUseCoreNum > threadNum) acutualUseCoreNum = threadNum;
        System.out.println("This node has " + numOfProcessors + " processors. Will use " + acutualUseCoreNum + " processors");
        System.out.println("TBT will be mapped by " + threadNum + " threads");
        System.out.println("Each core runs " + threadNumPerCore + " threads");
        System.out.println("Each TBT chunk has " + chunkSize + " tags, which will be split and mapped by the " + threadNum + " threads");
        System.out.println("Each thread will map " + chunkSize/threadNum + " tags");
        int left = tbt.getTagCount()%chunkSize;
        if (left != 0) {
            System.out.println("The last TBT chunk has " + left + " tags");
        }
        System.out.println("");
	}
    
    /**
     * Class for conducting mapping on one thread
     */
    class Task implements Runnable {
        int tagStartIndex;
        int tagEndIndex;
        TagsByTaxaByte subTBT;
		String[] result = null;

		Task (int tagStartIndex, int tagEndIndex) {
			this.buildSubTBT(tagStartIndex, tagEndIndex);
		}
        
        /**
         * Build TBTByte for each task/thread
         * @param tagStartIndex
         * @param tagEndIndex 
         */
        private void buildSubTBT (int tagStartIndex, int tagEndIndex) {
            this.tagStartIndex = tagStartIndex;
            this.tagEndIndex = tagEndIndex;
            int tagCount = tagEndIndex - tagStartIndex;
            long[][] tags = new long[tbt.getTagSizeInLong()][tagCount];
            byte[] tagLength = new byte[tagCount];
            byte[][] tagDist = new byte[tbt.getTaxaCount()][tagCount];
            String[] namesForTaxa = tbt.getTaxaNames();
            for (int i = 0; i < tagCount; i++) {
                int tagIndex = i + tagStartIndex;
                long[] t = tbt.getTag(tagIndex);
                for (int j = 0; j < tags.length; j++) {
                    tags[j][i] = t[j];
                }
                tagLength[i] = (byte)tbt.getTagLength(tagIndex);
                for (int j = 0; j < tbt.getTaxaCount(); j++) {
                    tagDist[j][i] = (byte)tbt.getReadCountForTagTaxon(tagIndex, j);
                }
                
            }
            subTBT = new TagsByTaxaByte (tags, tagLength, tagDist, namesForTaxa);
        }
        
        /**
         * Return result for output when all the tasks are done
         * @return 
         */
		public String[] getResult () {
			return result;
		}

        @Override
		public void run() {
            long lastTimePoint = getCurrentTimeNano();
			ArrayList<String> resultList = new ArrayList();          
            for (int i = 0; i < subTBT.getTagCount(); i++) {
                if(subTBT.getNumberOfTaxaWithTag(i) <minCount) continue;
                double[] bestR={-1,-1,-1, 1, -1};
                int blastChr=Integer.MIN_VALUE, blastPos=Integer.MIN_VALUE, bestAlignment=Integer.MIN_VALUE, refDiv=Integer.MIN_VALUE;
                long[] testTag = subTBT.getTag(i);
                long[] testTagDist;
                double[][] theResults=new double[chromosomeNumber.length][];
                testTagDist=getTagsInBits(subTBT, i, tbtRedirect, anchor.getSequenceCount());
                ScanChromosome[] scanOnChr = new ScanChromosome[chromosomeNumber.length];
                for (int j = 0; j < chromosomeNumber.length; j++) {
                    scanOnChr[j] = new ScanChromosome(testTagDist, theResults, j, chrStartIndex[j], chrEndIndex[j], pThresh, 1, blastPos);
                    scanOnChr[j].scan();
                }
                int countRealSig=0;
                double[] pRank=new double[chromosomeNumber.length];
                for (int j = 0; j < chromosomeNumber.length; j++) {
                    double[] r=theResults[j];
                    pRank[j]=r[3];
                    if(r[3]<bestR[3]) {bestR=r.clone(); bestAlignment=j;}
                    if(r[3]<pThresh) countRealSig++;
                }
                if (bestR[3] == 1) continue;
                Arrays.sort(pRank);
                double[][] bestResWithNewThreshold=new double[1][];
                ScanChromosome bestChrNewThres = new ScanChromosome(testTagDist, bestResWithNewThreshold, 0, chrStartIndex[bestAlignment], chrEndIndex[bestAlignment],pRank[1], 1, blastPos);
                bestChrNewThres.scan();
                int countOfSitesBetterThanNextBestChr=(int)bestResWithNewThreshold[0][4];
                String s=String.format("%s %d %d %d %d %d %d %d %g %d %d %d %g %g %d %d %d %n",BaseEncoder.getSequenceFromLong(testTag), subTBT.getReadCount(i), blastChr, blastPos, refDiv,
                        (int)bestR[0],(int)bestR[1],(int)bestR[2], bestR[3], (int)bestR[4], subTBT.getNumberOfTaxaWithTag(i),
                        countRealSig, Math.log10(pRank[1]/pRank[0]), Math.log10(pRank[chromosomeNumber.length/2]/pRank[0]), countOfSitesBetterThanNextBestChr,
                        bestChrNewThres.minSigPos, bestChrNewThres.maxSigPos);
                resultList.add(s);
            }
			
			result = resultList.toArray(new String[resultList.size()]);
			//System.out.println("Task with " + subTBT.getTagCount() + " tags indexing from " + tagStartIndex + " to " + tagEndIndex + " are finished in " + getTimeSpanSecond(lastTimePoint) + " seconds");
		}

	}
    
    /**
     * Class used to scan SNPs in one chromosome
     */
    private class ScanChromosome {
        int chrIndex;
        int chrStartIndex;
        int chrEndIndex;
        OpenBitSet obsTdist;
        Binomial binomFunc=new Binomial(5, 0.5, new RandomJava());
        double[][] resultReport;
        double sigThreshold;
        int minSigPos=Integer.MIN_VALUE, maxSigPos=Integer.MIN_VALUE;
        int step=1;
        double bionomialThreshold=0.2;
        int blockPosition=Integer.MIN_VALUE;  //position to block from testing.
        int blockWindow=100;
        
        /**
         * 
         * @param tdist tag distribution if bits across all taxa
         * @param resultReport
         * @param chrIndex 
         * @param chrStartIndex SNP start index of this chromosome
         * @param chrEndIndex SNP end index of this chromosome
         * @param sigThreshold threadhold of P-value, default 1E-6
         * @param step step of scanning
         * @param blockPosition 
         */
        public ScanChromosome(long[] tdist, double[][] resultReport, int chrIndex, int chrStartIndex, int chrEndIndex, double sigThreshold, int step, int blockPosition) {
            obsTdist=new OpenBitSet(tdist,tdist.length);
            this.resultReport=resultReport;
            this.chrIndex = chrIndex;
            this.chrStartIndex = chrStartIndex;
            this.chrEndIndex = chrEndIndex;
            this.sigThreshold=sigThreshold;
            this.blockPosition=blockPosition;
            this.step=step;
        }

        public void scan() {
            long tests=0;
            int bestSite=-1, countSig=0;
            double bestP=2;
            for (int i = chrStartIndex; i < chrEndIndex; i+=step) {
                if(Math.abs(anchor.getPositionInLocus(i)-blockPosition)<blockWindow) continue;
                OpenBitSet obsMajor = new OpenBitSet(anchor.getAllelePresenceForAllTaxa(i, 0).getBits());
                OpenBitSet obsMinor = new OpenBitSet(anchor.getAllelePresenceForAllTaxa(i, 1).getBits());
                if(obsMinor.cardinality()>4) {
                    double p=fastTestSites(obsTdist, obsMajor, obsMinor, anchorMaf[i], binomFunc);
                    if(p<bestP) {bestP=p; bestSite=i;}
                    if(p<sigThreshold) {
                        countSig++;
                        if(minSigPos==-1) minSigPos=anchor.getPositionInLocus(i);
                        maxSigPos=anchor.getPositionInLocus(i);
                    }
                }
                tests++;
            }
            int chr=Integer.parseInt(anchor.getLocus(bestSite).getChromosomeName());
            double[] result={chr, bestSite, anchor.getPositionInLocus(bestSite),bestP, countSig};
            resultReport[chrIndex]=result;
        }
    }
    
    /**
     * fast test, using ratio to reduce the calculation of binomial and p (reduce 73%), but the time is only saved by 3%
     * @param obsTdist
     * @param obsMajor
     * @param obsMinor
     * @param maf
     * @param binomFunc
     * @return 
     */
    public double fastTestSites(OpenBitSet obsTdist, OpenBitSet obsMajor, OpenBitSet obsMinor, double maf, Binomial binomFunc) {
        double result=1;
        int tagMinorCount=0, tagMajorCount=0;
        tagMinorCount=(int)OpenBitSet.intersectionCount(obsTdist, obsMinor);
        tagMajorCount=(int)OpenBitSet.intersectionCount(obsTdist,obsMajor);

        int sumTagAllele = tagMinorCount + tagMajorCount;
        if(sumTagAllele<4) return result;
        
        double ratio;
        double minorProb;
        int cdfCount;
        if (tagMinorCount < tagMajorCount) {
            ratio = (double)tagMinorCount/sumTagAllele;
            minorProb = maf;
            cdfCount = tagMinorCount;
        }
        else {
            ratio = (double)tagMajorCount/sumTagAllele;
            minorProb = 1-maf;
            cdfCount = tagMajorCount;
        }
        if ((ratio-minorProb) < -0.003) return result;
        binomFunc.setNandP(sumTagAllele,minorProb);
        try {
            result = binomFunc.cdf(cdfCount);
            //System.out.println(ratio+"\t"+minorProb+"\t"+result);
  
        } catch (Exception e) {
            System.err.println("Error in the BinomialDistributionImpl");
        }
        return result;
    }
    
    /**
     * Test association using binomial test, collecting P-value
     * @param obsTdist
     * @param obsMajor
     * @param obsMinor
     * @param binomFunc
     * @return 
     */
    public double testSites(OpenBitSet obsTdist, OpenBitSet obsMajor, OpenBitSet obsMinor, double maf, Binomial binomFunc) {
        double result=1;
        int tagMinorCount=0, tagMajorCount=0;
        tagMinorCount=(int)OpenBitSet.intersectionCount(obsTdist, obsMinor);
        tagMajorCount=(int)OpenBitSet.intersectionCount(obsTdist,obsMajor);

        int sumTagAllele = tagMinorCount + tagMajorCount;
        if(sumTagAllele<4) return result;
    
        double minorProb = (tagMinorCount<tagMajorCount)? maf:1-maf;
        binomFunc.setNandP(sumTagAllele,minorProb);
        try {
            result = (tagMinorCount<tagMajorCount)?binomFunc.cdf(tagMinorCount):binomFunc.cdf(tagMajorCount);
  
        } catch (Exception e) {
            System.err.println("Error in the BinomialDistributionImpl");
        }
        return result;
    }
    
    /**
     * Test association using binomial test, collecting P-value
     * @param obsTdist
     * @param obsMajor
     * @param obsMinor
     * @param binomFunc
     * @return 
     */
    public static double testSites(OpenBitSet obsTdist, OpenBitSet obsMajor, OpenBitSet obsMinor, Binomial binomFunc) {
        double result=1;
        int minorCount=0, tagMinorCount=0, majorCount=0, tagMajorCount=0;
        tagMinorCount=(int)OpenBitSet.intersectionCount(obsTdist, obsMinor);
        tagMajorCount=(int)OpenBitSet.intersectionCount(obsTdist,obsMajor);
        minorCount=(int)obsMinor.cardinality();
        majorCount=(int)obsMajor.cardinality();

        int sumAllele = minorCount + majorCount;
        int sumTagAllele = tagMinorCount + tagMajorCount;
        if(sumTagAllele<4) return result;
    
        double minorProb = (tagMinorCount<tagMajorCount)?(double)minorCount/(double)sumAllele:(double)majorCount/(double)sumAllele;
        binomFunc.setNandP(sumTagAllele,minorProb);
        try {
            result = (tagMinorCount<tagMajorCount)?binomFunc.cdf(tagMinorCount):binomFunc.cdf(tagMajorCount);
  
        } catch (Exception e) {
            System.err.println("Error in the BinomialDistributionImpl");
        }
        return result;
    }
    
    /**
     * Return bit version of tag distribution across taxa
     * @param aTBT
     * @param tagIndex
     * @param reDirect
     * @param anchorTaxa
     * @return 
     */
    private long[] getTagsInBits(TagsByTaxaByte aTBT, int tagIndex, int[] reDirect, int anchorTaxa) {
        int lgPerSite = (anchorTaxa / 64) + 1;
        long[] seq = new long[lgPerSite];
        for (int j = 0; j < aTBT.getTaxaCount(); j++) {
            if(reDirect[j]<0) continue;
            int index=reDirect[j]/64;
            int offset=reDirect[j]%64;
            if (aTBT.getReadCountForTagTaxon(tagIndex, j)>0) {  //reference alleles
                seq[index]=seq[index]|(1L<<offset);
            }
        }
        return seq;
	}
    
    private void redirect () {
        long lastTimePoint = this.getCurrentTimeNano();
        tbtRedirect = new int[tbt.getTaxaCount()];
        IdGroup g = anchor.getIdGroup();
        for (int i = 0; i < tbtRedirect.length; i++) {
            tbtRedirect[i] = g.whichIdNumber(tbt.getTaxaName(i));
        }
        System.out.println("Taxa redirection took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
    }
    
    /**
     * load HDF5 TBT
     * @param tbtHDF5 
     */
    public void loadTBT (String tbtHDF5) {
        long lastTimePoint = this.getCurrentTimeNano();
        tbt = new TagsByTaxaByteHDF5TagGroups (tbtHDF5);
        System.out.println("Loading TBT HDF5 took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
        System.out.println("TBT has " + tbt.getTagCount() + " tags and " + tbt.getTaxaCount() + " taxa");
    }
    
    /**
     * load HDF5 genotype and convert to BitAlignment
     * @param hapMapHDF5 
     */
    private void loadAnchorMap (String hapMapHDF5) {
        long lastTimePoint = this.getCurrentTimeNano();
        Alignment a = ImportUtils.readGuessFormat(hapMapHDF5, true);
        System.out.println("Loading hapmap HDF5 took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
        System.out.println("The anchor map has " + a.getSiteCount() + " sites and " + a.getTaxaCount() + " taxa");
        lastTimePoint = this.getCurrentTimeNano();
        int[] chrOffSet = a.getLociOffsets();
        chromosomeNumber = new int[a.getLoci().length];
        chrStartIndex = new int[chromosomeNumber.length];
        chrEndIndex = new int[chromosomeNumber.length];
        for (int i = 0; i < chromosomeNumber.length; i++) {
            chromosomeNumber[i] = a.getLoci()[i].getChromosomeNumber();
            chrStartIndex[i] = chrOffSet[i];
            chrEndIndex[i] = chrOffSet[i] + a.getLocusSiteCount(a.getLoci()[i]);
        }
        anchor = (BitAlignment)BitAlignment.getHomozygousNucleotideInstance(a, true);
        anchorMaf = new double[anchor.getSiteCount()];
        for (int i = 0; i < anchorMaf.length; i++) {
            anchorMaf[i] = anchor.getMinorAlleleFrequency(i);
        }
        System.out.println("Loading and converting to BitAlignment took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
        this.screenPrintGbMemoryCurrentUse();
        this.screenPrintGbMemoryAvailable();
        System.out.println();
    }
    
    private double getTimeSpanSecond (long lastTimePoint) {
        return (double)this.getTimeSpanNano(lastTimePoint)/1000000000;
    }
    
    private long getTimeSpanNano (long lastTimePoint) {
        return this.getCurrentTimeNano()- lastTimePoint;
    }
    
    private long getCurrentTimeNano () {
        return System.nanoTime();
    }
    
    private String getCurrentTimeHR() {
		SimpleDateFormat sd = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        Date date = new Date();
        String str = sd.format(date);
        return str;
	}
    
    private double getGbMemoryAvailable () {
        return (double)(Runtime.getRuntime().maxMemory()- (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()))/1024/1024/1024;
    }
    
    private double getGbMemoryCurrentUse () {
        return (double)(Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/1024/1024/1024;
    }
    
    private void screenPrintTimeSpanSecond (long lastTimePoint) {
        System.out.println("Time span is " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
    }
    
    private void screenPrintTimeSpanNano (long lastTimePoint) {
        System.out.println("Time span is " + String.valueOf(this.getTimeSpanNano(lastTimePoint)) + " ns");
    }
    
    private void screenPrintCurrentTimeHR () {
        System.out.println("Current time is " + this.getCurrentTimeHR());
    }
    
    private void screenPrintGbMemoryAvailable () {
        System.out.println("Available memory is " + String.valueOf(this.getGbMemoryAvailable()) + " GB");
    }
    
    private void screenPrintGbMemoryCurrentUse () {
        System.out.println("Current memory in use is " + String.valueOf(this.getGbMemoryCurrentUse()) + " GB");
    }
}
