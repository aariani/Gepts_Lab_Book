/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import ch.systemsx.cisd.hdf5.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import net.maizegenetics.gbs.maps.PETagsOnPhysicalMap;
import net.maizegenetics.gbs.maps.PETagsOnPhysicalMapV3;
import net.maizegenetics.gbs.maps.TagMappingInfoV3;
import net.maizegenetics.gbs.maps.TagsOnGeneticMap;
import net.maizegenetics.gbs.maps.TagsOnPhysMapHDF5;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMapV3;
import net.maizegenetics.gbs.tagdist.PETagCounts;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.util.OpenBitSet;
import org.apache.log4j.Logger;

/**
 *
 * @author Fei Lu
 */
public class FeiPipelines {
    private Logger myLogger = Logger.getLogger(FeiPipelines.class);
    
    public FeiPipelines () {
        //this.pipelinePE();
        //this.testPipeline();
        this.GBSV3Pipeline();
    }
    
    public static void main (String[] args) {
        new FeiPipelines();
        
    }
    
    public void GBSV3Pipeline () {
        //this.addTripsicum();
        //this.mkV3Alignment(); //take some time to align using blast+
        //this.initializeV3TOPM();
        //this.annotateV3TOPM();
        this.geneticMapping();
    
    }
    
    public void geneticMapping () {
        //String hapMapHDF5 = "/workdir/mingh/AllZeaGBSv27i3b.imp.hmp.h5";
        //String hapMapHDF5 = "/workdir/mingh/GBS27.small.imp.hmp.h5";
        String hapMapHDF5 = "M:/GBSV3/genotype/GBS27.small.imp.hmp.h5";
        //String tbtHDF5 = "/workdir/mingh/smallerTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String tbtHDF5 = "M:/GBSV3/tbt/smallerTBTHDF5_mergedtaxa_pivot_20120921.h5";
        //String outfileS = "/workdir/mingh/outfile.txt";
        String outfileS = "M:/GBSV3/mappingResult/outfile.txt";
        TagAgainstAnchor taa = new TagAgainstAnchor(hapMapHDF5, tbtHDF5, outfileS, 0.000001, 20, -1, 4096);
        //TagAgainstAnchor.getChunkNum(tbtHDF5, 1024);    
    }
    
    public void annotateV3TOPM () {
        String inputFileS = "M:/GBSV3/topm/v3.topm.h5";
        //String inputFileS = "/workdir/mingh/v3.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(inputFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        String bowtie2SamFileS = "M:/GBSV3/alignment/v3.bowtie2.sam.gz";
        //anno.annotateWithBowtie2(bowtie2SamFileS, 5);
        String bwaSamFileS = "M:/GBSV3/alignment/v3.bwa.sam.gz";
        //anno.annotateWithBWA(bwaSamFileS, 5);
        //String blastDirS = "M:/GBSV3/alignment/blastResult/";
        String blastDirS = "/workdir/mingh/blastResult/";
        //anno.annotateWithBlastFromDir(blastDirS, 5);
        String PETOPMFileS = "M:/af/ptopm/PE.topm";
        //anno.annotateWithPE(PETOPMFileS, 5);
        String TOGMFileS = "M:/pav/PhyGenMapping/v1.togm.txt";
        anno.annotateWithGM(TOGMFileS, 1);
    }
    
    public void initializeV3TOPM () {
        String inputFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        String outputFileS = "M:/GBSV3/topm/v3.topm.h5";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Bit);
        TagsOnPhysicalMapV3.createFile(tc, outputFileS);
    }
    
    /**
     * make alignment using bowtie2m, bwa and blast+, fasta file need to  be split to many small files to speed up blast (parallel)
     */
    public void mkV3Alignment () {
        String inputFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        String fastqFileS = "M:/GBSV3/alignment/v3.fastq";
        //new FeiUtils ().convertTagCount2Fastq(inputFileS, fastqFileS);
        
        String fastaFileS = "M:/GBSV3/alignment/v3.fasta";
        //TagCounts tc = new TagCounts(inputFileS, FilePacking.Bit);
        //tc.toFASTA(fastaFileS);
        
        String zipFastaFileS = "M:/GBSV3/alignment/v3.fasta.gz";
        String fastaDirS = "M:/GBSV3/alignment/splitFasta/";
        FeiUtils util = new FeiUtils();
        util.splitFastaFileS(zipFastaFileS, fastaDirS, 20000);  
    }
    
    public void addTripsicum () {
        String sourceDir = "M:/GBSV3/tagCount/source/";
        String mergeFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        String arguments = "-i " + sourceDir + " -o " + mergeFileS + " -c " + 10 ;
        String[] args = arguments.split(" ");
        MergeMultipleTagCountPlugin m = new MergeMultipleTagCountPlugin();
        m.setParameters(args);
		m.performFunction(null);  
    }
    
    public void testPipeline () {
        //this.mkSmallTagCount();
        //this.mkFastq();
        //this.mkFasta();
        //this.mkTOPM(); //old version
        //this.mkTOPMHDF5(); //old version
        //this.initializeTOPMHDF5();
        //this.annotateTOPMHDF5WithAligner();
        this.readTOPMHDF5();
        //this.mkSmallAnchorHDF5();
        //this.mkSmallTBTHDF5();
    }
    
    public void mkSmallTBTHDF5 () {
        String inputTBTS = "M:/GBSV3/tbt/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String outputTBTS = "M:/GBSV3/tbt/smallTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String ooutputTBTS = "M:/GBSV3/tbt/smallerTBTHDF5_mergedtaxa_pivot_20120921.h5";
        TagsByTaxaByteHDF5TagGroups tbt = new TagsByTaxaByteHDF5TagGroups (outputTBTS);
        int[] tagIndex = new int[4096];
        for (int i = 0; i < tagIndex.length; i++) {
            tagIndex[i] = i;
        }
        tbt.writeDistFile(ooutputTBTS, tagIndex);
        
    }
    
    public void mkSmallAnchorHDF5 () {
        //String hapMapInputS = "M:/GBSV3/genotype/AllZeaGBSv27i3b.imp.hmp.h5";
        //String hapMapOutputS = "M:/GBSV3/genotype/GBS27.small.imp.hmp.h5";
        String hapMapInputS = "./AllZeaGBSv27i3b.imp.hmp.h5";
        String hapMapOutputS = "./GBS27.small.imp.hmp.h5";
        Alignment a = ImportUtils.readGuessFormat(hapMapInputS, true);
        
        int[] snpIndex = new int[1024];
        int startIndex = 0;
        for (int i = 0; i < snpIndex.length; i++) {
            snpIndex[i] = i+startIndex;
        }

        ExportUtils.writeToMutableHDF5(a, hapMapOutputS, snpIndex);
    }
    
    public void readTOPMHDF5 () {
        //String inputFileS = "M:/GBStest/topm/ini.topm.h5";
        String inputFileS = "M:/GBSV3/topm/v3.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3 (inputFileS);
        System.out.println(topm.getTagCount());
 /*       
        for (int i = 0; i < 500000; i++) {
            TagMappingInfoV3[] tmi = new TagMappingInfoV3[topm.getMappingNum()];
            for (int j = 0; j < topm.getMappingNum(); j++) {
                tmi[j] = topm.getMappingInfo(i, j);
            }
            
            int last = topm.getMappingNum()-1;
            if (tmi[last].chromosome < 0) continue;
            if (tmi[0].chromosome < 0) continue;
            if (tmi[5].chromosome < 0) continue;
            if (tmi[10].chromosome < 0) continue;
            int dBowtie = Math.abs((tmi[0].chromosome*1000000000+tmi[0].startPosition) - (tmi[last].chromosome*1000000000+tmi[last].startPosition));
            int dBwa = Math.abs((tmi[5].chromosome*1000000000+tmi[5].startPosition) - (tmi[last].chromosome*1000000000+tmi[last].startPosition));
            int dBlast = Math.abs((tmi[10].chromosome*1000000000+tmi[10].startPosition) - (tmi[last].chromosome*1000000000+tmi[last].startPosition));
            System.out.println(dBowtie+"\t"+dBwa+"\t"+dBlast+"\t"+i);
        }
 */
    }
    
    public void annotateTOPMHDF5WithAligner () {
        String inputFileS = "M:/GBStest/topm/ini.topm.h5";

        
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(inputFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        String bowtie2SamFileS = "M:/GBStest/alignment/bowtie2-K5.sam";
        anno.annotateWithBowtie2(bowtie2SamFileS, 5);
        String bwaSamFileS = "M:/GBStest/alignment/bwa-N5.sam";
        anno.annotateWithBWA(bwaSamFileS, 5);
        //String blastFileS = "M:/GBStest/alignment/blast.m8.txt";
        //anno.annotateWithBLAST(blastFileS, 5);
        String blastDirS = "M:/GBStest/alignment/blastOut/";
        anno.annotateWithBlastFromDir(blastDirS, 5);
        String PETOPMFileS = "M:/af/ptopm/PE.topm";
        anno.annotateWithPE(PETOPMFileS, 5);
        String TOGMFileS = "M:/pav/PhyGenMapping/v1.togm.txt";
        anno.annotateWithGM(TOGMFileS, 1);
    }
    
    public void initializeTOPMHDF5 () {
        String inputFileS = "M:/GBStest/tagCount/small.cnt";
        String outputFileS = "M:/GBStest/topm/ini.topm.h5";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Bit);
        TagsOnPhysicalMapV3.createFile(tc, outputFileS);
    }
    
    public void mkTOPMHDF5 () {
        String infileS = "M:/GBStest/topm/bowtie2.topm.bin";
        String outfileS = "M:/GBStest/topm/bowtie2.topm.h5";
        TagsOnPhysicalMap topm = new TagsOnPhysicalMap (infileS, false);
        TagsOnPhysMapHDF5.createFile(topm, outfileS, 4, 16);
    }
    
    public void mkTOPM () {
        String inputFileS = "M:/GBStest/alignment/bowtie2.sam";
        String outputFileS = "M:/GBStest/topm/bowtie2.topm.bin";
        new FeiUtils ().convertSAM2TOPM(inputFileS, outputFileS);
    }
    
    /**
     * Make FASTA file and do alignment using Blast, -m 8 -e 1e-10
     */
    public void mkFasta () {
        String inputFileS = "M:/GBStest/tagCount/small.cnt";
        String outputFileS = "M:/GBStest/alignment/small.fa";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Bit);
        tc.toFASTA(outputFileS);
    }
    
    /**
     * Make Fastq file and do alignment using Bowtie2 and BWA
     * Alignment files should be in the alignment folder
     */
    public void mkFastq () {
        String inputFileS = "M:/GBStest/tagCount/small.cnt";
        String outputFileS = "M:/GBStest/alignment/small.fq";
        new FeiUtils ().convertTagCount2Fastq(inputFileS, outputFileS);
    }
    
    public void mkSmallTagCount () {
        String inputFileS = "M:/GBStest/tagCount/434GFAAXX_s_4.cnt";
        String outputFileS = "M:/GBStest/tagCount/small.cnt";
        new FeiUtils ().mkSmallTagCountsFile(inputFileS, outputFileS, 10001, 500);
    }
    
    public void pipelinePE () {
        //parseQseq();
        //this.parseFastq();
        //this.checkPETagCounts(); //for checking, not in the pipeline.
        //this.mergePETagCounts();
        //this.contigPETagCounts();
        //this.mkPEstatistics();//for presentation, not included in pipeline.
        //this.mkFastaFileFromPE();
        this.mkPEAlignment();
        
        //************************************
        //Deprecated
        //this.alignmentStep1();
        //this.alignmentStep2();
        //this.checkAlignmentOfPE(); //check alignment of longer sequence
        //************************************
    }
    
    public void checkAlignmentOfPE () {
        String TOGMFileS = "M:/pav/PhyGenMapping/v1.togm.txt";
        String topmFileS = "N:/Zea/AllZeaBuild_2.X/04_TOPM/2.6_production/02_MergedTOPM/AllZeaGBS_v2.6_MergedUnfiltProdTOPM_20130425.topm";
        String ptopmFileS = "M:/af/ptopm/merge.ptopm";
        String compareTableS = "E:/Research/af/alignmentImprovement/alignmentCompare.txt";
        FeiUtils fu = new FeiUtils ();
        fu.mkAlignmentCompareTable(TOGMFileS, topmFileS, ptopmFileS, compareTableS);
    }
    
    public void alignmentStep2 () {
        String fSamFileS = "M:/af/alignment/f.sam";
        String bSamFileS = "M:/af/alignment/b.sam";
        String contigSamFileS = "M:/af/alignment/c.sam";
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        PETagsOnPhysicalMap topm = new PETagsOnPhysicalMap (ptc, fSamFileS, bSamFileS, contigSamFileS);
        String ptopmFileS = "M:/af/ptopm/merge.ptopm";
        topm.writeDistFile(ptopmFileS, FilePacking.Text);
    }
    
    public void alignmentStep1 () {
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        String fastaFFileS = "M:/af/alignment/f.fasta.txt";
        String fastaBFileS = "M:/af/alignment/b.fasta.txt";
        String fastaCFileS = "M:/af/alignment/c.fasta.txt";
        ptc.mkFastaFile(fastaFFileS, fastaBFileS, fastaCFileS);
    }
    
    public void mkPEAlignment () {
        String fFastaFileS = "M:/af/alignment/f.fasta.txt";
        String bFastaFileS = "M:/af/alignment/b.fasta.txt";
        String fSamFileS = "M:/af/alignment/f.k5.sam.gz";
        String bSamFileS = "M:/af/alignment/b.k5.sam.gz";
        String a = null;
        PETagsOnPhysicalMapV3 ptopm = new PETagsOnPhysicalMapV3 (fFastaFileS, bFastaFileS, fSamFileS, bSamFileS);
        String PETOPM = "M:/af/ptopm/PE.topm";
        ptopm.writeBinaryFile(PETOPM);
        ptopm = new PETagsOnPhysicalMapV3(PETOPM);
    }
    
    public void mkFastaFileFromPE () {
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        String fFastaFileS = "M:/af/alignment/f.fasta.txt";
        String bFastaFileS = "M:/af/alignment/b.fasta.txt";
        ptc.mkFastaFile(fFastaFileS, bFastaFileS);
    }
    
    public void mkPEstatistics () {
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        System.out.println("Tag count is " + ptc.getTagCount());
        System.out.println("Read count is " + ptc.getTotalReadCount());
        System.out.println("Contig count is " + ptc.getContigCount());
        String staFileS = "E:/Research/af/PEstatistics/sta.txt";
        int maxCount = 1;
        for (int i = 0; i < ptc.getTagCount(); i++) {
            if (ptc.getReadCount(i) > maxCount) maxCount = ptc.getReadCount(i);
            //System.out.println(maxCount);            
        }
        System.out.println("Max read count is " + maxCount);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(staFileS), 65536);
            bw.write("FLength\tBLength\tContigLength\tReadCount");
            bw.newLine();
            for (int i = 0; i < 10000; i++) {
                //if (ptc.getContigLengthInLong(i) == 0) continue; 
                bw.write(String.valueOf(ptc.getTagFLength(i))+"\t"+String.valueOf(ptc.getTagBLength(i))+"\t"+String.valueOf(ptc.getContigLength(i))+"\t"+String.valueOf(ptc.getReadCount(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
     public void contigPETagCounts () {
        String infileS = "M:/af/mergePETagCounts/merge.pe.cnt";
        String outfileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        String arguments = "-i " + infileS + " -o " + outfileS;
        String[] args = arguments.split(" ");
        ContigPETagCountPlugin m = new ContigPETagCountPlugin();
        m.setParameters(args);
		m.performFunction(null);
    }
    
    public void mergePETagCounts () {
        String inputDirS = "M:/af/PETagCounts/";
        String outfileS = "M:/af/mergePETagCounts/merge.pe.cnt";
        String arguments = "-i " + inputDirS + " -o " + outfileS;
        String[] args = arguments.split(" ");
        MergePETagCountPlugin m = new MergePETagCountPlugin();
        m.setParameters(args);
		m.performFunction(null); 
    }

    public void parseFastq () {
        String infile1 = "M:/af/Illumina/ImputationP15_1_1_fastq.txt";
        String infile2 = "M:/af/Illumina/ImputationP15_1_2_fastq.txt";
        
        String keyfile = "M:/af/key/ImputationP15_key.txt";
        String outputDirS = "M:/af/PETagCounts/";
        String arguments = "-iF " + infile1 + " -iB " + infile2 + " -k " + keyfile + " -e ApekI -l 8 -o " + outputDirS;
        String[] args = arguments.split(" ");
        FastqToPETagCountPlugin q = new FastqToPETagCountPlugin();
        q.setParameters(args);
		q.performFunction(null);     
    }
    
    public void checkPETagCounts () {
        String outputDirS = "M:/af/PETagCounts/";
        String txtDirS = "M:/af/text/";
        File[] files = new File (outputDirS).listFiles();
        for (int i = 0; i < files.length; i++) {
            PETagCounts p = new PETagCounts (files[i].getAbsolutePath(), FilePacking.Bit);
            String out = txtDirS + "/" + files[i].getName();
            p.writeDistFile(out, FilePacking.Text, 0);
        }
    }
     
    public void parseQseq () {
        String infile1 = "M:/af/Illumina/81546ABXX_8_1_qseq.txt";
        String infile2 = "M:/af/Illumina/81546ABXX_8_2_qseq.txt";
        String keyfile = "M:/af/key/81546ABXX_key.txt";
        String outputDirS = "M:/af/PETagCounts/";
        String arguments = "-iF " + infile1 + " -iB " + infile2 + " -k " + keyfile + " -e ApekI -o " + outputDirS;
        String[] args = arguments.split(" ");
        QseqToPETagCountPlugin q = new QseqToPETagCountPlugin();
        q.setParameters(args);
		q.performFunction(null);     
    }
    
}
