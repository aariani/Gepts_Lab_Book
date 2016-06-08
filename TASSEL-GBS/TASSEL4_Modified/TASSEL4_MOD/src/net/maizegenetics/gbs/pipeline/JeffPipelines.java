/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.File;
import java.util.Arrays;
import net.maizegenetics.baseplugins.ExtractHapmapSubsetPlugin;
import net.maizegenetics.gbs.maps.AbstractTagsOnPhysicalMap;
import net.maizegenetics.gbs.maps.TOPMUtils;
import net.maizegenetics.gbs.maps.TagsOnPhysMapHDF5;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.prefs.TasselPrefs;

/**
 *
 * @author terry
 */
public class JeffPipelines {

    public static void main(String[] args) {
//        helloWorld();
//        runDiscoverySNPCallerPlugin();
//        runExtractHapmapSubsetPlugin();
//        runCompareGenosBetweenHapMapFilesPlugin();
//        expandVariantsInTOPM();
//        getChrsFromTOPM();
//        runTOPMSummaryPlugin();
//        convertTOPMtoHDF5();
        runProductionSNPCallerPlugin();
//        convertHDF5ToHapMap();
    }
    
    public static void helloWorld() {
        System.out.println("Hello world!");
    }

    public static void runDiscoverySNPCallerPlugin() {

        String baseDirMDPLowVol = "/Users/jcg233/Documents/GBS/MDP1_low_vol/";
        String[] MDPLowVolArgs = new String[]{
            "-i", baseDirMDPLowVol + "C08L7ACXX_6_min2.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirMDPLowVol + "tassel4/hapmap/biSNPWGap/MDP1_low_vol_WGap3rdAlleleUpdateTOPM2.chr+.hmp.txt",
//            "-vcf",
            "-m", baseDirMDPLowVol + "MGP1_low_vol_min2_wPosit.topm.bin",
            "-mUpd", baseDirMDPLowVol + "MGP1_low_vol_min2_wPosit_wVariants2.topm.bin",
            "-ref", "/Users/jcg233/Documents/GBS/refGenome/ZmB73_RefGen_v2.fa",  // baseDirMDPLowVol + "maize_agp_v2_chr10.fasta",
//            "-LocusBorder", "150",
            "-p", baseDirMDPLowVol + "MDP1_betterFake_ped.txt", 
            "-mnF", "0.8",
            "-mnMAF", "0.005",
            "-mnMAC", "99999", // this will never be satified: this way -mnMAF overrides it
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclRare",
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
            "-callBiSNPsWGap",  // call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-cF", // customFiltering (inbredCoverage & inbredHetScore)
            "-sC", "10", // Start chromosome
            "-eC", "10" // End chromosome
        };

        String baseDirJulyBuild = "/Users/jcg233/Documents/GBS/";
        String[] JulyBuildArgs = new String[]{
            "-i", baseDirJulyBuild + "20120701Build/tbt/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5",
            "-m", baseDirJulyBuild + "20120701Build/topm/AllZeaMasterTags_c10_20120703.topm",
            "-p", baseDirJulyBuild + "20120701Build/AllZeaPedigree2012oct01C.txt",
            "-o", baseDirJulyBuild + "20120701Build/hapmap",
            "-ref", baseDirJulyBuild + "refGenome/ZmB73_RefGen_v2.fa",
            "-mnF", "0.8",
            "-mnMAF", "0.001",
            "-mnMAC", "10",
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            "-mxSites", "2000",
//            "-y", // use TagsByTaxaByte
//            "-mUpd", baseDir+"",
//            "-LocusBorder", "150",
//            "-inclRare",
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-sC", "7", // Start chromosome
            "-eC", "7" // End chromosome
        };

        String baseDirGrape = "/Users/jcg233/Documents/GBS/grape/";
        String[] GrapeArgs = new String[]{
            "-i", baseDirGrape + "mergedtbt/paola_test_old.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirGrape + "hapmap/vcf/grapeTest_chr+.hmp.txt",
            "-vcf",
            "-m", baseDirGrape + "topm/paola_test_old.topm.bin",
//            "-mUpd", baseDirGrape + "",
//            "-ref", baseDirGrape + "Vitis_vinifera.IGGP_12x.dna.toplevel.clean_chr_name.fa",
//            "-LocusBorder", "150",
//            "-p", baseDirGrape + "MDP1_betterFake_ped.txt", 
//            "-mnF", "0.8",
            "-mnMAF", "0.005",
            "-mnMAC", "99999", // this will never be satified: this way -mnMAF overrides it
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclRare",
            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  // call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-sC", "1", // Start chromosome
            "-eC", "1" // End chromosome
        };

        DiscoverySNPCallerPlugin plugin = new DiscoverySNPCallerPlugin();
        plugin.setParameters(MDPLowVolArgs);
        plugin.performFunction(null);
    }
    
    public static void runMergeDuplicateSNPsPlugin() {
        String Dir = "";
        
        String[] Args = new String[] {
            ""
        };

        MergeDuplicateSNPsPlugin plugin = new MergeDuplicateSNPsPlugin();
        plugin.setParameters(Args);
        plugin.performFunction(null);
    }
    
    public static void runExtractHapmapSubsetPlugin() {
        String baseDir = "/Volumes/nextgen/Zea/build20120701/06_HapMap/RC2/04_BPECFilteredSNPs/";
        String outDir =  "/Users/jcg233/Documents/GBS/ShilpaNIL28FMJuly2012BuildRC2BPEC/";
//        for (int chr=1; chr<11; chr++) {
//            String[] args = new String[]{
//                "-h", baseDir+"rje22_BPEC_AllZea_GBS_Build_July_2012_RC-2_chr"+chr+".hmp.txt.gz",
//                "-o", outDir+"ShilpaNIL28FMJuly2012BuildRC2BPECOrig_chr"+chr+".hmp.txt.gz",
//                "-p", outDir+"ShilpaNIL28FMSamples20120701buildUTF8.txt"
//            };
//            ExtractHapmapSubsetPlugin plugin = new ExtractHapmapSubsetPlugin(null);
//            plugin.setParameters(args);
//            plugin.performFunction(null);
//        }

//        baseDir =  "/Users/jcg233/Documents/GBS/ShilpaNIL28FMJuly2012BuildRC2BPEC/";
//        for (int chr=1; chr<11; chr++) {
//            String[] args = new String[]{
//                "-h", baseDir+"ShilpaNIL28FMJuly2012BuildRC2BPECOrig_chr"+chr+".hmp.txt.gz",
//                "-o", baseDir+"ShilpaNIL28FMJuly2012BuildRC2BPECOrig216taxa_chr"+chr+".hmp.txt.gz",
//                "-p", baseDir+"ShilpaNIL28FMSamples20120701build216taxa.txt"
//            };
//            ExtractHapmapSubsetPlugin plugin = new ExtractHapmapSubsetPlugin(null);
//            plugin.setParameters(args);
//            plugin.performFunction(null);
//        }

//        baseDir = "/Volumes/nextgen/Zea/build20120701/06_HapMap/RC2.1/04_BPECFilteredSNPs/";
//        outDir =  "/Users/jcg233/Documents/GBS/ZakFMJuly2012BuildRC2BPEC/";
//        for (int chr=1; chr<11; chr++) {
//            String[] args = new String[]{
//                "-h", baseDir+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_RC-2.1_chr"+chr+".hmp.txt.gz",
//                "-o", outDir+"ZakKRNCulmFMJuly2012BuildRC2-1BPEC_chr"+chr+".hmp.txt.gz",
//                "-p", outDir+"ZakKRNCulmFMSamples01072012Build.txt"
//            };
//            ExtractHapmapSubsetPlugin plugin = new ExtractHapmapSubsetPlugin(null);
//            plugin.setParameters(args);
//            plugin.performFunction(null);
//        }

        baseDir = "/Users/jcg233/Documents/GBS/20120701BuildRC2-1BPEC/";
        outDir =  "/Users/jcg233/Documents/GBS/ZakFMJuly2012BuildRC2BPEC/";
        for (int chr=5; chr<11; chr++) {
            String[] args = new String[]{
                "-h", baseDir+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_RC-2.1_chr"+chr+".hmp.txt.gz",
                "-o", outDir+"ZakKRNCulmFMHighCovBC2S3July2012BuildRC2-1BPEC_chr"+chr+".hmp.txt.gz",
                "-p", outDir+"ZakKRNCulmFMHighCovBC2S3Samples01072012Build.txt",
                "-a", "2" // at least 2 "alleles" (actually, genotypes) = polymorphic
            };
            ExtractHapmapSubsetPlugin plugin = new ExtractHapmapSubsetPlugin(null);
            plugin.setParameters(args);
            plugin.performFunction(null);
        }
    }
    
    public static void runCompareGenosBetweenHapMapFilesPlugin() {
        String baseDir = "/Users/jcg233/Documents/GBS/MDP1_low_vol/";
        String[] tassel4WRefVsTassel3Args = new String[]{
            "-hmp1", baseDir+"tassel4/hapmap/withRef/MDP1_low_vol_mergeSNPs.c+.hmp.txt",
//            "-hmp2", baseDir+"hapmap/MDP1_low_vol_noRefOption_mergeSNPs.c+.hmp.txt",
            "-hmp2", baseDir+"hapmap/MDP1_low_vol_RefOptionWOutput2_mergeSNPs.c+.hmp.txt",
            "-sC",   "10",
            "-eC",   "10",
            "-syn",  baseDir+"282synsWRef.txt",
            "-o",    baseDir+"tassel4WRefVsTassel3WRefGenoCompare.txt",
        };

        String prodTestDir = "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/v2.6/";
        String[] prodTestArgs = new String[]{
            "-hmp1", prodTestDir+"ProductionTest/tassel3_RawReadsToHapMap/C08L7ACXX_6_c+.hmp.txt",
            "-hmp2", prodTestDir+"02_MergeDupSNPs/AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr+.hmp.txt.gz",
            "-sC",   "1",
            "-eC",   "10",
            "-syn",  prodTestDir+"ProductionTest/tassel3_RawReadsToHapMap/C08L7ACXX_6_c1_Synonyms.txt",
            "-o",    prodTestDir+"ProductionTest/tassel3_RawReadsToHapMap/productionTestMDPLowVolGenoCompare.txt",
        };

        String tassel3v4ProdDir = "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/v2.6/ProductionTest/";
        String[] tassel3v4ProdArgs = new String[]{
            "-hmp1", tassel3v4ProdDir+"tassel3_RawReadsToHapMap/C08L7ACXX_6_c+.hmp.txt",
            "-hmp2", tassel3v4ProdDir+"tassel4_RawReadsToHapMap/C08L7ACXX_6_c+.hmp.txt",
            "-sC",   "1",
            "-eC",   "10",
            "-syn",  tassel3v4ProdDir+"tassel4_RawReadsToHapMap/C08L7ACXX_6_SynonymsProd.txt",
            "-o",    tassel3v4ProdDir+"tassel4_RawReadsToHapMap/tassel3v4ProdMDPLowVolGenoCompare.txt",
        };

        String rawReadsVsSeqToGenosDir = "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/v2.6/ProductionTest/";
        String[] rawReadsVsSeqToGenosArgs = new String[]{
            "-hmp1", rawReadsVsSeqToGenosDir+"tassel4_RawReadsToHapMap/C08L7ACXX_6_c+.hmp.txt",
            "-hmp2", rawReadsVsSeqToGenosDir+"tassel4_SeqToGenos/C08L7ACXX_6_c+.hmp.txt",
            "-sC",   "1",
            "-eC",   "10",
            "-syn",  rawReadsVsSeqToGenosDir+"tassel4_SeqToGenos/C08L7ACXX_6_SynonymsProd.txt",
            "-o",    rawReadsVsSeqToGenosDir+"tassel4_SeqToGenos/rawReadsVsSeqToGenosMDPLowVolGenoCompare.txt",
        };

        String QualVsQuantDir = "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/v2.6/ProductionTest/";
        String[] QualVsQuantArgs = new String[]{
            "-hmp1", QualVsQuantDir+"tassel4_SeqToGenos/C08L7ACXX_6_c+.hmp.txt",
            "-hmp2", QualVsQuantDir+"tassel4_SeqToGenos/quant/MGP1_low_vol_2reps_chr+.hmp.txt.gz",
            "-sC",   "1",
            "-eC",   "10",
            "-syn",  QualVsQuantDir+"tassel4_SeqToGenos/C08L7ACXX_6_XXvZZSynonymsProd.txt",
            "-o",    QualVsQuantDir+"tassel4_SeqToGenos/quant/QualVsQuantGenoCompare.txt",
        };
        
        String QuantVsOneHapMapDir = "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/v2.6/ProductionTest/tassel4_SeqToGenos/";
        String[] QuantVsOneHapMapArgs = new String[]{
            "-hmp1", QuantVsOneHapMapDir+"quant/testFastq2/MGP1_low_vol_2smallReps_chr+.hmp.txt.gz",
            "-hmp2", QuantVsOneHapMapDir+  "quant/oneAlign/MGP1_low_vol_2smallReps_chr+.hmp.txt.gz",
            "-sC",   "1",
            "-eC",   "10",
            "-syn",  QuantVsOneHapMapDir+"C08L7ACXX_6_XXvZZSynonymsProd.txt",
            "-o",    QuantVsOneHapMapDir+"quant/oneAlign/Quant2VsOneAlignGenoCompare.txt",
        };

        TasselPrefs.putAlignmentRetainRareAlleles(false);
        CompareGenosBetweenHapMapFilesPlugin plugin = new CompareGenosBetweenHapMapFilesPlugin();
        plugin.setParameters(QuantVsOneHapMapArgs);
        plugin.performFunction(null);
    }
    
    public static void expandVariantsInTOPM() {
        String inTOPMFile =  "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/topm/AllZeaMasterTags_c10_20120703.topm";
        String outTOPMFile = "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/topm/AllZeaMasterTags_c10_16vars_20120703.topm";
        int newMaxVariants = 16;
        boolean loadBinary = (inTOPMFile.endsWith(".txt")) ? false : true;
        TagsOnPhysicalMap theTOPM = new TagsOnPhysicalMap(inTOPMFile, loadBinary);
        theTOPM.expandMaxVariants(newMaxVariants);
        theTOPM.writeBinaryFile(new File(outTOPMFile));
    }
    
    public static void getChrsFromTOPM() {
        String inTOPMFile =  "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/topm/AllZeaMasterTags_c10_20120703.topm";
        boolean loadBinary = (inTOPMFile.endsWith(".txt")) ? false : true;
        TagsOnPhysicalMap theTOPM = new TagsOnPhysicalMap(inTOPMFile, loadBinary);
        int[] chrs = theTOPM.getChromosomes();
        Arrays.sort(chrs);
        System.out.println("The chromosomes in this TOPM are:");
        for (int c=0; c<chrs.length; c++) {
            System.out.println(chrs[c]+"");
        }
    }
    
    public static void runTOPMSummaryPlugin() {
        String baseDir = "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/topm/";
        String[] TOPMSummaryArgs = new String[]{
            "-input", baseDir+"AllZeaMasterTags_c10_16vars_UnfiltProdTOPM_chr4_20120703.topm",
            "-output", baseDir+"AllZeaGBSBuild_v2.6_VariantModifiedTOPM_TOPM_Summary_chr4.txt",
        };

        TOPMSummaryPlugin plugin = new TOPMSummaryPlugin(null);
        plugin.setParameters(TOPMSummaryArgs);
        plugin.performFunction(null);
    }
    
    public static void convertTOPMtoHDF5() {
        String inTOPMFileStr =  "/Users/jcg233/largeFiles/topm/AllZeaGBS_v2.6_MergedUnfiltProdTOPM_20130425.topm";
        String outTOPMFileStr = "/Users/jcg233/largeFiles/topm/AllZeaGBSv2.6ProdTOPM_20130605.topm.h5";
        AbstractTagsOnPhysicalMap inTOPM=(AbstractTagsOnPhysicalMap)TOPMUtils.readTOPM(inTOPMFileStr);
        int maxMapping = 4;
        int maxVariants = 16;
        TagsOnPhysMapHDF5.createFile(inTOPM, outTOPMFileStr, maxMapping, maxVariants);
    }
    
    public static void runProductionSNPCallerPlugin() {
        String baseDir = "/Users/jcg233/Documents/GBS/";
        String[] MDPLowVolArgs = new String[]{
            "-i", "/Users/jcg233/largeFiles/testFastq/",
            "-k", baseDir+"MDP1_low_vol/MGP1_low_vol_2smallReps_key.txt",
            "-m", "/Users/jcg233/largeFiles/topm/AllZeaGBSv2.6ProdTOPM_20130605.topm.h5",
            "-e", "ApeKI",
//            "-vL", // VCF likelihood-based calling of hets
            "-o", baseDir+"AllZeaBuild2.X/v2.6/ProductionTest/tassel4_SeqToGenos/lastIndexOf",
        };

        String[] AnnaObrienArgs = new String[]{
            "-i", baseDir+"AllZeaBuild2.X/v2.7/AnnaOBrien/fastq",
            "-k", baseDir+"AllZeaBuild2.X/v2.7/AnnaOBrien/C08L7ACXX_5_key.txt",
            "-m", "/Users/jcg233/largeFiles/topm/AllZeaGBSv2.6ProdTOPM_20130605.topm.h5",
            "-e", "ApeKI",
//            "-vL", // VCF likelihood-based calling of hets
            "-o", baseDir+"AllZeaBuild2.X/v2.7/AnnaOBrien/genos",
        };

        String[] DallasArgs = new String[]{
            "-i", "/Volumes/group/E/SoftwareTesting/GBS/ProductionPipeline/MGP1LowVol2SmallReps/testFastq",
            "-k", "/Volumes/group/E/SoftwareTesting/GBS/ProductionPipeline/MGP1LowVol2SmallReps/MGP1_low_vol_2smallReps_key.txt",
            "-m", "/Volumes/nextgen/Zea/AllZeaBuild_2.X/04_TOPM/2.6_production/02_MergedTOPM/AllZeaGBSv2.6ProdTOPM_20130605.topm.h5",
            "-e", "ApeKI",
//            "-vL", // VCF likelihood-based calling of hets
            "-o", baseDir+"AllZeaBuild2.X/v2.6/ProductionTest/tassel4_SeqToGenos/Dallas",
        };

        String[] MooseP1Args = new String[]{
            "-i", "/Users/jcg233/largeFiles/fastq/Moose",
            "-k", "/Users/jcg233/largeFiles/fastq/Moose/C12B2ACXX_5_key.txt",
            "-m", "/Users/jcg233/largeFiles/topm/AllZeaGBSv2.6ProdTOPM_20130605.topm.h5",
            "-e", "ApeKI",
//            "-vL", // VCF likelihood-based calling of hets
            "-o", baseDir+"AllZeaBuild2.X/v2.7/MooseP1",
        };

        ProductionSNPCallerPlugin plugin = new ProductionSNPCallerPlugin(null);
        plugin.setParameters(MooseP1Args);
        plugin.performFunction(null);
    }
    
    public static void convertHDF5ToHapMap() {
        String hdf5File = "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/v2.6/ProductionTest/tassel4_SeqToGenos/lastIndexOf/MGP1_low_vol_2smallReps.hmp.h5";
        String hmpFile  = "/Users/jcg233/Documents/GBS/AllZeaBuild2.X/v2.6/ProductionTest/tassel4_SeqToGenos/lastIndexOf/MGP1_low_vol_2smallRepsCompare2.hmp.txt.gz";
        System.out.println("Converting from *.hmp.h5 to *.hmp.txt.gz");
        System.out.println("   Reading hmp.h5 file:     "+hdf5File);
        Alignment a = ImportUtils.readGuessFormat(hdf5File, false);
        System.out.println("   Writing hmp.txt.gz file: "+hmpFile);
        if (a == null) {
            System.out.println("Better luck next time!");
        } else {
            ExportUtils.writeToHapmap(a, false, hmpFile, '\t', null);
        }
        System.out.println("Done!");
    }
    
    public static void testAllZeaGBSv27Genos() {
        // compare the h5 output to the hmp.txt.gz output
        String baseDir="/Users/jcg233/Documents/GBS/AllZeaBuild2.X/v2.7/genos/part";
        String genoBase="AllZeaGBS_v2.7_SeqToGenos_part";
        int nParts = 28;
        for (int part = 1; part<=nParts; part++) {
            String partS = part<10? "0"+part : ""+part;
            String hmpFile = baseDir+partS+"/"+genoBase+partS+".hmp.txt.gz";
            Alignment hmp = ImportUtils.readFromHapmap(hmpFile, null);
            String hdf5File = baseDir+partS+"/"+genoBase+partS+".hmp.h5";
            Alignment hdf5 = ImportUtils.readGuessFormat(hdf5File, false);
 //           AlignmentTestingUtils.alignmentsEqual();
            
        }
    }
}
