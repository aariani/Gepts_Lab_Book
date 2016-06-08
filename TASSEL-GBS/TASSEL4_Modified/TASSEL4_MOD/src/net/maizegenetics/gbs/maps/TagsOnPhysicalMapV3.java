/*
 * TagsOnPhysMapHDF5
 */
package net.maizegenetics.gbs.maps;

import ch.systemsx.cisd.hdf5.*;
import java.io.File;
import java.util.Arrays;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.util.GBSHDF5Constants;
import org.apache.log4j.Logger;

/**
 * HDF5 version of TagsOnPhysical Map.  This is the preferred version of physical map as it uses less
 * memory, loads faster, and is more flexible with mapping positions.
 * <p>
 * Multiple mapping positions can be stored for each Tag.  For example, separate aligners could record their
 * positions in the {@link TagMappingInfoV3} objects.  Then the genetic mapping algorithm could be used to resolve,
 * which is the true mapping position.  MapPosition0 is used as the best mapping position, and used by the SNP caller.
 * The fields in {@link TagMappingInfoV3} may change for each aligner and genetic mapping, since part of the field might be missing different aligner and genetic mapping
 * <p>
 * TODO: createFile that includes a Locus filter, only exports positions on the same locus and position range
 * TODO: Resort map positions by quality (by model training)
 * 
 * 
 * @author Ed Buckler, Terry Casstevens, Fei Lu
 */
public class TagsOnPhysicalMapV3 extends AbstractTagsOnPhysicalMap implements TOPMInterface {

    private static final Logger myLogger = Logger.getLogger(TagsOnPhysicalMapV3.class);
    /**Shift 2^16*/
    private static final int BITS_TO_SHIFT_FOR_CHUNK = 16;
    /**65536 tags in a truck*/
    private static final int CHUNK_SIZE = 1 << BITS_TO_SHIFT_FOR_CHUNK;
    /**gzip format compression for TagMappingInfo*/
    private static HDF5GenericStorageFeatures genoFeatures = HDF5GenericStorageFeatures.createDeflation(5); //used by mapping object
    /**gzip format compression for other datasets*/
    private static HDF5IntStorageFeatures vectorFeatures = HDF5IntStorageFeatures.createDeflation(5); //used by vectors
    
    /**Number of physical positions from different aligner or aligner with different parameters*/
    protected int mappingNum = 0;
    /**Writer/Reader of TagsOnPhysicalMapV3*/
    private IHDF5Writer myHDF5 = null;
    /**Tag index in block*/
    private int cachedTagIndex = -1;
    private int cachedMapIndex = -1;
    private TagMappingInfoV3 cachedTMI = null;
    /**Block index where the tag belongs to. Max: TagCount>>BITS_TO_SHIFT_FOR_CHUNK+1*/
    private int cachedChunkIndex = -1;
    private String[] mapNames = null;
    private TagMappingInfoV3[][] cachedTMIChunk = null;
    private int chunkStartTagIndex;
    private int chunkEndTagIndex;
    private boolean cleanMap = true;
    private boolean cacheAllMappingBlocks = false;
    private HDF5CompoundType<TagMappingInfoV3> tmiType = null;
    private boolean hasDetailedMapping=false;

     
    /**
     * Initialize HDF5 TOPM from TagCounts file. The "MAXMAPPING"(set to 0), "TAGLENGTHINLONG", "tags" and "tagLength" are added.
     * @param inTags
     * @param newHDF5file 
     */
    public static void createFile (Tags inTags, String newHDF5file) {
        int tagLengthInLong = inTags.getTagSizeInLong();
        int tagCount = inTags.getTagCount();
        long[][] tags = new long[tagLengthInLong][tagCount];
        byte[] tagLength = new byte[tagCount];
        for (int i = 0; i < tagCount; i++) {
            long[] ct = inTags.getTag(i);
            for (int j = 0; j < tagLengthInLong; j++) {
                tags[j][i] = ct[j];
            }
            tagLength[i] = (byte) inTags.getTagLength(i);
        }
        IHDF5Writer h5 = null;
        try {
            myLogger.info("Creating HDF5 File: " + newHDF5file);
            IHDF5WriterConfigurator config = HDF5Factory.configure(new File(newHDF5file));
            config.overwrite();
            config.useUTF8CharacterEncoding();
            h5 = config.writer();
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXMAPPING, 0);
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGLENGTHINLONG, tagLengthInLong);
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGCOUNT, tagCount);
            h5.createLongMatrix(GBSHDF5Constants.TAGS, inTags.getTagSizeInLong(), tagCount, inTags.getTagSizeInLong(), tagCount, vectorFeatures);
            h5.writeLongMatrix(GBSHDF5Constants.TAGS, tags, vectorFeatures);
            tags=null;
            System.out.println("...Tags written");
            h5.createByteArray(GBSHDF5Constants.TAGLENGTH, tagCount, vectorFeatures);
            h5.writeByteArray(GBSHDF5Constants.TAGLENGTH, tagLength, vectorFeatures);
            tagLength=null;
            System.out.println("...Tags lengths written");
            System.gc();
        }
        catch (Exception e) {
            e.printStackTrace();
            h5.close();
            System.exit(1);
        }
    }
    
    /**
     * Constructor from a HDF5 TOPM file
     * @param theHDF5file 
     */
    public TagsOnPhysicalMapV3 (String theHDF5file) {
        System.out.println("Opening :" + theHDF5file);
        myHDF5 = HDF5Factory.open(theHDF5file);
        tmiType = myHDF5.compounds().getInferredType(TagMappingInfoV3.class);
        this.tagLengthInLong = myHDF5.getIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGLENGTHINLONG);
        this.myNumTags = myHDF5.getIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGCOUNT);
        this.tags = myHDF5.readLongMatrix(GBSHDF5Constants.TAGS);
        this.tagLength = myHDF5.readByteArray(GBSHDF5Constants.TAGLENGTH);
        this.mappingNum = myHDF5.getIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXMAPPING);
        if (mappingNum != 0) {
            this.renameMapNames();
        }
        if(myHDF5.exists(GBSHDF5Constants.BEST_STRAND)) {
            bestStrand=myHDF5.readByteArray(GBSHDF5Constants.BEST_STRAND);
            bestChr= myHDF5.readIntArray(GBSHDF5Constants.BEST_CHR);
            bestStartPos= myHDF5.readIntArray(GBSHDF5Constants.BEST_STARTPOS);  
        }
        if (myHDF5.exists(GBSHDF5Constants.MULTIMAPS)) {
            multimaps = myHDF5.readByteArray(GBSHDF5Constants.MULTIMAPS);
        }
        if(myHDF5.exists(GBSHDF5Constants.VARIANTDEF)) {
            loadVariantsIntoMemory();
        }
        if(myHDF5.exists(GBSHDF5Constants.BEST_STRAND)) {
            this.populateChrAndVarPositions();
        }
        System.gc();
    }
    
    /**
     * Creat datasets in HDF5 holding mapping information, which is used to annotate the TOPM with multiple alignment hypothesis
     * @param startIndex Start index of tag mapping information. This is essentially the current mappingNum
     * @param size the number of datasets which will be created
     * @return names of the datasets
     */
    public String[] creatTagMappingInfoDatasets (int startIndex, int size) {
        int chunkCount = this.getChunkCount(); 
        int chunkSize = this.getChunkSize();
        String[] dataSetNames = new String[size];
        for (int i = 0; i < size; i++) {
            dataSetNames[i] = GBSHDF5Constants.MAPBASE + this.getThreeFigureString(i+startIndex);
            myHDF5.compounds().createArray(dataSetNames[i], tmiType, chunkSize*chunkCount, chunkSize, genoFeatures);
            
        }
        System.out.println("Created new TagMappingInfo datasets. They are");
        for (int i = 0; i < dataSetNames.length; i++) {
            System.out.println(dataSetNames[i]);
        }
        return dataSetNames;
    }
    
    /**
     * Write TMI buffer/chunk to HDF5 datasets, which is used to annotate the TOPM with multiple alignment hypothesis
     * @param dataSetNames
     * @param tmiChunk TMI chunk [dataSetNames.length]*[chunk_size]
     * @param chunkIndex index of this chunk
     */
    public void writeTagMappingInfoDataSets (String[] dataSetNames, TagMappingInfoV3[][] tmiChunk, int chunkIndex) {
        for (int i = 0; i < dataSetNames.length; i++) {   
            myHDF5.compounds().writeArrayBlock(dataSetNames[i], tmiType, tmiChunk[i], chunkIndex);
        }
 
    }
    
    /**
     * Set mappingNum attribute in HDF5
     * @param maxMapping 
     */
    public void setMappingNum (int mappingNum) {
        this.mappingNum = mappingNum;
        myHDF5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXMAPPING, mappingNum);
        this.renameMapNames();
        System.out.println("TOPM maxMapping attibute was set to " + String.valueOf(mappingNum));
    }
    
    /**
     * Rename the mapNames based on the number of mapping
     */
    private void renameMapNames () {
        if (mappingNum == 0) {
            mapNames = null;
            return;
        }
        this.mapNames = new String[mappingNum];
        for (int i = 0; i < mappingNum; i++) {
            mapNames[i] = GBSHDF5Constants.MAPBASE + this.getThreeFigureString(i);
        }
    }
    
    private boolean loadVariantsIntoMemory() {
        int howManyDef=0;
        int readBlock=4096*16;
        variantDefs=new byte[myNumTags][];
        variantOffsets=new byte[myNumTags][];
        if(!myHDF5.exists(GBSHDF5Constants.VARIANTDEF)) return false;
        for (int blockStep = 0; blockStep < myNumTags; blockStep+=readBlock) {
            int blockSize=(myNumTags-blockStep<readBlock)?myNumTags-blockStep:readBlock;
            byte[][] vd=myHDF5.readByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTDEF,blockSize,myMaxVariants,blockStep,0);
//            System.out.println(Arrays.toString(vd[0]));
            byte[][] vo=myHDF5.readByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTPOSOFF,blockSize,myMaxVariants,blockStep,0);
            for (int j = 0; j < blockSize; j++) {
                int cnt=0;
                for (byte bs : vd[j]) {if (bs!=TOPMInterface.BYTE_MISSING) cnt++;}
                if(cnt==0) continue;
                byte[] vdReDim=new byte[cnt];
                byte[] voReDim=new byte[cnt];
                for (int i = 0; i < cnt; i++) {
                    vdReDim[i]=vd[j][i];
                    voReDim[i]=vo[j][i];
                    howManyDef++;
                }
                variantDefs[blockStep+j]=vdReDim;
                variantOffsets[blockStep+j]=voReDim;
            }
            
            //byte[] vd=myHDF5.readByteArrayBlockWithOffset(null, i, i)
        }
        System.out.println("Real Variant Defs:"+howManyDef);
        return true;
    }
    
    private static boolean writeVariantsToHDF5(IHDF5Writer aHDF5, AbstractTagsOnPhysicalMap aTOPM) {
        int howManyDef=0;
        int readBlock=4096*16;
        
        int myNumTags=aTOPM.myNumTags;
        int myMaxVariants=aTOPM.myMaxVariants;
        aHDF5.createByteMatrix(GBSHDF5Constants.VARIANTDEF, myNumTags, myMaxVariants);
        aHDF5.createByteMatrix(GBSHDF5Constants.VARIANTPOSOFF, myNumTags, myMaxVariants);
//        variantDefs=new byte[myNumTags][];
//        variantOffsets=new byte[myNumTags][];
        if(!aHDF5.exists(GBSHDF5Constants.VARIANTDEF)) return false;
        byte[][] vd=new byte[readBlock][myMaxVariants];
        byte[][] vo=new byte[readBlock][myMaxVariants];
        for (int blockStep = 0; blockStep < myNumTags; blockStep+=readBlock) {
            int blockSize=(myNumTags-blockStep<readBlock)?myNumTags-blockStep:readBlock;
            vd=new byte[blockSize][myMaxVariants];
            vo=new byte[blockSize][myMaxVariants];
            for (int j = 0; j < blockSize; j++) {
                for (int v = 0; v < vo[0].length; v++) {
                    vd[j][v]=aTOPM.getVariantDef(blockStep+j, v);
                    vo[j][v]=aTOPM.getVariantPosOff(blockStep+j, v);
                }
            }
            aHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTDEF,vd,blockStep,0);
            aHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTPOSOFF,vo,blockStep,0);
            //byte[] vd=myHDF5.readByteArrayBlockWithOffset(null, i, i)
        }
        System.out.println("Real Variant Defs:"+howManyDef);
        return true;
    }
    
    /**
     * Load mapping information in a chunk to memory and reset tag start and end index
     * @param chunkIndex 
     */
    private void cacheMappingInfoChunk (int chunkIndex) {
        this.cachedTMIChunk = new TagMappingInfoV3[this.getMappingNum()][this.getChunkSize()];
        for (int i = 0; i < mappingNum; i++) {
            cachedTMIChunk[i] = myHDF5.compounds().readArrayBlock(mapNames[i], tmiType, this.getChunkSize(), chunkIndex);
        }
        this.cachedChunkIndex = chunkIndex;
        this.chunkStartTagIndex = chunkIndex*this.getChunkSize();
        this.chunkEndTagIndex = chunkStartTagIndex+this.getChunkSize();
        if (chunkEndTagIndex > this.getTagCount()) chunkEndTagIndex = this.getTagCount();
    }
    
    /**
     * Update current cachedTMI
     * @param tagIndex
     * @param mapIndex 
     */
    private void cacheMappingInfo(int tagIndex, int mapIndex) {
        if (tagIndex == cachedTagIndex && mapIndex == cachedMapIndex) {
            return;
        }
        int chunkIndex = tagIndex >> BITS_TO_SHIFT_FOR_CHUNK;
        if (chunkIndex != this.cachedChunkIndex) {
            this.cacheMappingInfoChunk(chunkIndex); 
        }
        cachedTMI = this.cachedTMIChunk[mapIndex][tagIndex%this.getChunkSize()];
        cachedTagIndex = tagIndex;
        cachedMapIndex = mapIndex;
    }

    private void saveCacheBackToFile() {
        int block = cachedTagIndex >> BITS_TO_SHIFT_FOR_CHUNK;
        if (cachedChunkIndex != block) {
            for (int mi = 0; mi < mappingNum; mi++) {
                if (cleanMap == false) {
//f                    myHDF5.compounds().writeArrayBlock(GBSHDF5Constants.MAPBASE+this.getThreeFigureString(mi), tmiType, cachedTMIBlocks[mi], block);
                }
                if (cacheAllMappingBlocks == false) {
                    break;
                }
            }
            if (cleanMap == false) {
                myHDF5.writeByteArray("multimaps", multimaps);
            }  //this could be made more efficient by just writing the block
            cleanMap = true;
        }
    }

    /**
     * Return the total count of chunks
     * @return 
     */
    public int getChunkCount () {
        return (this.getTagCount()>>BITS_TO_SHIFT_FOR_CHUNK) + 1;
    }
    
    /**
     * Return the chunk size (Number of tags in a chunk)
     * @return 
     */
    public int getChunkSize () {
        return this.CHUNK_SIZE;
    }
    
    public void getFileReadyForClosing() {
//        writeCachedVariantDefs();
//        writeCachedVariantOffsets();
//        saveCacheBackToFile();
    }

    /**
     * Return number of mapping result
     * @return 
     */
    public int getMappingNum () {
        return this.mappingNum; 
    }


    @Override
    public int addVariant(int tagIndex, byte offset, byte base) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getDcoP(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedTagIndex != index) {
            //cacheMappingInfo(index);
        }
        return cachedTMI.dcoP;
    }
    
    @Override
    public byte getDivergence(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedTagIndex != index) {
            //cacheMappingInfo(index);
        }
        return cachedTMI.divergence;
    }

    @Override
    public int getEndPosition(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedTagIndex != index) {
            //cacheMappingInfo(index);
        }
        return cachedTMI.endPosition;
    }

    @Override
    public byte getMapP(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedTagIndex != index) {
            //cacheMappingInfo(index);
        }
        return cachedTMI.mapP;
    }

    @Override
    public int[] getPositionArray(int index) {
        if (cachedTagIndex != index) {
            //cacheMappingInfo(index);
        }
        int[] r = {cachedTMI.chromosome, cachedTMI.strand, cachedTMI.startPosition};
        return r;
    }

    @Override
    public int getReadIndexForPositionIndex(int posIndex) {
        return indicesOfSortByPosition[posIndex];
    }
    
    /**
     * Return tag mapping information of a tag in one map
     * @param tagIndex
     * @param mapIndex
     * @return 
     */
    public TagMappingInfoV3 getMappingInfo (int tagIndex, int mapIndex) {
        this.cacheMappingInfo(tagIndex, mapIndex);
        return this.cachedTMI;
    }
    
    private String getThreeFigureString (int number) {
        String s = String.valueOf(number);
        int length = s.length();
        for (int i = 0; i < 3-length; i++) {
            s = "0"+s;
        }
        return s;
    }
    
    @Override
    public int[] getUniquePositions(int chromosome) {
        if(myUniquePositions==null) populateChrAndVarPositions();
        return myUniquePositions[chromosome];
    }

    public void setMultimaps(int index, byte multimaps) {
        this.multimaps[index] = multimaps;
    }

    @Override
    public void setChromoPosition(int index, int chromosome, byte strand, int positionMin, int positionMax) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setDivergence(int index, byte divergence) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setMapP(int index, byte mapP) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setMapP(int index, double mapP) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public synchronized void setVariantDef(int tagIndex, int variantIndex, byte def) {
        myHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTDEF, new byte[][]{{def},}, tagIndex, variantIndex);
        variantDefs[tagIndex][variantIndex]=def;
    }

    @Override
    public synchronized void setVariantPosOff(int tagIndex, int variantIndex, byte offset) {
        myHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTPOSOFF, new byte[][]{{offset},}, tagIndex, variantIndex);
        variantOffsets[tagIndex][variantIndex]=offset;
    }
    
    /**
     * Preferred method for setting variant information
     * @param tagIndex
     * @param defAndOffset Two dimension [0=definition, 1=offset][upto 16 bytes for each SNP]
     */
    public synchronized void setAllVariantInfo(int tagIndex, byte[][] defAndOffset) {
        myHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTDEF, new byte[][]{defAndOffset[0]}, tagIndex, 0);
        myHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTPOSOFF, new byte[][]{defAndOffset[1]}, tagIndex, 0);
        variantDefs[tagIndex]=defAndOffset[1];
        variantOffsets[tagIndex]=defAndOffset[1];
    }
    

    @Override
    public void clearVariants() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
