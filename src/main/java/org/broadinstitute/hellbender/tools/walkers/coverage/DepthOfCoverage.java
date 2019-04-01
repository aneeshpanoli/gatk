package org.broadinstitute.hellbender.tools.walkers.coverage;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;


/**
 * TODO
 */

@CommandLineProgramProperties(
        summary = "TODO",
        oneLineSummary = "TODO",
        programGroup = CoverageAnalysisProgramGroup.class)
@BetaFeature
@DocumentedFeature
public class DepthOfCoverage extends LocusWalker {
    private final static Logger logger = Logger.getLogger(DepthOfCoverage.class);
    private CoverageOutputWriter writer;

    /**
     * Warning message for when the incompatible arguments --calculateCoverageOverGenes and --omitIntervalStatistics are used together.
     */
    private static final String incompatibleArgsMsg = "The arguments --calculateCoverageOverGenes and --omitIntervalStatistics are incompatible. Using them together will result in an empty gene summary output file.";
// TODO figure out output multiplexing of files/types
//    @Output
//    @Multiplex(value=DoCOutputMultiplexer.class,arguments={"partitionTypes","refSeqGeneList","omitDepthOutput","omitIntervals","omitSampleSummary","omitLocusTable"})
//    Map<DoCOutputType,PrintStream> out;
    /**
     * Reads with mapping quality values lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", optional = true, minValue = 0, maxValue = Integer.MAX_VALUE)
    int minMappingQuality = -1;
    /**
     * Reads with mapping quality values higher than this threshold will be skipped. The default value is the largest number that can be represented as an integer by the program.
     */
    @Argument(fullName = "maxMappingQuality", doc = "Maximum mapping quality of reads to count towards depth", optional = true, minValue = 0, maxValue = Integer.MAX_VALUE)
    int maxMappingQuality = Integer.MAX_VALUE;
    /**
     * Bases with quality scores lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth", optional = true, minValue = 0, maxValue = Byte.MAX_VALUE)
    byte minBaseQuality = -1;
    /**
     * Bases with quality scores higher than this threshold will be skipped. The default value is the largest number that can be represented as a byte.
     */
    @Argument(fullName = "maxBaseQuality", doc = "Maximum quality of bases to count towards depth", optional = true, minValue = 0, maxValue = Byte.MAX_VALUE)
    byte maxBaseQuality = Byte.MAX_VALUE;

    @Argument(fullName = "countType", doc = "How should overlapping reads from the same fragment be handled?", optional = true)
    CoverageUtils.CountPileupType countType = CoverageUtils.CountPileupType.COUNT_READS;

    /**
     * Instead of reporting depth, the program will report the base pileup at each locus
     */
    @Argument(fullName = "printBaseCounts", shortName = "baseCounts", doc = "Add base counts to per-locus output", optional = true)
    boolean printBaseCounts = false;

    /**
     * Disabling the tabulation of locus statistics (# loci covered by sample by coverage) should speed up processing.
     */
    @Argument(fullName = "omitLocusTable", shortName = "omitLocusTable", doc = "Do not calculate per-sample per-depth counts of loci", optional = true)
    boolean omitLocusTable = false;

    /**
     * Disabling the tabulation of interval statistics (mean, median, quartiles AND # intervals by sample by coverage) should speed up processing. This option is required in order to use -nt parallelism.
     */
    @Argument(fullName = "omit-interval-statistics", doc = "Do not calculate per-interval statistics", optional = true, mutex = "calculate-coverage-over-genes")
    boolean omitIntervals = false;
    /**
     * Disabling the tabulation of total coverage at every base should speed up processing.
     */
    @Argument(fullName = "omitDepthOutputAtEachBase", shortName = "omitBaseOutput", doc = "Do not output depth of coverage at each base", optional = true)
    boolean omitDepthOutput = false;

    /**
     * Specify a RefSeq file for use in aggregating coverage statistics over genes.
     *
     * This argument is incompatible with --calculateCoverageOverGenes and --omitIntervalStatistics. A warning will be logged and no output file will be produced for the gene list if these arguments are enabled together.
     *
     */
    @Argument(fullName = "calculate-coverage-over-genes", shortName = "gene-list", doc = "Calculate coverage statistics over this list of genes", optional = true, mutex = "omit-interval-statistics")
    File refSeqGeneList = null;

    /**
     * Output file format (e.g. csv, table, rtable); defaults to r-readable table.
     */
    @Argument(fullName = "output-format", doc = "The format of the output file", optional = true)
    DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat = DEPTH_OF_COVERAGE_OUTPUT_FORMAT.RTABLE;


    // ---------------------------------------------------------------------------
    //
    // Advanced arguments
    //
    // ---------------------------------------------------------------------------

    /**
     * Normally, sites where the reference is N (or another non-canonical base) are skipped. If this option is enabled, these sites will be included in DoC calculations if there is coverage from neighboring reads.
     */
    @Advanced
    @Argument(fullName = "includeRefNSites", doc = "Include sites where the reference is N", optional = true)
    boolean includeRefNBases = false;
    /**
     * TODO this is a debug option that seems a bit
     * Use this option to calibrate what bins you want before performing full calculations on your data.
     */
    @Advanced
    @Argument(fullName = "printBinEndpointsAndExit", doc = "Print the bin values and exit immediately", optional = true)
    boolean printBinEndpointsAndExit = false;
    /**
     * Sets the low-coverage cutoff for granular binning. All loci with depth < START are counted in the first bin.
     */
    @Advanced
    @Argument(fullName = "start", doc = "Starting (left endpoint) for granular binning", optional = true, minValue = 0)
    int start = 1;
    /**
     * Sets the high-coverage cutoff for granular binning. All loci with depth > STOP are counted in the last bin.
     */
    @Advanced
    @Argument(fullName = "stop", doc = "Ending (right endpoint) for granular binning", optional = true, minValue = 1)
    int stop = 500;
    /**
     * Sets the number of bins for granular binning
     */
    @Advanced
    @Argument(fullName = "nBins", doc = "Number of bins to use for granular binning", optional = true, minValue = 0, minRecommendedValue = 1)
    int nBins = 499;

    /**
     * This option simply disables writing separate files for per-sample summary statistics (total, mean, median, quartile coverage per sample). These statistics are still calculated internally, so enabling this option will not improve runtime.
     */
    @Argument(fullName = "omitPerSampleStats", shortName = "omitSampleSummary", doc = "Do not output the summary files per-sample", optional = true)
    boolean omitSampleSummary = false;
    /**
     * By default, coverage is partitioning by sample, but it can be any combination of sample, readgroup and/or library.
     */
    @Argument(fullName = "partitionType", shortName = "pt", doc = "Partition type for depth of coverage", optional = true)
    Set<DoCOutputType.Partition> partitionTypes = EnumSet.of(DoCOutputType.Partition.sample);

    /**
     * Consider a spanning deletion as contributing to coverage. Also enables deletion counts in per-base output.
     */
    @Advanced
    @Argument(fullName = "includeDeletions", shortName = "dels", doc = "Include information on deletions", optional = true)
    boolean includeDeletions = false;

    @Advanced
    @Argument(fullName = "ignoreDeletionSites", doc = "Ignore sites consisting only of deletions", optional = true)
    boolean ignoreDeletionSites = false;

    /**
     * For summary file outputs, report the percentage of bases covered to an amount equal to or greater than this number  (e.g. % bases >= CT for each sample). Defaults to 15; can take multiple arguments.
     */
    @Advanced
    @Argument(fullName = "summaryCoverageThreshold", shortName = "ct", doc = "Coverage threshold (in percent) for summarizing statistics", optional = true)
    int[] coverageThresholds = {15};

    String[] OUTPUT_FORMATS = {"table","rtable","csv"};
    String separator = "\t";
    Map<DoCOutputType.Partition, List<String>> orderCheck = new HashMap<DoCOutputType.Partition,List<String>>();

    public enum DEPTH_OF_COVERAGE_OUTPUT_FORMAT {
        TABLE,
        RTABLE,
        CSV
    }


    @Override
    public void onTraversalStart() {
        if ( printBinEndpointsAndExit ) {
            int[] endpoints = DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins);
            System.out.print("[ ");
            for ( int e : endpoints ) {
                System.out.print(e+" ");
            }
            System.out.println("]");
            System.exit(0);
        }

        try {
            writer = new FastaReferenceWriterBuilder()
                    .setFastaFile(path)
                    .setBasesPerLine(basesPerLine)
                    .build();
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Couldn't create " + output + ", encountered exception: " + e.getMessage(), e);
        }

//        // Check the output format
//        boolean goodOutputFormat = false;
//        for ( String f : OUTPUT_FORMATS ) {
//            goodOutputFormat = goodOutputFormat || f.equals(outputFormat);
//        }
//
//        if ( ! goodOutputFormat ) {
//            throw new IllegalArgumentException("Improper output format. Can be one of table,rtable,csv. Was "+outputFormat);
//        }
        //TODO write an output parser
        if ( outputFormat.equals("csv") ) {
            separator = ",";
        }

        //TODO this becomes a header writer
        if ( ! omitDepthOutput ) { // print header
            PrintStream out = getCorrectStream(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary);
            out.printf("%s\t%s","Locus","Total_Depth");
            for (DoCOutputType.Partition type : partitionTypes ) {
                out.printf("\t%s_%s","Average_Depth",type.toString());
            }

            // get all the samples
            HashSet<String> allSamples = getSamplesFromToolKit(partitionTypes);
            ArrayList<String> allSampleList = new ArrayList<String>(allSamples.size());
            for ( String s : allSamples ) {
                allSampleList.add(s);
            }
            Collections.sort(allSampleList);

            for ( String s : allSampleList) {
                out.printf("\t%s_%s","Depth_for",s);
                if ( printBaseCounts ) {
                    out.printf("\t%s_%s",s,"base_counts");
                }
            }

            out.printf("%n");

        } else {
            logger.info("Per-Locus Depth of Coverage output was omitted");
        }

        for (DoCOutputType.Partition type : partitionTypes ) {
            orderCheck.put(type,new ArrayList<String>());
            for ( String id : getSamplesFromToolKit(type) ) {
                orderCheck.get(type).add(id);
            }
            Collections.sort(orderCheck.get(type));
        }
    }






    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (includeRefNBases || BaseUtils.isRegularBase(referenceContext.getBase())) {
            Map<DoCOutputType.Partition, Map<String, int[]>> countsByPartition = CoverageUtils.getBaseCountsByPartition(alignmentContext, minMappingQuality, maxMappingQuality, minBaseQuality, maxBaseQuality, countType, partitionTypes);

            if ( ! omi)

        }
    }


    @Override
    public Object onTraversalSuccess() {
        if ( ! omitSampleSummary ) {
            logger.info("Outputting summary info");
            for (DoCOutputType.Partition type : partitionTypes ) {
                outputSummaryFiles(coverageProfiles,type);
            }
        }

        if ( ! omitLocusTable ) {
            logger.info("Outputting locus summary");
            for (DoCOutputType.Partition type : partitionTypes ) {
                outputLocusFiles(coverageProfiles,type);
            }
        }
        return null;
    }

    private void outputLocusFiles(CoveragePartitioner coverageProfiles, DoCOutputType.Partition type ) {
        printPerLocus(getCorrectStream(type, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_counts),
                getCorrectStream(type, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_proportions),
                coverageProfiles.getCoverageByAggregationType(type),type);
    }

    private void outputSummaryFiles(CoveragePartitioner coverageProfiles, DoCOutputType.Partition type ) {
        printPerSample(getCorrectStream(type, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.statistics),coverageProfiles.getCoverageByAggregationType(type));
        printSummary(getCorrectStream(type, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.summary),coverageProfiles.getCoverageByAggregationType(type));
    }








    private HashSet<String> getSamplesFromToolKit( Collection<DoCOutputType.Partition> types ) {
        HashSet<String> partitions = new HashSet<String>(); // since the DOCS object uses a HashMap, this will be in the same order
        for (DoCOutputType.Partition t : types ) {
            partitions.addAll(getSamplesFromToolKit(t));
        }

        return partitions;
    }

    private HashSet<String> getSamplesFromToolKit(DoCOutputType.Partition type) {
        final HashSet<String> partition = new HashSet<String>();
        final SAMFileHeader header = getHeaderForReads();
        if (type == DoCOutputType.Partition.sample) {
            partition.addAll(ReadUtils.getSamplesFromHeader(header));
        } else {
            for ( SAMReadGroupRecord rg : header.getReadGroups() ) {
                partition.add(CoverageUtils.getTypeID(rg, type));
            }
        }
        return  partition;
    }
}



class CoveragePartitioner {
    private Collection<DoCOutputType.Partition> types;
    private Map<DoCOutputType.Partition,DepthOfCoverageStats> coverageProfiles;
    private Map<DoCOutputType.Partition,List<String>> identifiersByType;
    private Set<String> allIdentifiers;
    public CoveragePartitioner(Collection<DoCOutputType.Partition> typesToUse, int start, int stop, int nBins) {
        coverageProfiles = new TreeMap<DoCOutputType.Partition,DepthOfCoverageStats>();
        identifiersByType = new HashMap<DoCOutputType.Partition,List<String>>();
        types = typesToUse;
        for ( DoCOutputType.Partition type : types ) {
            coverageProfiles.put(type,new DepthOfCoverageStats(DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins)));
            identifiersByType.put(type,new ArrayList<String>());
        }
        allIdentifiers = new HashSet<String>();
    }

    public void merge(CoveragePartitioner otherAggregator) {
        for ( DoCOutputType.Partition type : types ) {
            this.coverageProfiles.get(type).merge(otherAggregator.coverageProfiles.get(type));
        }
    }

    public DepthOfCoverageStats getCoverageByAggregationType(DoCOutputType.Partition t) {
        return coverageProfiles.get(t);
    }

    public void addIdentifiers(DoCOutputType.Partition t, Set<String> ids) {
        for ( String s : ids ) {
            coverageProfiles.get(t).addSample(s);
            identifiersByType.get(t).add(s);
            allIdentifiers.add(s);
        }
        Collections.sort(identifiersByType.get(t));
    }

    public void initialize(boolean useDels, boolean omitLocusTable) {
        for ( DoCOutputType.Partition t : types ) {
            if ( useDels ) {
                coverageProfiles.get(t).initializeDeletions();
            }
            if ( ! omitLocusTable ) {
                coverageProfiles.get(t).initializeLocusCounts();
            }
        }
    }

    public void update(Map<DoCOutputType.Partition,Map<String,int[]>> countsByIdentifierByType) {
        for ( DoCOutputType.Partition t : types ) {
            coverageProfiles.get(t).update(countsByIdentifierByType.get(t));
        }
    }

    public Set<String> getAllIdentifiers() {
        return allIdentifiers;
    }

    public Map<DoCOutputType.Partition,List<String>> getIdentifiersByType() {
        return identifiersByType;
    }
}