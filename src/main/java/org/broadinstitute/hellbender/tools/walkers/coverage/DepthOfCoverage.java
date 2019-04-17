package org.broadinstitute.hellbender.tools.walkers.coverage;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.IOException;
import java.nio.file.Path;
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
public class DepthOfCoverage extends LocusWalkerByInterval {
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
     * Base file name about which to create the coverage information
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Base file location to which to write coverage summary information")
    private String baseFileName = null;

    /**
     * Reads with mapping quality values lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", optional = true, minValue = 0, maxValue = Integer.MAX_VALUE)
    private int minMappingQuality = -1;
    /**
     * Reads with mapping quality values higher than this threshold will be skipped. The default value is the largest number that can be represented as an integer by the program.
     */
    @Argument(fullName = "maxMappingQuality", doc = "Maximum mapping quality of reads to count towards depth", optional = true, minValue = 0, maxValue = Integer.MAX_VALUE)
    private int maxMappingQuality = Integer.MAX_VALUE;
    /**
     * Bases with quality scores lower than this threshold will be skipped. This is set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth", optional = true, minValue = 0, maxValue = Byte.MAX_VALUE)
    private byte minBaseQuality = -1;
    /**
     * Bases with quality scores higher than this threshold will be skipped. The default value is the largest number that can be represented as a byte.
     */
    @Argument(fullName = "maxBaseQuality", doc = "Maximum quality of bases to count towards depth", optional = true, minValue = 0, maxValue = Byte.MAX_VALUE)
    private byte maxBaseQuality = Byte.MAX_VALUE;

    @Argument(fullName = "countType", doc = "How should overlapping reads from the same fragment be handled?", optional = true)
    private CoverageUtils.CountPileupType countType = CoverageUtils.CountPileupType.COUNT_READS;

    /**
     * Instead of reporting depth, the program will report the base pileup at each locus
     */
    @Argument(fullName = "printBaseCounts", shortName = "baseCounts", doc = "Add base counts to per-locus output", optional = true)
    private boolean printBaseCounts = false;

    /**
     * Disabling the tabulation of locus statistics (# loci covered by sample by coverage) should speed up processing.
     */
    @Argument(fullName = "omitLocusTable", shortName = "omitLocusTable", doc = "Do not calculate per-sample per-depth counts of loci", optional = true)
    private boolean omitLocusTable = false;

    /**
     * Disabling the tabulation of interval statistics (mean, median, quartiles AND # intervals by sample by coverage) should speed up processing. This option is required in order to use -nt parallelism.
     */
    @Argument(fullName = "omit-interval-statistics", doc = "Do not calculate per-interval statistics", optional = true, mutex = "calculate-coverage-over-genes")
    private boolean omitIntervals = false;
    /**
     * Disabling the tabulation of total coverage at every base should speed up processing.
     */
    @Argument(fullName = "omitDepthOutputAtEachBase", shortName = "omitBaseOutput", doc = "Do not output depth of coverage at each base", optional = true)
    private boolean omitDepthOutput = false;

    /**
     * Specify a RefSeq file for use in aggregating coverage statistics over genes.
     * <p>
     * This argument is incompatible with --calculateCoverageOverGenes and --omitIntervalStatistics. A warning will be logged and no output file will be produced for the gene list if these arguments are enabled together.
     */
    @Argument(fullName = "calculate-coverage-over-genes", shortName = "gene-list", doc = "Calculate coverage statistics over this list of genes", optional = true, mutex = "omit-interval-statistics")
    private List<String> refSeqGeneListFiles = new ArrayList<>();
    //MAJOR TODO, make sure this ends up behaving correctly with provided interval list
    //TODO this should suppport other file output formats
    // TODO this is the right way to define this private List<FeatureInput<RefSeqFeature>> refSeqGeneList = null;
    /**
     * Remove genes from the gene summary output file if all of its exon bases were not completely covered by traversal.
     */
    @Argument(fullName = "omit-genes-not-entirely-covered-by-traversal", doc = "Do not output gene summary if it was not completely covered by traversal intervals", optional = true)
    private boolean omitPartiallyCoveredGenes = false;

    /**
     * Output file format (e.g. csv, table, rtable); defaults to r-readable table.
     */
    @Argument(fullName = "output-format", doc = "The format of the output file", optional = true)
    CoverageOutputWriter.DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat = CoverageOutputWriter.DEPTH_OF_COVERAGE_OUTPUT_FORMAT.RTABLE;


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
    private boolean includeRefNBases = false;
    /**
     * TODO this is a debug option that seems a bit
     * Use this option to calibrate what bins you want before performing full calculations on your data.
     */
    // TODO worry about this feature later
//    @Advanced
//    @Argument(fullName = "printBinEndpointsAndExit", doc = "Print the bin values and exit immediately", optional = true)
//    private boolean printBinEndpointsAndExit = false;
    /**
     * Sets the low-coverage cutoff for granular binning. All loci with depth < START are counted in the first bin.
     */
    @Advanced
    @Argument(fullName = "start", doc = "Starting (left endpoint) for granular binning", optional = true, minValue = 0)
    private int start = 1;
    /**
     * Sets the high-coverage cutoff for granular binning. All loci with depth > STOP are counted in the last bin.
     */
    @Advanced
    @Argument(fullName = "stop", doc = "Ending (right endpoint) for granular binning", optional = true, minValue = 1)
    private int stop = 500;
    /**
     * Sets the number of bins for granular binning
     */
    @Advanced
    @Argument(fullName = "nBins", doc = "Number of bins to use for granular binning", optional = true, minValue = 0, minRecommendedValue = 1)
    private int nBins = 499;

    /**
     * This option simply disables writing separate files for per-sample summary statistics (total, mean, median, quartile coverage per sample). These statistics are still calculated internally, so enabling this option will not improve runtime.
     */
    @Argument(fullName = "omitPerSampleStats", shortName = "omitSampleSummary", doc = "Do not output the summary files per-sample", optional = true)
    private boolean omitSampleSummary = false;
    /**
     * By default, coverage is partitioning by sample, but it can be any combination of sample, readgroup and/or library.
     */
    @Argument(fullName = "partitionType", shortName = "pt", doc = "Partition type for depth of coverage", optional = true)
    private EnumSet<DoCOutputType.Partition> partitionTypes = EnumSet.of(DoCOutputType.Partition.sample);

    /**
     * Consider a spanning deletion as contributing to coverage. Also enables deletion counts in per-base output.
     */
    @Advanced
    @Argument(fullName = "includeDeletions", shortName = "dels", doc = "Include information on deletions", optional = true)
    private boolean includeDeletions = false;

    @Advanced
    @Argument(fullName = "ignoreDeletionSites", doc = "Ignore sites consisting only of deletions", optional = true)
    boolean ignoreDeletionSites = false;

    // We explicitly handle reference N bases in the code
    public boolean includeNs() {
        return true;
    }

    // We want to make sure to still generate coverage information over uncovered bases.
    public boolean emitEmptyLoci() {
        return true;
    }

    /**
     * For summary file outputs, report the percentage of bases covered to an amount equal to or greater than this number  (e.g. % bases >= CT for each sample). Defaults to 15; can take multiple arguments.
     */
    @Advanced
    @Argument(fullName = "summaryCoverageThreshold", shortName = "ct", doc = "Coverage threshold (in percent) for summarizing statistics", optional = true)
    private List<Integer> coverageThresholds = new ArrayList<>(Collections.singleton(15));

    //    String separator = "\t";
    Map<DoCOutputType.Partition, List<String>> orderCheck = new HashMap<DoCOutputType.Partition, List<String>>();

    //TODO comment this
    //Map of the running intervals to be computed over
    private DepthOfCoveragePartitionedDataStore coverageTotalsForEntireTraversal;
    private Map<DoCOutputType.Partition, List<String>> globalIdentifierMap;

    @Override
    public boolean includeDeletions() {
        return includeDeletions && ! ignoreDeletionSites;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>(2);
        defaultFilters.add(new WellformedReadFilter());
        defaultFilters.add(new ReadFilterLibrary.NotDuplicateReadFilter());
        defaultFilters.add(new ReadFilterLibrary.PrimaryLineReadFilter());
        defaultFilters.add(new ReadFilterLibrary.MappedReadFilter());
        return defaultFilters;
    }


    @Override
    public void onTraversalStart() {

        try {
            writer = new CoverageOutputWriter(outputFormat,
                    partitionTypes,
                    baseFileName,
                    !refSeqGeneListFiles.isEmpty(),
                    printBaseCounts,
                    omitDepthOutput,
                    omitIntervals,
                    omitSampleSummary,
                    omitLocusTable);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Couldn't create output file, encountered exception: " + e.getMessage(), e);
        }

        globalIdentifierMap = makeGlobalIdentifierMap(partitionTypes);
        writer.writeCoverageOutputHeaders(globalIdentifierMap, getSamplesByPartitionFromReadHeader(partitionTypes), coverageThresholds);

        ReadUtils.getSamplesFromHeader(getHeaderForReads());

        // TODO this will probably end up getting
        for (DoCOutputType.Partition type : partitionTypes) {
            orderCheck.put(type, new ArrayList<String>());
            for (String id : getSamplesByPartitionFromReadHeader(type)) {
                orderCheck.get(type).add(id);
            }
            Collections.sort(orderCheck.get(type));
        }

        coverageTotalsForEntireTraversal = createCoveragePartitioner();
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext, Set<Locatable> overlappingIntervals) {
        // TODO evaluate consequences of supporting nonexistant references
        if (includeRefNBases || (hasReference() && BaseUtils.isRegularBase(referenceContext.getBase()))) {
            final Map<DoCOutputType.Partition, Map<String, int[]>> countsByPartition = CoverageUtils.getBaseCountsByPartition(alignmentContext, minMappingQuality, maxMappingQuality, minBaseQuality, maxBaseQuality, countType, partitionTypes, getHeaderForReads());

            if (!omitDepthOutput) {
                writer.writePerLocusDepthSummary(referenceContext.getInterval(), countsByPartition, globalIdentifierMap, includeDeletions);
            }

            // Update the traversing partitioners with this locus data:
            coverageTotalsForEntireTraversal.update(countsByPartition);

            for (Locatable loc : activeCoveragePartitioners.keySet()) {
                // For genes, we don't want to update the interval for non-exon bases
                if (loc.contains(alignmentContext)) {
                    activeCoveragePartitioners.get(loc).update(countsByPartition);
                }
            }


            //TODO figure out what I want to do about this... it appears the old behavior
//            for ( Pair<SimpleInterval, DepthOfCoveragePartitionedDataStore> targetStats : statsByTarget ) {
//                List<String> genes = getGeneNames(targetStats.first,refseqIterator);
//                for (String gene : genes) {
//                    if ( geneNamesToStats.keySet().contains(gene) ) {
//                        logger.debug("Merging "+geneNamesToStats.get(gene).toString()+" and "+targetStats.second.getCoverageByAggregationType(DoCOutputType.Partition.sample).toString());
//                        geneNamesToStats.get(gene).merge(targetStats.second.getCoverageByAggregationType(DoCOutputType.Partition.sample));
//                    } else {
//                        DepthOfCoverageStats merger = new DepthOfCoverageStats(targetStats.second.getCoverageByAggregationType(DoCOutputType.Partition.sample));
//                        geneNamesToStats.put(gene,merger);
//                        statsByGene.add(new Pair<String,DepthOfCoverageStats>(gene,merger));
//                    }
//                }
//            }
        }
    }


    private Map<Locatable, DepthOfCoveragePartitionedDataStore> activeCoveragePartitioners = new HashMap<>();
    private Map<DoCOutputType.Partition, int[][]> statisticsAggrigationByPartitioning = new HashMap<>();

    public void onIntervalStart(Locatable activeInterval) {
        DepthOfCoveragePartitionedDataStore partitioner = createCoveragePartitioner();
        //NOTE: we don't populate statisticsAggrigationByPartitioning here because that gets populated on the fly later
        activeCoveragePartitioners.put(activeInterval, partitioner);
    }

    // TODO this should have logic for handling if the gene didn't get completely covered

    /**
     * When this method is called, we expect at some point in the past onIntervalStart() was called on the same Locatable
     * object and that we have been continuously tracking and adding per-base information for each base in the pileup
     * that overlaps those objects. Given the interval being closed out, we perform one of the following tasks:
     *
     *  - Write out per-interval summary information to each of the partition interval summary files and update the target
     *    summary table corresponding to each track.
     *
     *  - If it corresponded to a gene interval, write out to the gene summary tables the line corresponding the coverage
     *    over targets contained within the genes exome.
     *
     * @param activeInterval
     */
    public void onIntervalEnd(Locatable activeInterval) {
        final DepthOfCoveragePartitionedDataStore partitionerToRemove = activeCoveragePartitioners.remove(activeInterval);
        if (activeInterval instanceof SimpleInterval) {
            // For each partition type we are managing, make sure we update the active statistics table
            for (DoCOutputType.Partition p : partitionTypes) {
                // Write the per-interval depth information as necessary
                final DepthOfCoverageStats coverageByAggregationPartitionType = partitionerToRemove.getCoverageByAggregationType(p);
                writer.writePerIntervalDepthInformation(p, (SimpleInterval) activeInterval, coverageByAggregationPartitionType);

                // Create a new table if necessary
                if (!statisticsAggrigationByPartitioning.containsKey(p)) {
                    statisticsAggrigationByPartitioning.put(p, new int[coverageByAggregationPartitionType.getHistograms().size()][coverageByAggregationPartitionType.getEndpoints().length + 1]);
                }
                // TODO this MUST be switched to a per-locus addition
                CoverageUtils.updateTargetTable(statisticsAggrigationByPartitioning.get(p), coverageByAggregationPartitionType);
            }


        } else if (activeInterval instanceof RefSeqFeature) {
            DepthOfCoverageStats coverageBySample = partitionerToRemove.getCoverageByAggregationType(DoCOutputType.Partition.sample);

            if ( ! omitPartiallyCoveredGenes || ((RefSeqFeature)activeInterval).getTotalExonLength() >= coverageBySample.getNumberOfLociCovered()) {
                writer.writePerGeneDepthInformation((RefSeqFeature) activeInterval, coverageBySample);
            }
        } else {
            throw new GATKException("Unrecognized Locatable object supplied for traversal, only RefSeqFeature and SimpleInterval are supported: "+activeInterval.toString());
        }
    }

    /**
     * When we finish traversal we should have the following held in state:
     *  1. {@link #statisticsAggrigationByPartitioning} that contains a map from partition to a doubly indexed integer array
     *     that corresponds to a count of the number of data sources with at least the given depth. Importantly this is
     *     expected to have had CoverageUtils.updateTargetTable() called on it once per interval (i.e. once per locus).
     *
     *  2. {@link #coverageTotalsForEntireTraversal} Which we expect to have had update() called on it exactly once for
     *     every locus traversed by the locus walker. Consequently we expect partitioner to contain partitioned counts
     *     corresponding to the entire window covered by the tool.
     *
     */
    @Override
    public Object onTraversalSuccess() {
        // Write out accumulated interval summary statistics
        for (DoCOutputType.Partition partition : statisticsAggrigationByPartitioning.keySet()) {
            writer.writeOutputIntervalStatistics(partition, statisticsAggrigationByPartitioning.get(partition), DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins));
        }

        if (!omitSampleSummary) {
            logger.info("Outputting summary info");
            for (DoCOutputType.Partition type : partitionTypes) {
                writer.writeCumulativeOutputSummaryFiles(coverageTotalsForEntireTraversal, type);
            }
        }

        if (!omitLocusTable) {
            logger.info("Outputting locus summary");
            for (DoCOutputType.Partition type : partitionTypes) {
                writer.writePerTraversalLocusFiles(coverageTotalsForEntireTraversal, type);
            }
        }

        writer.close();
        return null;
    }





    private DepthOfCoveragePartitionedDataStore createCoveragePartitioner() {
        DepthOfCoveragePartitionedDataStore aggro = new DepthOfCoveragePartitionedDataStore(partitionTypes, start, stop, nBins);
        for (DoCOutputType.Partition t : partitionTypes) {
            aggro.addIdentifiers(t, getSamplesByPartitionFromReadHeader(t));
        }
        aggro.initialize(includeDeletions, omitLocusTable);
        checkOrder(aggro);
        return aggro;
    }

    private LinkedHashSet<String> getSamplesByPartitionFromReadHeader(Collection<DoCOutputType.Partition> types) {
        LinkedHashSet<String> partitions = new LinkedHashSet<String>(); // since the DOCS object uses a HashMap, this will be in the same order
        for (DoCOutputType.Partition t : types) {
            partitions.addAll(getSamplesByPartitionFromReadHeader(t));
        }
        return partitions;
    }

    /**
     * TODO ensure comments
     * This makes a global sample map by identifier. That is, per every partition type the lists are lexographically
     * sorted internally.
     *
     * @param types Partition types to generate samples for
     * @return A map of {@link DoCOutputType.Partition} to list of sorted samples
     */
    private Map<DoCOutputType.Partition, List<String>> makeGlobalIdentifierMap(Collection<DoCOutputType.Partition> types) {
        Map<DoCOutputType.Partition, List<String>> partitions = new HashMap<>(); // since the DOCS object uses a HashMap, this will be in the same order
        for (DoCOutputType.Partition t : types) {
            List<String> samplesForType = new ArrayList<>(getSamplesByPartitionFromReadHeader(t));
            samplesForType.sort(String::compareTo);
            partitions.put(t, Collections.unmodifiableList(samplesForType));
        }
        return partitions;
    }

    private HashSet<String> getSamplesByPartitionFromReadHeader(DoCOutputType.Partition type) {
        final HashSet<String> partition = new HashSet<String>();
        final SAMFileHeader header = getHeaderForReads();
        if (type == DoCOutputType.Partition.sample) {
            partition.addAll(ReadUtils.getSamplesFromHeader(header));
        } else {
            for (SAMReadGroupRecord rg : header.getReadGroups()) {
                partition.add(CoverageUtils.getTypeID(rg, type));
            }
        }
        return partition;
    }


    /**
     * We want to keep track of per-interval information for both the user specified intervals and for user specified genes
     *
     * @return A combined list of the unmerged user specified intervals and any specified refSeqFeatures specified
     */
    @Override
    public List<Locatable> getIntervalObjectsToQueryOver() {
        List<? extends Locatable> userProvidedIntervals = intervalArgumentCollection.getSpecifiedIntervalsWithoutMerging(getBestAvailableSequenceDictionary());

        final List<Locatable> refSeqInputs = new ArrayList<>(userProvidedIntervals);
        for (String input : refSeqGeneListFiles) {
            FeatureDataSource<RefSeqFeature> refSeqReader = new FeatureDataSource<>(input);
            for (final RefSeqFeature vcfRecord : refSeqReader) {
                refSeqInputs.add(vcfRecord);
            }
        }
        return refSeqInputs;
    }



    //TODO do we really need to do this going forwards?
    //TODO I well understand this is a very crude way of keeping track of this information.... Ideally this should live
    //TODO as a global immutible list that gets passed arround everywhere necessary....
    private void checkOrder(DepthOfCoveragePartitionedDataStore ag) {
        // make sure the ordering stored at initialize() is propagated along reduce
        for (DoCOutputType.Partition t : partitionTypes ) {
            List<String> order = orderCheck.get(t);
            List<String> namesInAg = ag.getIdentifiersByType().get(t);

            // todo -- chris check me
            Set<String> namesInDOCS = ag.getCoverageByAggregationType(t).getAllSamples();
            int index = 0;
            for ( String s : namesInAg ) {
                if ( ! s.equalsIgnoreCase(order.get(index)) ) {
                    throw new GATKException("IDs are out of order for type "+t+"! Aggregator has different ordering");
                }
                index++;
            }
        }
    }
}


/**
 * A class for storing running intervalPartition data.
 *
 * This class is responsible for holding a running
 */
class DepthOfCoveragePartitionedDataStore {
    private Collection<DoCOutputType.Partition> types;
    private Map<DoCOutputType.Partition,DepthOfCoverageStats> coverageProfiles;
    private Map<DoCOutputType.Partition,List<String>> identifiersByType; //TODO do away with this
    private Set<String> allIdentifiers; //TODO do away with this
    public DepthOfCoveragePartitionedDataStore(Collection<DoCOutputType.Partition> typesToUse, int start, int stop, int nBins) {
        coverageProfiles = new TreeMap<>();
        identifiersByType = new HashMap<>();
        types = typesToUse;
        for ( DoCOutputType.Partition type : types ) {
            coverageProfiles.put(type,new DepthOfCoverageStats(DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins)));
            identifiersByType.put(type,new ArrayList<>());
        }
        allIdentifiers = new HashSet<String>();
    }

    public void merge(DepthOfCoveragePartitionedDataStore otherAggregator) {
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