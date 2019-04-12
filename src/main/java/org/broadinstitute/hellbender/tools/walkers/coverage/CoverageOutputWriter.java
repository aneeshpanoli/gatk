package org.broadinstitute.hellbender.tools.walkers.coverage;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Tuple;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.*;


/**
 * This is a class for managing the output formatting/files for DepthOfCoverage.
 *
 * Output for {@link DepthOfCoverage} can be organized as a combination Partition, Aggrigation, and OutputType, with this writer
 * storing an internal list of streams to which to write output organized by {@link DoCOutputType} objects. Generally
 * speaking there are three patterns that this writer is responsible for managing:
 *
 *  1. writePerLocusDepthSummary() - This should be called over every locus and is responsible for producing the locus
 *                                   coverage output table summarizing the coverage (and possibly base counts) for each
 *                                   for every partition type x relevant samples. This takes as input the output from
 *                                   {@link CoverageUtils#getBaseCountsByPartition}
 *
 *  2a. writePerIntervalDepthInformation() - This should be called once per traversal interval once it is finished and takes
 *                                        as input a {@link DepthOfCoveragePartitioner} object corresponding to the coverage
 *                                        information over the whole interval. This is used to write both "_interval_summary"
 *                                        and "_interval_statistics" files for each partition.
 *
 *  2b. writePerGeneDepthInformation() - Similarly to 2a, this should be called once per gene in the interval traversal and
 *                                       it takes the same gen-interval coverage summary as 2a. TODO Note that this may not always cover the gene at hand because of reasons
 *
 *  3. writeTraversalCumulativeCoverage() - This method takes a {@link DepthOfCoveragePartitioner} object that should correspond
 *                                         to the partitioned counts for every base traversed by DepthOfCoverage aggregated.
 *                                         This method is responsible for outputting the "_cumulative_coverage_counts",
 *                                         "_cumulative_coverage_proportions", "_statistics", and "_summary" file writing.
 *
 */
public class CoverageOutputWriter implements Closeable {

    private final Set<DoCOutputType.Partition> partitions;
    private final Path refSeqGeneList;
    private final String separator;
    private boolean printBaseCounts;
    private Map<DoCOutputType,PrintStream> outputs;
    private DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat;
    private boolean omitDepthOutput;
    private int[] coverageThresholds;

    public enum DEPTH_OF_COVERAGE_OUTPUT_FORMAT {
        TABLE,
        RTABLE,
        CSV
    }

    public CoverageOutputWriter(final DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat,
                                final Set<DoCOutputType.Partition> partitions,
                                final Path refSeqGeneList, //TODO this is a placeholder until we figure out what the hell to do in order to parse these inputs
                                final String outputBaseName,
                                final boolean printBaseCounts,
                                final boolean omitDepthOutput,
                                final boolean omitIntervals,
                                final boolean omitSampleSummary,
                                final boolean omitLocusTable)  throws IOException{
        this.outputFormat = outputFormat;
        this.partitions = partitions;
        this.refSeqGeneList = refSeqGeneList;
        this.printBaseCounts = printBaseCounts;
        this.omitDepthOutput = omitDepthOutput;

        //TODO more elegit
        if (outputFormat == DEPTH_OF_COVERAGE_OUTPUT_FORMAT.CSV) {
            separator = ",";
        } else {
            separator = "\t";
        }

        outputs = new HashMap<>();
        if(!omitDepthOutput) {
            // Create a depth summary output sink
            DoCOutputType depthSummaryByLocus = new DoCOutputType(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary);
            outputs.put(depthSummaryByLocus, getOutputStream(outputBaseName, depthSummaryByLocus));
        }

        if(!omitIntervals) {
            for(DoCOutputType.Partition partition: partitions) {
                // Create an interval summary output sink
                DoCOutputType intervalSummaryBypartition = new DoCOutputType(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.summary);
                outputs.put(intervalSummaryBypartition, getOutputStream(outputBaseName, intervalSummaryBypartition));

                // Create an interval summary output sink
                DoCOutputType intervalStatisticsByPartition = new DoCOutputType(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.statistics);
                outputs.put(intervalStatisticsByPartition, getOutputStream(outputBaseName, intervalStatisticsByPartition));
            }
        }

        if(refSeqGeneList != null && partitions.contains(DoCOutputType.Partition.sample)) {
            DoCOutputType geneSummaryOut = new DoCOutputType(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
            outputs.put(geneSummaryOut, getOutputStream(outputBaseName, geneSummaryOut));
        }

        if(!omitSampleSummary) {
            for(DoCOutputType.Partition partition: partitions) {
                DoCOutputType cumulativeSummaryOut = new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.summary);
                outputs.put(cumulativeSummaryOut, getOutputStream(outputBaseName, cumulativeSummaryOut));


                DoCOutputType cumulativeStatisticsOut = new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.statistics);
                outputs.put(cumulativeStatisticsOut, getOutputStream(outputBaseName, cumulativeStatisticsOut));
            }
        }

        if(!omitLocusTable) {
            for(DoCOutputType.Partition partition: partitions) {
                DoCOutputType cumulativeCoverageCountsOut = new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_counts);
                outputs.put(cumulativeCoverageCountsOut, getOutputStream(outputBaseName, cumulativeCoverageCountsOut));


                DoCOutputType cumulativeCoverageProportionsOut = new DoCOutputType(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_proportions);
                outputs.put(cumulativeCoverageProportionsOut, getOutputStream(outputBaseName, cumulativeCoverageProportionsOut));
            }
        }
    }

    // Helper method to generate an output stream given a DoCOutputType Object and the base name
    protected static PrintStream getOutputStream(String outputBaseName, DoCOutputType depthSummaryByLocus) throws IOException {
        return new PrintStream(new BufferedOutputStream(Files.newOutputStream(IOUtils.getPath(depthSummaryByLocus.getFilePath(outputBaseName)))));
    }

    //TODO wildly clean every part of this up incredibly
    public void writeCoverageOutputHeaders(final LinkedHashSet<String> allSamples, final int[] coverageThresholds) {
        this.coverageThresholds = coverageThresholds;

        if ( ! omitDepthOutput ) { // print header
            writeDepthOutputSummaryHeader(allSamples);
        }

        // write
        if ( refSeqGeneList != null && partitions.contains(DoCOutputType.Partition.sample) ) {
            PrintStream geneSummaryOut = getCorrectStream(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
            printLocusSummaryHeader("Gene", geneSummaryOut, allSamples);
        }

    }

    protected void writeDepthOutputSummaryHeader(LinkedHashSet<String> allSamples) {
        PrintStream out = getCorrectStream(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary);
        out.printf("%s\t%s","Locus","Total_Depth");
        for (DoCOutputType.Partition type : partitions ) {
            out.printf("\t%s_%s","Average_Depth",type.toString());
        }

        // get all the samples
        ArrayList<String> allSampleList = new ArrayList<>(allSamples);
        Collections.sort(allSampleList);

        for ( String s : allSampleList) {
            out.printf("\t%s_%s","Depth_for",s);
            if ( printBaseCounts ) {
                out.printf("\t%s_%s",s,"base_counts");
            }
        }
        out.printf("%n");
    }


    private PrintStream getCorrectStream(DoCOutputType.Partition partition, DoCOutputType.Aggregation aggregation, DoCOutputType.FileType fileType) {
        DoCOutputType outputType = new DoCOutputType(partition,aggregation,fileType);
        if(!outputs.containsKey(outputType)) {
            throw new UserException(String.format("Unable to find appropriate stream for partition = %s, aggregation = %s, file type = %s", partition, aggregation, fileType));
        }
        return outputs.get(outputType);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outward facing writer methods
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    final static DecimalFormat doubleFormat = new DecimalFormat("#.2#");

    public void writePerLocusDepthSummary(SimpleInterval locus, Map<DoCOutputType.Partition, Map<String,int[]>> countsBySampleByType, Map<DoCOutputType.Partition,List<String>> identifiersByType) {
        PrintStream stream = getCorrectStream(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary);
        // get the depths per sample and build up the output string while tabulating total and average coverage
        StringBuilder perSampleOutput = new StringBuilder();
        int tDepth = 0;
        boolean depthCounted = false;
        for (DoCOutputType.Partition type : partitions ) {
            Map<String,int[]> countsByID = countsBySampleByType.get(type);
            for ( String s : identifiersByType.get(type) ) {
                perSampleOutput.append(separator);
                long dp = (countsByID != null && countsByID.keySet().contains(s)) ? MathUtils.sum(countsByID.get(s)) : 0 ;
                perSampleOutput.append(dp);
                if ( printBaseCounts ) {
                    perSampleOutput.append(separator);
                    perSampleOutput.append(baseCounts(countsByID != null ? countsByID.get(s) : null ));
                }
                if ( ! depthCounted ) {
                    tDepth += dp;
                }
            }
            depthCounted = true; // only sum the total depth once
        }

        // remember -- genome locus was printed in map()
        stream.print(locus);
        stream.print(separator+tDepth);
        for (DoCOutputType.Partition type : partitions ) {
            stream.print(separator + doubleFormat.format((double) tDepth / identifiersByType.get(type).size() ) );
        }
        stream.printf("%s%n",perSampleOutput);
    }


    public void writePerIntervalDepthInformation(SimpleInterval locus, DepthOfCoverageStats intervalStats) {

    }
    public void writePerGeneDepthInformation(Locatable, DepthOfCoverageStats intervalStats) {

    }
    public void writeTraversalCumulativeCoverage(DepthOfCoverageStats) {}


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outward facing writer methods
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    private DepthOfCoverageStats printIntervalStats(List<Tuple<Locatable, DepthOfCoveragePartitioner>> statsByInterval, PrintStream summaryOut, PrintStream statsOut, DoCOutputType.Partition type) {
        Tuple<Locatable, DepthOfCoveragePartitioner> firstPair = statsByInterval.get(0);
        DepthOfCoveragePartitioner firstAggregator = firstPair.b;
        DepthOfCoverageStats firstStats = firstAggregator.getCoverageByAggregationType(type);

        printLocusSummaryHeader(summaryOut, firstStats);

        int[][] nTargetsByAvgCvgBySample = new int[firstStats.getHistograms().size()][firstStats.getEndpoints().length+1];

        for ( Tuple<Locatable, DepthOfCoveragePartitioner> targetAggregator : statsByInterval ) {

            Tuple<Locatable,DepthOfCoverageStats> targetStats = new Tuple<Locatable,DepthOfCoverageStats>(
                    targetAggregator.a, targetAggregator.b.getCoverageByAggregationType(type));
            printTargetSummary(summaryOut,targetStats);
            updateTargetTable(nTargetsByAvgCvgBySample,targetStats.b);
        }

        printIntervalTable(statsOut,nTargetsByAvgCvgBySample,firstStats.getEndpoints());

        return firstStats;
    }

    private void printLocusSummaryHeader(final String title, final PrintStream summaryOut, final List<String> allSamples) {
        StringBuilder summaryHeader = new StringBuilder();
        summaryHeader.append(title);
        summaryHeader.append(separator);
        summaryHeader.append("total_coverage");
        summaryHeader.append(separator);
        summaryHeader.append("average_coverage");

        //TODO this was probably where the order in the header went wrong compared to gatk3
        for ( String s : allSamples ) {
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_total_cvg");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_mean_cvg");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_granular_Q1");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_granular_median");
            summaryHeader.append(separator);
            summaryHeader.append(s);
            summaryHeader.append("_granular_Q3");
            // TODO right here we should really add some space for the new site score variable
            for ( int thresh : coverageThresholds ) {
                summaryHeader.append(separator);
                summaryHeader.append(s);
                summaryHeader.append("_%_above_");
                summaryHeader.append(thresh);
            }
        }

        summaryOut.printf("%s%n",summaryHeader);
    }

    private void printGeneStats(List<Tuple<String,DepthOfCoverageStats>> statsByGene) {
        logger.debug("statsByTarget size is "+Integer.toString(statsByTarget.size()));
        logger.debug("Initializing refseq...");
        LocationAwareSeekableRODIterator refseqIterator = initializeRefSeq(); //TODO big time this must be a feature context
        logger.debug("Refseq init done.");
        List<Pair<String,DepthOfCoverageStats>> statsByGene = new ArrayList<Pair<String,DepthOfCoverageStats>>();// maintains order
        Map<String,DepthOfCoverageStats> geneNamesToStats = new HashMap<String,DepthOfCoverageStats>(); // allows indirect updating of objects in list


        //TODO THIS is the merging code that needs to change, this whole approach looks inadequate, this will probably have to ha[pen over time

        PrintStream geneSummaryOut = getCorrectStream(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
        StringBuilder summaryHeader = new StringBuilder();
        summaryHeader.append("Gene");
        summaryHeader.append(separator);
        summaryHeader.append("total_coverage");
        summaryHeader.append(separator);
        summaryHeader.append("average_coverage");

        if ( !statsByGene.isEmpty() ) {
            // Only need to get the first item in statsByGene since all have the same samples
            for (String s : statsByGene.get(0).second.getAllSamples()) {
                summaryHeader.append(separator);
                summaryHeader.append(s);
                summaryHeader.append("_total_cvg");
                summaryHeader.append(separator);
                summaryHeader.append(s);
                summaryHeader.append("_mean_cvg");
                summaryHeader.append(separator);
                summaryHeader.append(s);
                summaryHeader.append("_granular_Q1");
                summaryHeader.append(separator);
                summaryHeader.append(s);
                summaryHeader.append("_granular_median");
                summaryHeader.append(separator);
                summaryHeader.append(s);
                summaryHeader.append("_granular_Q3");
                for (int thresh : coverageThresholds) {
                    summaryHeader.append(separator);
                    summaryHeader.append(s);
                    summaryHeader.append("_%_above_");
                    summaryHeader.append(thresh);
                }
            }
        }

        geneSummaryOut.printf("%s%n",summaryHeader);

        for ( Pair<String,DepthOfCoverageStats> geneStats : statsByGene ) {
            printTargetSummary(geneSummaryOut,geneStats);
        }
    }

    private void printTargetSummary(PrintStream output, Pair<?,DepthOfCoverageStats> intervalStats) {
        DepthOfCoverageStats stats = intervalStats.second;
        int[] bins = stats.getEndpoints();

        StringBuilder targetSummary = new StringBuilder();
        targetSummary.append(intervalStats.first.toString());
        targetSummary.append(separator);
        targetSummary.append(stats.getTotalCoverage());
        targetSummary.append(separator);
        targetSummary.append(String.format("%.2f",stats.getTotalMeanCoverage()));

        for ( String s : stats.getAllSamples() ) {
            targetSummary.append(separator);
            targetSummary.append(stats.getTotals().get(s));
            targetSummary.append(separator);
            targetSummary.append(String.format("%.2f", stats.getMeans().get(s)));
            targetSummary.append(separator);
            int median = getQuantile(stats.getHistograms().get(s),0.5);
            int q1 = getQuantile(stats.getHistograms().get(s),0.25);
            int q3 = getQuantile(stats.getHistograms().get(s),0.75);
            targetSummary.append(formatBin(bins,q1));
            targetSummary.append(separator);
            targetSummary.append(formatBin(bins,median));
            targetSummary.append(separator);
            targetSummary.append(formatBin(bins,q3));
            for ( int thresh : coverageThresholds ) {
                targetSummary.append(String.format("%s%.1f",separator,getPctBasesAbove(stats.getHistograms().get(s),stats.value2bin(thresh))));
            }

        }

        output.printf("%s%n", targetSummary);
    }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
     * @updateTargetTable
     * The idea is to have counts for how many *targets* have at least K samples with
     * median coverage of at least X.
     * To that end:
     * Iterate over the samples the DOCS object, determine how many there are with
     * median coverage > leftEnds[0]; how many with median coverage > leftEnds[1]
     * and so on. Then this target has at least N, N-1, N-2, ... 1, 0 samples covered
     * to leftEnds[0] and at least M,M-1,M-2,...1,0 samples covered to leftEnds[1]
     * and so on.
     */
    private void updateTargetTable(int[][] table, DepthOfCoverageStats stats) {
        int[] cutoffs = stats.getEndpoints();
        int[] countsOfMediansAboveCutoffs = new int[cutoffs.length+1]; // 0 bin to catch everything
        for ( int i = 0; i < countsOfMediansAboveCutoffs.length; i ++) {
            countsOfMediansAboveCutoffs[i]=0;
        }

        for ( String s : stats.getAllSamples() ) {
            int medianBin = getQuantile(stats.getHistograms().get(s),0.5);
            for ( int i = 0; i <= medianBin; i ++) {
                countsOfMediansAboveCutoffs[i]++;
            }
        }

        for ( int medianBin = 0; medianBin < countsOfMediansAboveCutoffs.length; medianBin++) {
            for ( ; countsOfMediansAboveCutoffs[medianBin] > 0; countsOfMediansAboveCutoffs[medianBin]-- ) {
                table[countsOfMediansAboveCutoffs[medianBin]-1][medianBin]++;
                // the -1 is due to counts being 1-based and offsets being 0-based
            }
        }
    }

    private int getQuantile(long[] histogram, double prop) {
        int total = 0;

        for ( int i = 0; i < histogram.length; i ++ ) {
            total += histogram[i];
        }

        int counts = 0;
        int bin = -1;
        while ( counts < prop*total ) {
            counts += histogram[bin+1];
            bin++;
        }

        return bin == -1 ? 0 : bin;
    }

    private String baseCounts(int[] counts) {
        if ( counts == null ) {
            counts = new int[6];
        }
        StringBuilder s = new StringBuilder();
        int nbases = 0;
        for ( byte b : BaseUtils.BASES_EXTENDED ) {
            nbases++;
            if ( includeDeletions || b != BaseUtils.Base.D.base ) {
                s.append((char)b);
                s.append(":");
                s.append(counts[BaseUtils.extendedBaseToBaseIndex(b)]);
                if ( nbases < 6 ) {
                    s.append(" ");
                }
            }
        }

        return s.toString();
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    @Override
    public void close() {
        try {
            for (OutputStream stream : outputs.values()) {
                stream.close();
            }
        } catch (IOException e) {
            throw new GATKException("Error closing output files", e);
        }
    }
}
