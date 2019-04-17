package org.broadinstitute.hellbender.tools.walkers.coverage;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.refseq.RefSeqFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.*;
import java.nio.file.Files;
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
 *                                        as input a {@link org.broadinstitute.hellbender.tools.walkers.coverage.DepthOfCoverage.DepthOfCoveragePartitionedDataStore} object corresponding to the coverage
 *                                        information over the whole interval. This is used to write both "_interval_summary"
 *                                        and "_interval_statistics" files for each partition.
 *
 *  2b. writePerGeneDepthInformation() - Similarly to 2a, this should be called once per gene in the interval traversal and
 *                                       it takes the same gen-interval coverage summary as 2a. TODO Note that this may not always cover the gene at hand because of reasons
 *
 *  3. writeTraversalCumulativeCoverage() - This method takes a {@link org.broadinstitute.hellbender.tools.walkers.coverage.DepthOfCoverage.DepthOfCoveragePartitionedDataStore} object that should correspond
 *                                         to the partitioned counts for every base traversed by DepthOfCoverage aggregated.
 *                                         This method is responsible for outputting the "_cumulative_coverage_counts",
 *                                         "_cumulative_coverage_proportions", "_statistics", and "_summary" file writing.
 *
 */
//TODO this will need to be converted into a class holding csvwriteres
public class CoverageOutputWriter implements Closeable {

    private final EnumSet<DoCOutputType.Partition> partitions;
    private final boolean includeGeneOutput;
    private final boolean omitIntervals;
    private final String separator;
    private boolean printBaseCounts;
    private Map<DoCOutputType,PrintStream> outputs;
    private DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat;
    private boolean omitDepthOutput;
    private List<Integer> coverageThresholds;

    //TODO this corresponds to truncating to 2 decimal places, this is legacy string format behavior from the previous incarnation of DoC, this should probably round instead...
    final static DecimalFormat DOUBLE_FORMAT = new DecimalFormat(".00");

    public enum DEPTH_OF_COVERAGE_OUTPUT_FORMAT {
        TABLE,
        RTABLE,
        CSV
    }

    public CoverageOutputWriter(final DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat,
                                final EnumSet<DoCOutputType.Partition> partitions,
                                final String outputBaseName,
                                final boolean includeGeneOutput,
                                final boolean printBaseCounts,
                                final boolean omitDepthOutput,
                                final boolean omitIntervals,
                                final boolean omitSampleSummary,
                                final boolean omitLocusTable)  throws IOException{
        this.outputFormat = outputFormat;
        this.partitions = partitions;
        this.includeGeneOutput = includeGeneOutput;
        this.omitIntervals = omitIntervals;
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

        if(includeGeneOutput && partitions.contains(DoCOutputType.Partition.sample)) {
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
    private static PrintStream getOutputStream(String outputBaseName, DoCOutputType depthSummaryByLocus) throws IOException {
        return new PrintStream(new BufferedOutputStream(Files.newOutputStream(IOUtils.getPath(depthSummaryByLocus.getFilePath(outputBaseName)))));
    }

    //TODO wildly clean every part of this up incredibly
    public void writeCoverageOutputHeaders(final Map<DoCOutputType.Partition, List<String>> allSamplesPartitionedByStart, final LinkedHashSet<String> allSamples, final List<Integer> coverageThresholds) {
        this.coverageThresholds = coverageThresholds;

        if ( ! omitDepthOutput ) { // print header
            writeDepthOutputSummaryHeader(allSamples);
        }

        // write
        if ( includeGeneOutput && partitions.contains(DoCOutputType.Partition.sample) ) {
            PrintStream geneSummaryOut = getCorrectStream(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
            printLocusSummaryHeader("Gene", geneSummaryOut, allSamplesPartitionedByStart.get(DoCOutputType.Partition.sample));
        }

        //TODO this is a horrible mess, allsamples is wrong here...
        if(!omitIntervals) {
            for (DoCOutputType.Partition partition : partitions) {
                PrintStream intervalSummaryOut = getCorrectStream(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.summary);
                printLocusSummaryHeader("Target", intervalSummaryOut, allSamplesPartitionedByStart.get(partition));
            }
        }
    }

    private void writeDepthOutputSummaryHeader(LinkedHashSet<String> allSamples) {
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

    public void writePerLocusDepthSummary(SimpleInterval locus, Map<DoCOutputType.Partition, Map<String,int[]>> countsBySampleByType, Map<DoCOutputType.Partition,List<String>> identifiersByType, boolean includeDeletions) {
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
                    perSampleOutput.append(getBaseCountsString(countsByID != null ? countsByID.get(s) : null , includeDeletions));
                }
                if ( ! depthCounted ) {
                    tDepth += dp;
                }
            }
            depthCounted = true; // only sum the total depth once
        }

        // remember -- genome locus was printed in map()
        stream.print(locus.getContig()+":"+locus.getStart());
        stream.print(separator+tDepth);
        for (DoCOutputType.Partition type : partitions ) { //Note that this is a deterministic traversal since the underlying set is an EnumSet
            stream.print(separator + DOUBLE_FORMAT.format((double) tDepth / identifiersByType.get(type).size() ) );
        }
        stream.printf("%s%n",perSampleOutput);
    }


    public void writePerIntervalDepthInformation(final DoCOutputType.Partition partition, final SimpleInterval locus, final DepthOfCoverageStats intervalStats) {
        PrintStream output = getCorrectStream(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.summary);
        printTargetSummary(output, locus.toString(), intervalStats);
    }
    public void writePerGeneDepthInformation(final RefSeqFeature gene, final DepthOfCoverageStats intervalStats) {
        PrintStream output = getCorrectStream(DoCOutputType.Partition.sample, DoCOutputType.Aggregation.gene, DoCOutputType.FileType.summary);
        printTargetSummary(output, gene.getGeneName(), intervalStats);
    }
    public void writePerTraversalLocusFiles(final DepthOfCoveragePartitionedDataStore coverageProfilesForEntireTraversal, final DoCOutputType.Partition partition) {
        outputPerLocusCumulativeSummaryAndStatistics(getCorrectStream(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_counts),
                getCorrectStream(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.coverage_proportions),
                coverageProfilesForEntireTraversal.getCoverageByAggregationType(partition),partition);
    }
    public void writeCumulativeOutputSummaryFiles(final DepthOfCoveragePartitionedDataStore coverageProfilesForEntireTraversal, final DoCOutputType.Partition partition) {
        outputPerSampleCumulativeStatisticsForPartition(getCorrectStream(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.statistics), coverageProfilesForEntireTraversal.getCoverageByAggregationType(partition));
        outputCumulativeSummaryForPartition(getCorrectStream(partition, DoCOutputType.Aggregation.cumulative, DoCOutputType.FileType.summary), coverageProfilesForEntireTraversal.getCoverageByAggregationType(partition));
    }
    public void writeOutputIntervalStatistics(final DoCOutputType.Partition partition, final int[][] nTargetsByAvgCvgBySample, final int[] binEndpoints) {
        PrintStream output = getCorrectStream(partition, DoCOutputType.Aggregation.interval, DoCOutputType.FileType.statistics);
        printIntervalTable(output, nTargetsByAvgCvgBySample, binEndpoints);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outward facing writer methods
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //TODO splitme out into two method calls reusing code for readability sake
    private void outputPerLocusCumulativeSummaryAndStatistics(PrintStream output, PrintStream coverageOut, DepthOfCoverageStats stats, DoCOutputType.Partition partitionType) {
        int[] endpoints = stats.getEndpoints();
        int samples = stats.getHistograms().size();

        long[][] baseCoverageCumDist = stats.getLocusCounts();

        // rows - # of samples
        // columns - depth of coverage

        boolean printSampleColumnHeader = outputFormat == DEPTH_OF_COVERAGE_OUTPUT_FORMAT.CSV || outputFormat == DEPTH_OF_COVERAGE_OUTPUT_FORMAT.TABLE;

        StringBuilder header = new StringBuilder();
        if ( printSampleColumnHeader ) {
            // mhanna 22 Aug 2010 - Deliberately force this header replacement to make sure integration tests pass.
            // TODO: Update integration tests and get rid of this.
            header.append(partitionType == DoCOutputType.Partition.readgroup ? "read_group" : partitionType.toString());
        }
        header.append(String.format("%sgte_0",separator));
        for ( int d : endpoints ) {
            header.append(String.format("%sgte_%d",separator,d));
        }
        header.append(String.format("%n"));

        output.print(header);
        coverageOut.print(header);

        for ( int row = 0; row < samples; row ++ ) {
            output.printf("%s_%d","NSamples",row+1);
            for ( int depthBin = 0; depthBin < baseCoverageCumDist[0].length; depthBin ++ ) {
                output.printf("%s%d",separator,baseCoverageCumDist[row][depthBin]);
            }
            output.printf("%n");
        }

        for ( String sample : stats.getAllSamples() ) {
            coverageOut.printf("%s",sample);
            double[] coverageDistribution = stats.getCoverageProportions(sample);
            for ( int bin = 0; bin < coverageDistribution.length; bin ++ ) {
                coverageOut.printf("%s%.2f",separator,coverageDistribution[bin]);
            }
            coverageOut.printf("%n");
        }
    }




    private void printLocusSummaryHeader(final String title, final PrintStream summaryOut, final Collection<String> allSamples) {
        StringBuilder summaryHeader = new StringBuilder();
        summaryHeader.append(title);
        summaryHeader.append(separator);
        summaryHeader.append("total_coverage");
        summaryHeader.append(separator);
        summaryHeader.append("average_coverage");

        //TODO this was probably where the order in the header went wrong compared to gatk3
        for (String s : allSamples) {
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
            for (int thresh : coverageThresholds) {
                summaryHeader.append(separator);
                summaryHeader.append(s);
                summaryHeader.append("_%_above_");
                summaryHeader.append(thresh);
            }
        }

        summaryOut.printf("%s%n", summaryHeader);
    }

    private void printTargetSummary(final PrintStream output, final String locusName, final DepthOfCoverageStats stats) {
        int[] bins = stats.getEndpoints();

        StringBuilder targetSummary = new StringBuilder();
        targetSummary.append(locusName);
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
            int median = CoverageUtils.getQuantile(stats.getHistograms().get(s),0.5);
            int q1 = CoverageUtils.getQuantile(stats.getHistograms().get(s),0.25);
            int q3 = CoverageUtils.getQuantile(stats.getHistograms().get(s),0.75);
            targetSummary.append(formatBin(bins, q1));
            targetSummary.append(separator);
            targetSummary.append(formatBin(bins, median));
            targetSummary.append(separator);
            targetSummary.append(formatBin(bins, q3));
            for ( int thresh : coverageThresholds ) {
                targetSummary.append(String.format("%s%.1f", separator, CoverageUtils.getPctBasesAbove(stats.getHistograms().get(s),stats.value2bin(thresh))));
            }

        }

        output.printf("%s%n", targetSummary);
    }

    private void outputPerSampleCumulativeStatisticsForPartition(PrintStream output, DepthOfCoverageStats stats) {
        int[] leftEnds = stats.getEndpoints();

        StringBuilder hBuilder = new StringBuilder();
        if ( ! outputFormat.equals("rTable")) {
            hBuilder.append("Source_of_reads");
        }
        hBuilder.append(separator);
        hBuilder.append(String.format("from_0_to_%d)%s",leftEnds[0],separator));
        for ( int i = 1; i < leftEnds.length; i++ )
            hBuilder.append(String.format("from_%d_to_%d)%s",leftEnds[i-1],leftEnds[i],separator));
        hBuilder.append(String.format("from_%d_to_inf%n",leftEnds[leftEnds.length-1]));
        output.print(hBuilder.toString());
        Map<String,long[]> histograms = stats.getHistograms();

        for ( Map.Entry<String, long[]> p : histograms.entrySet() ) {
            StringBuilder sBuilder = new StringBuilder();
            sBuilder.append(String.format("sample_%s",p.getKey()));
            for ( long count : p.getValue() ) {
                sBuilder.append(String.format("%s%d",separator,count));
            }
            sBuilder.append(String.format("%n"));
            output.print(sBuilder.toString());
        }
    }

    //TODO this method is abjectly terrible
    private void outputCumulativeSummaryForPartition(PrintStream output, DepthOfCoverageStats stats) {
        if ( outputFormat != DEPTH_OF_COVERAGE_OUTPUT_FORMAT.CSV ) {
            output.printf("%s\t%s\t%s\t%s\t%s\t%s","sample_id","total","mean","granular_third_quartile","granular_median","granular_first_quartile");
        } else {
            output.printf("%s,%s,%s,%s,%s,%s","sample_id","total","mean","granular_third_quartile","granular_median","granular_first_quartile");
        }

        for ( int thresh : coverageThresholds ) {
            output.printf("%s%s%d",separator,"%_bases_above_",thresh);
        }

        output.printf("%n");

        Map<String,long[]> histograms = stats.getHistograms();
        Map<String,Double> means = stats.getMeans();
        Map<String,Long> totals = stats.getTotals();
        int[] leftEnds = stats.getEndpoints();

        for ( Map.Entry<String, long[]> p : histograms.entrySet() ) {
            String s = p.getKey();
            long[] histogram = p.getValue();
            int median = CoverageUtils.getQuantile(histogram,0.5);
            int q1 = CoverageUtils.getQuantile(histogram,0.25);
            int q3 = CoverageUtils.getQuantile(histogram,0.75);
            // if any of these are larger than the higest bin, put the median as in the largest bin
            median =  median == histogram.length-1 ? histogram.length-2 : median;
            q1 = q1 == histogram.length-1 ? histogram.length-2 : q1;
            q3 = q3 == histogram.length-1 ? histogram.length-2 : q3;
            if ( ! outputFormat.equals("csv") ) {
                output.printf("%s\t%d\t%.2f\t%d\t%d\t%d",s,totals.get(s),means.get(s),leftEnds[q3],leftEnds[median],leftEnds[q1]);
            } else {
                output.printf("%s,%d,%.2f,%d,%d,%d",s,totals.get(s),means.get(s),leftEnds[q3],leftEnds[median],leftEnds[q1]);
            }

            for ( int thresh : coverageThresholds ) {
                output.printf("%s%.1f", separator, CoverageUtils.getPctBasesAbove(histogram,stats.value2bin(thresh)));
            }

            output.printf("%n");
        }

        if ( ! outputFormat.equals("csv") ) {
            output.printf("%s\t%d\t%.2f\t%s\t%s\t%s%n","Total",stats.getTotalCoverage(),stats.getTotalMeanCoverage(),"N/A","N/A","N/A");
        } else {
            output.printf("%s,%d,%.2f,%s,%s,%s%n","Total",stats.getTotalCoverage(),stats.getTotalMeanCoverage(),"N/A","N/A","N/A");
        }
    }

    private void printIntervalTable(PrintStream output, int[][] intervalTable, int[] cutoffs) {
        String colHeader = outputFormat.equals("rtable") ? "" : "Number_of_sources";
        output.printf(colHeader + separator+"depth>=%d",0);
        for ( int col = 0; col < intervalTable[0].length-1; col ++ ) {
            output.printf(separator+"depth>=%d",cutoffs[col]);
        }

        output.printf(String.format("%n"));
        for ( int row = 0; row < intervalTable.length; row ++ ) {
            output.printf("At_least_%d_samples",row+1);
            for ( int col = 0; col < intervalTable[0].length; col++ ) {
                output.printf(separator+"%d",intervalTable[row][col]);
            }
            output.printf(String.format("%n"));
        }
    }

    private String formatBin(int[] bins, int quartile) {
        if ( quartile >= bins.length ) {
            return String.format(">%d",bins[bins.length-1]);
        } else if ( quartile < 0 ) {
            return String.format("<%d",bins[0]);
        } else {
            return String.format("%d",bins[quartile]);
        }
    }

    // Creates the output string for base counts
    private String getBaseCountsString(int[] counts, boolean includeDeletions) {
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

    @Override
    public void close() {
        try {
            for (OutputStream stream : outputs.values()) {
                stream.close();
            }
        } catch (IOException e) {
            throw new GATKException("Error closing output files:", e);
        }
    }
}
