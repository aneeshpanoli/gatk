package org.broadinstitute.hellbender.tools.walkers.coverage;

import htsjdk.samtools.util.Tuple;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.*;

public class CoverageOutputWriter implements Closeable {

    private final Set<DoCOutputType.Partition> partitions;
    private final Path refSeqGeneList;
    private final String separator;
    private boolean printBaseCounts;
    private Map<DoCOutputType,PrintStream> outputs;
    private DEPTH_OF_COVERAGE_OUTPUT_FORMAT outputFormat;

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


    private OutputStream getCorrectStream(DoCOutputType.Partition partition, DoCOutputType.Aggregation aggregation, DoCOutputType.FileType fileType) {
        DoCOutputType outputType = new DoCOutputType(partition,aggregation,fileType);
        if(!outputs.containsKey(outputType))
            throw new UserException(String.format("Unable to find appropriate stream for partition = %s, aggregation = %s, file type = %s",partition,aggregation,fileType));
        return outputs.get(outputType);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outward facing writer methods
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    final static DecimalFormat doubleFormat = new DecimalFormat("#.2#");

    public void writeDepths(Set<DoCOutputType.Partition> partitionTypes, Map<DoCOutputType.Partition,Map<String,int[]>> countsBySampleByType, Map<DoCOutputType.Partition,List<String>> identifiersByType) {
        PrintStream stream = getCorrectStream(null, DoCOutputType.Aggregation.locus, DoCOutputType.FileType.summary)
        // get the depths per sample and build up the output string while tabulating total and average coverage
        StringBuilder perSampleOutput = new StringBuilder();
        int tDepth = 0;
        boolean depthCounted = false;
        for (DoCOutputType.Partition type : partitionTypes ) {
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
        stream.print(separator+tDepth);
        for (DoCOutputType.Partition type : partitionTypes ) {
            stream.print(separator + doubleFormat.format((double) tDepth / identifiersByType.get(type).size() ) );
        }
        stream.write("%s%n",perSampleOutput);
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

    private DepthOfCoverageStats printIntervalStats(List<Pair<GenomeLoc, CoveragePartitioner>> statsByInterval, PrintStream summaryOut, PrintStream statsOut, DoCOutputType.Partition type) {
        Pair<GenomeLoc, CoveragePartitioner> firstPair = statsByInterval.get(0);
        CoveragePartitioner firstAggregator = firstPair.second;
        DepthOfCoverageStats firstStats = firstAggregator.getCoverageByAggregationType(type);

        StringBuilder summaryHeader = new StringBuilder();
        summaryHeader.append("Target");
        summaryHeader.append(separator);
        summaryHeader.append("total_coverage");
        summaryHeader.append(separator);
        summaryHeader.append("average_coverage");

        for ( String s : firstStats.getAllSamples() ) {
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

        int[][] nTargetsByAvgCvgBySample = new int[firstStats.getHistograms().size()][firstStats.getEndpoints().length+1];

        for ( Pair<GenomeLoc, CoveragePartitioner> targetAggregator : statsByInterval ) {

            Pair<GenomeLoc,DepthOfCoverageStats> targetStats = new Pair<GenomeLoc,DepthOfCoverageStats>(
                    targetAggregator.first, targetAggregator.second.getCoverageByAggregationType(type));
            printTargetSummary(summaryOut,targetStats);
            updateTargetTable(nTargetsByAvgCvgBySample,targetStats.second);
        }

        printIntervalTable(statsOut,nTargetsByAvgCvgBySample,firstStats.getEndpoints());

        return firstStats;
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
