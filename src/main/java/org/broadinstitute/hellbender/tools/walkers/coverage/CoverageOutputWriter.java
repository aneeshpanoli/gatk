package org.broadinstitute.hellbender.tools.walkers.coverage;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class CoverageOutputWriter {

    public CoverageOutputWriter(final Set<DoCOutputType.Partition> partitions,
                                final File refSeqGeneList,
                                final boolean omitDepthOutput,
                                final boolean omitIntervals,
                                final boolean omitSampleSummary,
                                final boolean omitLocusTable) {
        this.partitions = partitions;
        this.refSeqGeneList = refSeqGeneList;
        this.omitDepthOutput = omitDepthOutput;
        this.omitIntervals = omitIntervals;
        this.omitSampleSummary = omitSampleSummary;
        this.omitLocusTable = omitLocusTable;
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

    private void printGeneStats(List<Pair<GenomeLoc, CoveragePartitioner>> statsByTarget) {
        logger.debug("statsByTarget size is "+Integer.toString(statsByTarget.size()));
        logger.debug("Initializing refseq...");
        LocationAwareSeekableRODIterator refseqIterator = initializeRefSeq();
        logger.debug("Refseq init done.");
        List<Pair<String,DepthOfCoverageStats>> statsByGene = new ArrayList<Pair<String,DepthOfCoverageStats>>();// maintains order
        Map<String,DepthOfCoverageStats> geneNamesToStats = new HashMap<String,DepthOfCoverageStats>(); // allows indirect updating of objects in list

        for ( Pair<GenomeLoc, CoveragePartitioner> targetStats : statsByTarget ) {
            List<String> genes = getGeneNames(targetStats.first,refseqIterator);
            for (String gene : genes) {
                if ( geneNamesToStats.keySet().contains(gene) ) {
                    logger.debug("Merging "+geneNamesToStats.get(gene).toString()+" and "+targetStats.second.getCoverageByAggregationType(DoCOutputType.Partition.sample).toString());
                    geneNamesToStats.get(gene).merge(targetStats.second.getCoverageByAggregationType(DoCOutputType.Partition.sample));
                } else {
                    DepthOfCoverageStats merger = new DepthOfCoverageStats(targetStats.second.getCoverageByAggregationType(DoCOutputType.Partition.sample));
                    geneNamesToStats.put(gene,merger);
                    statsByGene.add(new Pair<String,DepthOfCoverageStats>(gene,merger));
                }
            }
        }

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
}
