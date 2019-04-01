
package org.broadinstitute.hellbender.tools.walkers.coverage;

/**
 * Models a single output file in the DoC walker.
 *
 * @author mhanna
 * @version 0.1
 */
public class DoCOutputType {
    public enum Partition { readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center }
    public enum Aggregation { locus, interval, gene, cumulative }
    public enum FileType { summary, statistics, coverage_counts, coverage_proportions }

    private final Partition partition;
    private final Aggregation aggregation;
    private final FileType fileType;

    public DoCOutputType(final Partition partition,
                         final Aggregation aggregation,
                         final FileType fileType) {
        this.partition = partition;
        this.aggregation = aggregation;
        this.fileType = fileType;
    }

    public String getFileName(final String baseName) {
        // main output
        if(partition == null)
            return baseName;

        if(baseName.trim().equals("/dev/null"))
            return "/dev/null";

        // mhanna 22 Aug 2010 - Deliberately force this header replacement to make sure integration tests pass.
        // TODO: Update integration tests and get rid of this.
        String partitionType = (partition == DoCOutputType.Partition.readgroup ? "read_group" : partition.toString());

        if(fileType == FileType.coverage_counts || fileType == FileType.coverage_proportions) {
            // coverage counts / proportions files always include aggregation.
            return baseName + "." +
                    partitionType + "_" +
                    aggregation + "_" +
                    fileType;
        }

        return  baseName + "." +
                partitionType + "_" +
                (aggregation == Aggregation.interval || aggregation == Aggregation.gene ? aggregation + "_" : "") +
                fileType;
    }

    public int hashCode() {
        return (partition!=null?partition.ordinal()+1:0) * aggregation.ordinal() * fileType.ordinal();
    }

    public boolean equals(Object other) {
        if(!(other instanceof DoCOutputType))
            return false;
        DoCOutputType otherOutputType = (DoCOutputType)other;
        return partition == otherOutputType.partition &&
                aggregation == otherOutputType.aggregation &&
                fileType == otherOutputType.fileType;
    }
}
