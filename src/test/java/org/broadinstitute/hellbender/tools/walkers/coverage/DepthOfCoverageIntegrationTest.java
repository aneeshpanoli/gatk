package org.broadinstitute.hellbender.tools.walkers.coverage;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

public class DepthOfCoverageIntegrationTest extends CommandLineProgramTest {
    private boolean RUN_TESTS = true;
    private String root = "-T DepthOfCoverage ";

    private String buildRootCmd(String ref, List<String> bams, List<String> intervals) {
        StringBuilder bamBuilder = new StringBuilder();
        do {
            bamBuilder.append(" -I ");
            bamBuilder.append(bams.remove(0));
        } while ( bams.size() > 0 );

        StringBuilder intervalBuilder = new StringBuilder();
        do {
            intervalBuilder.append(" -L ");
            intervalBuilder.append(intervals.remove(0));
        } while ( intervals.size() > 0 );


        return root + "-R "+ref+bamBuilder.toString()+intervalBuilder.toString();
    }


    @Test
    public void testBaseOutputNoFiltering() {
        final String[] intervals = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/fhs_jhs_30_targts.interval_list"};
        final String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/FHS_indexed_subset.bam"};
        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final String expectedDir = "/Users/emeryj/hellbender/gatk3TestData/expectedOutputsTest1";

        String cmd = "-R /Users/emeryj/hellbender/gatk3TestData/Homo_sapiens_assembly18.fasta -I /Users/emeryj/hellbender/gatk3TestData/FHS_indexed_subset.bam -L /Users/emeryj/hellbender/gatk3TestData/fhs_jhs_30_targts.interval_list -mmq 0 -mbq 0 -dels -baseCounts -pt readgroup -pt sample -pt library --outputFormat csv -ct 10 -ct 15 -ct 20 -ct 25";
        cmd += " -O "+output.getAbsolutePath();
        final IntegrationTestSpec spec = new IntegrationTestSpec(cmd,0, null);

        // our base file
        final File baseOutputFile = createTempFile("depthofcoveragenofiltering",".tmp");

        runCommandLine(cmd.split(" "));

        // TODO add tests for the output files?
    }

}