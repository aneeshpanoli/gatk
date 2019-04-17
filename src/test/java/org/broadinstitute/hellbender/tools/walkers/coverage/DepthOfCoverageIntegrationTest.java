package org.broadinstitute.hellbender.tools.walkers.coverage;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import sun.nio.ch.IOUtil;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static org.testng.Assert.*;

public class DepthOfCoverageIntegrationTest extends CommandLineProgramTest {

    private String getDoCExtensionFromFile(File docOutputFile, String basename) {
        String[] split = docOutputFile.getAbsolutePath().split(basename);
        if (split.length == 2) {
            return split[1];
        } else if (split.length == 1) {
            return "";
        }
        Assert.fail("There was a problem loading DoCExtension: "+basename);
        return null;
    }

    private File getExpectedDataDir() {
        return getTestFile( "expected/");
    }

    @Test
    public void testBaseOutputNoFiltering() throws IOException {
        final String expectedBaseName = "depthofcoveragenofiltering";
        final File baseOutputFile = createTempDir("testVCFModeIsConsistentWithPastResults");
        final File output = IOUtils.createTempFileInDirectory( "depthofcoveragenofiltering", ".csv", baseOutputFile);

        String cmd = "-R /Users/emeryj/hellbender/gatk3TestData/Homo_sapiens_assembly18.fasta " +
                "-I /Users/emeryj/hellbender/gatk3TestData/FHS_indexed_subset.bam " +
                "-L /Users/emeryj/hellbender/gatk3TestData/fhs_jhs_30_targts.interval_list " +
                "-mmq 0 -mbq 0 -dels -baseCounts -pt readgroup -pt sample -pt library --output-format CSV -ct 10 -ct 15 -ct 20 -ct 25";
        cmd += " -O "+output.getAbsolutePath();
        runCommandLine(cmd.split(" "));

        File[] actualFiles = baseOutputFile.listFiles();

        compareOutputDirectories(expectedBaseName, output.getName(), actualFiles);
    }


    private void compareOutputDirectories(final String expectedBaseName, final String actualFileBaseName, final File[] actual) throws IOException {
        List<File> expectedFiles = Arrays.stream(Objects.requireNonNull(getExpectedDataDir().listFiles())).filter(f -> f.getName().contains(expectedBaseName)).sorted().collect(Collectors.toList());
        List<File> actualFiles = Arrays.stream(actual).sorted().collect(Collectors.toList());

        // First check that all the expected files correspond to actual files
        List<String> expectedFilenames = expectedFiles.stream().map(f -> getDoCExtensionFromFile(f, expectedBaseName)).sorted().collect(Collectors.toList());
        List<String> actualFilenames = actualFiles.stream().map(f -> getDoCExtensionFromFile(f, actualFileBaseName)).sorted().collect(Collectors.toList());
        Assert.assertEquals(actualFilenames, expectedFilenames);

        // Now assert that the outputs exactly match with the expected outputs
        for (int i = 0; i < actual.length; i++) {
            IntegrationTestSpec.assertEqualTextFiles(actualFiles.get(i), expectedFiles.get(i));
        }
    }

}