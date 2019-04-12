package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import javafx.collections.ObservableList;
import javafx.collections.transformation.SortedList;
import org.apache.commons.collections4.SetUtils;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentContextIteratorBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;

public abstract class LocusWalkerByInterval extends LocusWalker {

    private OverlapDetector<Locatable> intervalsToTrack = null;
    private Set<Locatable> previousIntervals = new HashSet<>();

    /**
     * Implementation of locus-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     *
     * The default implementation iterates over all positions in the reference covered by reads (filtered and transformed)
     * for all samples in the read groups, using the downsampling method provided by {@link #getDownsamplingInfo()}
     * and including deletions only if {@link #includeDeletions()} returns {@code true}.
     */
    @Override
    public void traverse() {
        final SAMFileHeader header = getHeaderForReads();
        // get the samples from the read groups
        final Set<String> samples = header.getReadGroups().stream()
                .map(SAMReadGroupRecord::getSample)
                .collect(Collectors.toSet());
        final CountingReadFilter countedFilter = makeReadFilter();
        // get the filter and transformed iterator
        final Iterator<GATKRead> readIterator = getTransformedReadStream(countedFilter).iterator();

        intervalsToTrack = OverlapDetector.create(getIntervalObjectsToQueryOver());

        final AlignmentContextIteratorBuilder alignmentContextIteratorBuilder = new AlignmentContextIteratorBuilder();
        alignmentContextIteratorBuilder.setDownsamplingInfo(getDownsamplingInfo());
        alignmentContextIteratorBuilder.setEmitEmptyLoci(emitEmptyLoci());
        alignmentContextIteratorBuilder.setIncludeDeletions(includeDeletions());
        alignmentContextIteratorBuilder.setKeepUniqueReadListInLibs(keepUniqueReadListInLibs());
        alignmentContextIteratorBuilder.setIncludeNs(includeNs());

        final Iterator<AlignmentContext> iterator = alignmentContextIteratorBuilder.build(
                readIterator, header, userIntervals, getBestAvailableSequenceDictionary(),
                hasReference());

        // iterate over each alignment, and apply the function
        iterator.forEachRemaining(alignmentContext -> {
                    final SimpleInterval alignmentInterval = new SimpleInterval(alignmentContext);
                    apply(alignmentContext, new ReferenceContext(reference, alignmentInterval), new FeatureContext(features, alignmentInterval));
                    progressMeter.update(alignmentInterval);
                }
        );
        for (Locatable l : previousIntervals) {
            onIntervalEnd(l);
        }
        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Tool specified list of Locatable objects (which have been read into memory) that will have overlaps queried at each locus
     *
     * @return A list of Locatable Objects that will be queried over each genomic location
     */
    public abstract List<Locatable> getIntervalObjectsToQueryOver();

    @Override
    public final void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        Set<Locatable> currentIntervals = intervalsToTrack.getOverlaps(alignmentContext);
        Set<Locatable> passedIntervals = SetUtils.difference(previousIntervals, currentIntervals);
        Set<Locatable> newIntervals = SetUtils.difference(currentIntervals, previousIntervals);
        previousIntervals = currentIntervals;

        // Close out existing intervals that have been passed
        for(Locatable l : passedIntervals) {
            onIntervalEnd(l);
        }
        // Open new intervals that have been reached
        for(Locatable l : newIntervals) {
            onIntervalStart(l);
        }

        apply(alignmentContext, referenceContext, featureContext, currentIntervals);
    }

    /**
     * //TODO these comments need to be updated
     * Process an individual AlignmentContext (with optional contextual information). Must be implemented by tool authors.
     * Note that apply() may be called multiple times for the same reference coordinate in the even that there are multiple
     * intervals overlapping the same region. In this case activeInterval will point to the associated interval context.
     *
     * @param alignmentContext current alignment context
     * @param referenceContext Reference bases spanning the current locus. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current locus
     *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current locus. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext, Set<Locatable> activeInterval);

    public abstract void onIntervalStart(Locatable activeInterval);

    public abstract void onIntervalEnd(Locatable activeInterval);
}
