import collections
import logging
import pysam

import os

from genomeview.track import Track
from genomeview.intervaltrack import Interval, IntervalTrack
from genomeview.quickconsensus import MismatchCounts
from genomeview.utilities import match_chrom_format
from genomeview.graphtrack import GraphTrack
import numpy

COLORS = ["blue", "red", "green", "black"]

class Series:
    def __init__(self, x, y, color=None, label=None):
        x, y = zip(*sorted(zip(x,y)))
        self.x = x
        self.y = y
        self.color = color
        self.label = label

def allreads(read):
    return True

def color_by_strand(interval):
    # brightness = 0.2 + (cur_reads[0].mapq/40.0*0.8)

    if interval.strand == "-":
        color = "red"
        if interval.read.is_secondary:
            color = "pink"
    else:
        color = "purple"
        if interval.read.is_secondary:
            color = "blue"
    return color
    
class SingleEndBAMTrack(IntervalTrack):
    """
    Displays bam as single-ended reads
    
    Attributes:
        nuc_colors (dict): defines the SVG colors used to display mismatched nucleotides
        insertion_color, deletion_color, clipping_color (str): SVG colors for insertions, 
             deletions and soft/hard clipping

        quick_consensus (bool): specify whether the quick consensus mode should be used. When activated, 
            mismatches wrt the reference genome are only shown when at least several reads support
            a variant at that position (useful when displaying high-error rate data types eg 
            pacbio). Only relevant if draw_mismatches is also True. (default: True)
        draw_mismatches (bool): whether to show mismatches with respect to the reference genome.
            (default: True).

        include_secondary (bool): whether to draw alignments specified as "secondary" in the BAM flags 
            (default: True).
        
        include_read_fn: callback function used to specify which reads should be included in 
            the display. The function takes as its only argument a read (pysam.AlignedSegment) 
            and returns True (yes, display the read) or False (no, don't display). If this 
            function is not specified, by default all reads are shown.
    """
    def __init__(self, bam_path, name=None, bam_type="mormal"):
        """
        Args:
            bam_path (str): path of the bam file to display
            name (str): name of the track (optional - use None if you don't want to specify a name)

        """
        super().__init__([], name=name, bam_type=bam_type)

        self.bam_path = bam_path
        self.bam = pysam.AlignmentFile(bam_path)
        self.intervals = self
        self.mismatch_counts = None
        self.bam_type = bam_type
        
        self.nuc_colors = {"A":"blue", "C":"orange", "G":"green", "T":"black", "N":"gray"}
        self.insertion_color = "purple"
        self.clipping_color = "cyan"
        self.deletion_color = "cyan"

        self.quick_consensus = True
        self.draw_mismatches = True
        self.include_secondary = True

        self.min_indel_size = 0
        self.min_insertion_label_size = 5
        self.min_cigar_line_width = 2
        
        self.draw_read_labels = False

        self.include_read_fn = allreads
        self.max_window_height = 50
        self.window_height_count=0
        # self.color_fn = color_by_strand
        
    def fetch(self):
        """
        Iterator over reads from the bam file

        Overload this method in subclasses to feed this track reads from a different source
        (for example, reads that are already in memory, rather than being read from a file).
        """
        chrom = self.match_chrom_format(self.scale.chrom)
        start, end = self.scale.start, self.scale.end
        
        for read in self.bam.fetch(chrom, start, end):
            
            if not self.include_read_fn or self.include_read_fn(read):
                # print(read)

                yield read
        
    def __iter__(self):
        c = 0
        for i, read in enumerate(self.fetch()):
            c += 1
            if self.bam_type == "normal":
                if read.is_unmapped: continue
                if read.is_secondary and not self.include_secondary: continue
            id_ = read.query_name + str(i)
            interval = Interval(id_, self.scale.chrom, read.reference_start, read.reference_end, 
                                not read.is_reverse)
            print(interval)
            interval.read = read
            if self.draw_read_labels:
                interval.label = read.query_name
            yield interval

    def match_chrom_format(self, chrom):
        """
        Ensures that the input argument `chrom` matches the chromosome name formatting in
        the bam file being visualized (ie "chr14" vs "14").
        """
        return match_chrom_format(chrom, self.bam.references)
        
    def layout(self, scale):
        super().layout(scale)
        self.reset_mismatch_counts()

    def reset_mismatch_counts(self):
        if self.quick_consensus and self.draw_mismatches:
            print("tallying reads")
            self.mismatch_counts = MismatchCounts(
                self.scale.chrom, self.scale.start, self.scale.end)
            print("done tallying reads")
        

            # workaround for some quirk of pysam with crams and large cigars
            # (or something like that, opening a fresh file handle seems to fix the issue)
            bam = pysam.AlignmentFile(self.bam_path)
            self.mismatch_counts.tally_reads(bam)
            print("done tallying reads")
            # for i in range(len(self.mismatch_counts.counts)):
            #     for j in range(len(self.mismatch_counts.counts[i])):
            #         print(i, j, self.mismatch_counts.counts[i][j])


    def layout_interval(self, interval):
        super().layout_interval(interval)

    def draw_interval(self, renderer, interval):
        """
        Draw a read and then, if ``self.draw_mismatches`` is True, draw mismatches/indels 
        on top.
        """

        yield from super().draw_interval(renderer, interval)

        if self.draw_mismatches:
            yield from self._draw_cigar(renderer, interval)

    def _draw_mismatch(self, renderer, length, genome_position, sequence_position, yoffset, alnseq):
        extras = {"stroke":"none"}
        # print(length, genome_position, sequence_position, yoffset, alnseq)

        for i in range(length):
            if genome_position+i < self.scale.start: continue
            if genome_position+i >= self.scale.end: break

            alt = alnseq[sequence_position+i]
            try:
                ref = self.scale.get_seq(genome_position+i, genome_position+i+1)
                # if i==81:
                #     print(alt, ref)
            except AssertionError:
                logging.warn("Unable to get reference sequence; will not draw mismatches")
                return
            

            if alt != ref:
                curstart = self.scale.topixels(genome_position+i)
                curend = self.scale.topixels(genome_position+i+1)

                color = self.nuc_colors[alnseq[sequence_position+i]]
                # if not self.mismatch_counts or alt=="N" or self.mismatch_counts.query(alt, genome_position+i):
                width = max(curend-curstart, self.min_cigar_line_width)
                midpoint = (curstart+curend)/2
                yield from renderer.rect(midpoint-width/2, yoffset, width, self.row_height, fill=color, 
                                            **extras)

    def _draw_deletion(self, renderer, length, genome_position, yoffset):
        extras = {"stroke":"none"}

        if length > self.min_indel_size:
            curstart = self.scale.topixels(genome_position)
            curend = self.scale.topixels(genome_position+length+1)

            if genome_position > self.scale.end: return
            if genome_position+length < self.scale.start: return
            
            width = max(curend-curstart, self.min_cigar_line_width*2)
            midpoint = (curstart+curend)/2
            ymid = yoffset+self.row_height/2

            yield from renderer.rect(midpoint-width/2, yoffset, width, self.row_height, fill="white", 
                                     **extras)
            yield from renderer.line(midpoint-width/2, ymid, midpoint+width/2, ymid, 
                                     stroke="black", **{"stroke-width":1})

    def _draw_insertion(self, renderer, length, genome_position, yoffset):
        if length > self.min_indel_size:
            curstart = self.scale.topixels(genome_position-0.5)
            curend = self.scale.topixels(genome_position+0.5)

            if genome_position > self.scale.end: return
            if genome_position < self.scale.start: return

            midpoint = (curstart+curend)/2

            stroke_width = 1
            width = stroke_width
            ibeam_extension = 1.0
            
            font_size = self.row_height * 0.8
            draw_label = False
            length_string = str(length)
            label_width = len(length_string) * font_size * 0.9

            if length >= self.min_insertion_label_size:
                draw_label = True
            else:
                if label_width < self.scale.relpixels(length*1.5):
                    draw_label = True
            if draw_label:
                width = label_width

            yield from renderer.line(
               midpoint-width/2-ibeam_extension, yoffset+stroke_width/2, 
               midpoint+width/2+ibeam_extension, yoffset+stroke_width/2, stroke=self.insertion_color, 
               **{"stroke-width":stroke_width})

            yield from renderer.rect(
               midpoint-width/2, yoffset, width, self.row_height, 
               fill=self.insertion_color, **{"stroke":"none"})

            yield from renderer.line(
               midpoint-width/2-ibeam_extension, yoffset+self.row_height-stroke_width/2, 
               midpoint+width/2+ibeam_extension, yoffset+self.row_height-stroke_width/2, 
               stroke=self.insertion_color, **{"stroke-width":stroke_width})

            if draw_label:
                yield from renderer.text(midpoint, yoffset+self.row_height*0.75, length_string,
                    size=font_size, fill="white", **{"font-weight":"bold"})

    def _draw_clipping(self, renderer, length, genome_position, yoffset):
        extras = {"stroke":"none"}

        if length >= 5:
            # always draw clipping, irrespective of consensus sequence or mode
            curstart = self.scale.topixels(genome_position-0.5)
            curend = self.scale.topixels(genome_position+0.5)

            width = max(curend-curstart, self.min_cigar_line_width*2)
            midpoint = (curstart+curend)/2

            yield from renderer.rect(midpoint-width/2, yoffset, width, self.row_height, fill=self.clipping_color,
                                     **extras)

    def _draw_cigar(self, renderer, interval):
        """
        draw mismatches/insertions/deletions and clipping
        """
        read = interval.read
        if read.is_secondary: return
        
        # min_width = 2

        row = self.intervals_to_rows[interval.id]
        yoffset = row*(self.row_height+self.margin_y)+20

        genome_position = read.reference_start
        sequence_position = 0
        alnseq = read.query_sequence

        for code, length in read.cigartuples:
            length = int(length)
            if code == 0: #"M":
                yield from self._draw_mismatch(renderer, length, genome_position, sequence_position, yoffset, alnseq)

                sequence_position += length
                genome_position += length
            elif code in [2,3]: #in "D":
                # if not self.mismatch_counts or self.mismatch_counts.query("DEL", genome_position, genome_position+length+1):
                yield from self._draw_deletion(renderer, length, genome_position, yoffset)

                genome_position += length
            elif code == 1: # I
                # if not self.mismatch_counts or self.mismatch_counts.query("INS", genome_position-2, genome_position+2):
                yield from self._draw_insertion(renderer, length, genome_position, yoffset)

                sequence_position += length
            elif code in [4, 5]: #"HS":
                yield from self._draw_clipping(renderer, length, genome_position, yoffset)

                if code == 4:
                    sequence_position += length

    def render_label(self, renderer):
        if self.name is not None:
            yield from renderer.text_with_background(5, 14, self.name, anchor="start", size=18, bg_opacity=0.9)

    




class PairedEndBAMTrack(SingleEndBAMTrack):
    """
    Displays paired-end reads together (otherwise, same as :py:class:`genomeview.SingleEndBAMTrack`).

    Attributes:
        overlap_color: color used to highlight portions of read pairs that are overlapping one another
    """
    def __init__(self, bam_path, name=None):
        super().__init__(bam_path, name)

        self.overlap_color = "lime"

    def layout(self, scale):
        if scale == self.scale: return
        
        self.scale = scale
        self.reset_mismatch_counts()
        
        # for chrom, start, end in self.scale.regions():
        chrom, start, end = self.scale.chrom, self.scale.start, self.scale.end
        cur_read_coords = collections.defaultdict(list)

        for read in self.fetch():
            if read.is_unmapped: continue
            cur_read_coords[read.query_name].append(
                (read.reference_start, read.reference_end, read.next_reference_start, read.is_proper_pair))
        
        # a bit of hocus-pocus to deal with reads whose mates map outside of our region of interest
        for pair in cur_read_coords.values():
            if len(pair) == 1:
                read_end = pair[0]
                if read_end[3]:
                    pair.append((read_end[2], read_end[2]))
                pair.sort()

        for read_name, coords in sorted(cur_read_coords.items(), key=lambda x: x[1]):
            pair_start = coords[0][0]
            pair_end = coords[-1][1]
            interval = Interval(read_name, chrom, pair_start, pair_end)
            if self.draw_read_labels:
                interval.label = read_name
            self.layout_interval(interval)
            #self.intervals.append(interval)
                
        self.height = (len(self.rows)+1) * (self.row_height+self.margin_y)
    
    
    def draw_read_pair(self, renderer, reads):
        reads = [read for read in reads if not read.is_unmapped]
        if len(reads) == 0: return

        chrom = reads[0].reference_name
        row = self.intervals_to_rows[reads[0].query_name]
        
        pair_start = None
        if len(reads) > 1:
            pair_start = reads[0].reference_end
            pair_end = reads[-1].reference_start
        elif reads[0].is_proper_pair:
            # some more hocus-pocus to deal with reads whose mates map outside of our region of interest
            if reads[0].next_reference_start < reads[0].reference_start:
                pair_start = reads[0].next_reference_start
                pair_end = reads[0].reference_start
            else:
                pair_start = reads[0].reference_start
                pair_end = reads[0].next_reference_start

        if pair_start is not None:
            x1 = self.scale.topixels(pair_start)
            x2 = self.scale.topixels(pair_end)
            y = row*(self.row_height+self.margin_y) + self.row_height/2 # refactor

            yield from renderer.line(x1, y, x2, y, **{"stroke-width":1, "stroke":"gray"})
        
        for i, read_end in enumerate(reads):
            interval = Interval(read_end.query_name, chrom, read_end.reference_start,
                                read_end.reference_end, not read_end.is_reverse)
            
            if self.draw_read_labels:
                interval.label = "{}_{}".format(read_end.query_name, 1 if read_end.is_read1 else 2)

            interval.read = read_end

            yield from self.draw_interval(renderer, interval)

            if i == 1 and self.draw_read_labels:
                end = self.scale.topixels(read_end.reference_end)
                top = row*(self.row_height+self.margin_y)

                yield from renderer.text(end+self.label_distance, top+self.row_height,
                                         read_end.query_name, anchor="start")

            
    def render(self, renderer):
        read_buffer = {}
        for read in self.fetch():

            if self.bam_type == "normal":
                if read.is_unmapped: continue
            if read.query_name in read_buffer:
                other_read = read_buffer.pop(read.query_name)
                cur_reads = [other_read, read]
            else:
                read_buffer[read.query_name] = read
                continue
            
            for read_end in cur_reads:
                yield from self.draw_read_pair(renderer, cur_reads)
                
        for read_name, read in read_buffer.items():
            yield from self.draw_read_pair(renderer, [read])
        
        for x in self.render_label(renderer):
            yield x



def _get_filter_fn(keyfn, value):
    def filter_fn(read):
        return keyfn(read) == value
    return filter_fn

def get_group_by_tag_fn(tag):
    """
    creates a grouping function based on the values of "tag", to be used by GroupedBAMTrack
    for example, use tag="HP" with 10x genomics data to split the view into reads from 
    haplotype 1, haplotype 2, and those missing haplotype information
    """
    def group_by_tag(read):
        if not read.has_tag(tag):
            return "missing"
        return str(read.get_tag(tag))
    return group_by_tag

class GroupedBAMTrack(Track):
    """
    Displays reads from a BAM, separated out into groups based on a feature of the reads. 
    For example, group reads based on the value of tag.

    Attributes:
        keyfn: the function used to specify the groupings of reads. Takes as input a read 
            (:py:class:`pysam.AlignedSegment`).
        bam_track_class: the class used to display each group of reads, should probably be
            either :class:`genomeview.bamtrack.PairedEndBAMTrack` or 
            :class:`genomeview.bamtrack.SingleEndBAMTrack`

        space_between (float): the amount of space (pixels) between groups. (Default: 10)
        category_label_fn: a function that nicely formats the category labels. Takes as argument
            the result of the keyfn and should return a string. (Default: render as string)
    """
    def __init__(self, bam_path, keyfn, bam_track_class, name=None):
        """
        """
        super().__init__(name)
        self.keyfn = keyfn
        self.bam_track_class = bam_track_class
        self.bam_path = bam_path
        self.bam = pysam.AlignmentFile(bam_path)
        self.subtracks = []
        
        self.space_between = 10
        self.category_label_fn = str
        
    def layout(self, scale):
        self.scale = scale
        
        categories = set()
        chrom = match_chrom_format(self.scale.chrom, self.bam.references)
        
        for read in self.bam.fetch(chrom, self.scale.start, self.scale.end):
            category = self.keyfn(read)
            categories.add(category)
        
        categories = sorted(categories)
        self.height = 0
        for category in categories:
            cur_track = self.bam_track_class(self.bam_path, name=self.category_label_fn(category))
            cur_track.include_read_fn = _get_filter_fn(self.keyfn, category)
            cur_track.layout(scale)
            self.height += cur_track.height + self.space_between
            
            self.subtracks.append(cur_track)
            
    def render(self, renderer):
        cury = 0
        for subtrack in self.subtracks:
            subrenderer = renderer.subrenderer(y=cury, height=subtrack.height)
            yield from subrenderer.render(subtrack)
            cury += subtrack.height + self.space_between


class BAMCoverageTrack(GraphTrack):
    def __init__(self, bam_path, name=None):
        if name is None:
            name = os.path.basename(bam_path.split(".")[0])
        super().__init__(name=name)
        
        self.bam_path = bam_path
        self.bam = pysam.AlignmentFile(bam_path)
        height = 50
        
    def layout(self, scale):
        import numpy as np
        import pandas as pd

        super().layout(scale)

        chrom = match_chrom_format(scale.chrom, self.bam.references)
        counts = collections.defaultdict(int)
        
        for read in self.bam.fetch(chrom, scale.start, scale.end):
            for i in read.get_reference_positions():
                counts[i] += 1
        
        x = np.arange(scale.start, scale.end+1)
        y = []
        for i, curx in enumerate(x):
            y.append(counts[curx])
        y = np.array(y)
        
        s = pd.Series(y, index=x).sort_index()
        s = s[(s!=s.shift(-1))|(s!=s.shift(1))]
        
        x = s.index
        y = s.values
        
        if len(x):
            self.add_series(x, y)

    def add_series(self, x, y, color=None, label=None):
        """
        Add a dataset corresponding to a single line in the track (ie, a "series"). Note that
        while a single GraphTrack can visualize multiple datasets, they are all plotted on 
        the same y-axis and so should share the same units.

        Arguments:
            x: a list of genomic coordinates
            y: a list of data values; each y-value must correspond to a single x-value
            color: an SVG color for the line being plotted
            label: an optional text label for the graph being plotted (currently unused)
        """
        if label is None:
            label = "series_{}".format(len(self.series))
            
        assert label not in self.series

        x = numpy.asarray(x)
        y = numpy.asarray(y)

        if color is None:
            color = COLORS[len(self.series) % len(COLORS)]
            
        self.series[label] = Series(x, y, color, label)

        self.min_y = min(self.min_y, y[numpy.isfinite(y)].min())
        self.max_y = max(self.max_y, y[numpy.isfinite(y)].max())

    def render_label(self, renderer):
        if self.name is not None:
            yield from renderer.text_with_background(5, 14, self.name, anchor="start", size=18, bg_opacity=0.9)
            
    def render(self, renderer):
        # for label, series in self.series.items():
        #     for i in range(len(series.x)-1):
        #         if any(numpy.isnan(series.x[i:i+2])) or any(numpy.isnan(series.y[i:i+2])):
        #             continue
        #         x1 = self.scale.topixels(series.x[i])
        #         x2 = self.scale.topixels(series.x[i+1])
        #         y1 = self.ytopixels(series.y[i]) +20
        #         y2 = self.ytopixels(series.y[i+1])+20
                
        #         yield from renderer.line(x1, y1, x2, y2, 
        #             **{"stroke-width":1, "stroke":series.color, "stroke-linecap":"square"})
        # plot histogram
        # plot polygon
        points = []
        for i in range(len(self.series["series_0"].x)):
            x = self.scale.topixels(self.series["series_0"].x[i])
            y = self.ytopixels(self.series["series_0"].y[i])+20
            points.append((x, y))
        points.append((self.scale.topixels(self.series["series_0"].x[-1]), self.height))
        points.append((self.scale.topixels(self.series["series_0"].x[0]), self.height))
        yield from renderer.polygon(points, fill="#E86349", stroke="none")

        # since the labels are drawn at the top of the ticks, let's make sure the top tick/label is 
        # more than 12 pixels from the top of the track so it doesn't get clipped
        # TODO: this ignores the margin, as of now
        axis_max_y = self.min_y + (self.max_y - self.min_y) * (1-7/self.height)

        # ticks = get_ticks(self.min_y, axis_max_y, 4)
        ticks = numpy.linspace(self.min_y, axis_max_y, 4)

        yield from renderer.line(0, self.ytopixels(ticks[0])+20, 1, self.ytopixels(ticks[-1])+20, 
                                 **{"stroke-width":2, "stroke":"gray", "stroke-linecap":"square"})
        for tick in ticks:
            if self.max_y > 1_000:
                label = "{:.1g}".format(tick)
            elif self.max_y < 1:
                label = "{:.1f}".format(tick)
            else:
                label = "{:,.0f}".format(tick)

            y = self.ytopixels(tick)+20
            yield from renderer.line(0, y, 10, y, 
                                     **{"stroke-width":2, "stroke":"gray", "stroke-linecap":"square"})
            yield from renderer.text(14, y, label, anchor="start", fill="gray")
            
        for x in self.render_label(renderer):
            yield x