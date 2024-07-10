import pysam


from genomeview.intervaltrack import IntervalTrack, Interval
from genomeview.utilities import match_chrom_format


def color_by_nuc(interval):
    colors = {"A":"blue", "C":"orange", "G":"green", "T":"black", "N":"gray"}
    return colors.get(str(interval.variant.alts[0]).upper(), "gray")

class VCFTrack(IntervalTrack):
    def __init__(self, vcf_path, name=None):
        super().__init__([], name=name)
        
        self.vcf_path = vcf_path
        self.vcf = pysam.VariantFile(vcf_path)
        
        # ordinarily, IntervalTrack looks for intervals in a pre-loaded list of intervals;
        # here, we instead make this instance into an iterator to load variants
        # from the vcf as needed
        # self.intervals = self
        
        self.color_fn = color_by_nuc
        self.row_height = 25
        self.min_variant_pixel_width = 2
        self._interval_id = 0  # Counter for unique interval IDs
        self.var_cnt=0
        self.vcf_ids=set()
        _interval_id = 0
        self.intervals = []
        for rec in self.vcf:
            interval_id = f"variant_{_interval_id}"
            self.vcf_ids.add(interval_id)
            _interval_id += 1
            interval = Interval(interval_id, rec.chrom, rec.start, rec.stop+1, None)
            interval.variant = rec
            self.intervals.append(interval)
    
        
    def __iter__(self):
        chrom = match_chrom_format(self.scale.chrom, self.vcf.header.contigs)
        start, end = self.scale.start, self.scale.end
        var_cnt=0
        
        # for variant in self.vcf.fetch():
        #     self._interval_id += 1
        #     interval_id = f"variant_{self._interval_id}"
        #     print(f"Variant: {interval_id}: {variant}")
        #     interval = genomeview.Interval(interval_id, variant.chrom, variant.start, variant.stop+1, None)
        #     interval.variant = variant
        #     print(f"Interval: {interval.id}")
        #     var_cnt+=1
        #     self.var_cnt+=1
        #     yield interval
        # print(f"Total variants: {var_cnt}")
        for interval in self.intervals:
            yield interval
            
    def layout(self, scale):
        self.scale = scale
        self.intervals_to_rows = {}
        current_row = 0
        
        for interval in self.intervals:
            print(f"Layout interval: {interval}")
            self.intervals_to_rows[interval.id] = current_row
            current_row += 1
        
        print(f"Intervals to rows: {self.intervals_to_rows}")
        
    def draw_interval(self, renderer, interval):
        # overriding this method isn't really necessary - we're just going to make
        # sure that every variant is at least several screen pixels wide, even
        # if we're zoomed pretty far out
        # print(f"Drawing interval: {interval}")
        # print(f"var_cnt: {self.var_cnt}")
        start = self.scale.topixels(interval.start)
        end = self.scale.topixels(interval.end)
        
        if end - start < self.min_variant_pixel_width:
            mid = (end + start) / 2
            start = mid - self.min_variant_pixel_width / 2
            end = mid + self.min_variant_pixel_width / 2
        # print(f"Start: {interval.start}, End: {interval.end}")
        # print(f"Variant: {interval.id}")
        row = self.intervals_to_rows[interval.id]
        top = row * (self.row_height + self.margin_y)
        top=30
        
        # print(f"Variant: {interval.variant.alts[0]}")
        color = self.color_fn(interval)

        print(f"render rect: {start}, {top}, {end - start}, {self.row_height}, fill={color}, **{{'stroke': 'none'}}")

        yield from renderer.rect(start, top, end - start, self.row_height, fill=color, 
                                 **{"stroke": "none"})