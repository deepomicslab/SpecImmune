import igv_notebook
igv_notebook.init()

import igv_notebook
igv_notebook.init()
igv_browser= igv_notebook.Browser(
    {
        "genome": "hg19",
        "locus": "chr22:24,376,166-24,376,456",
        "tracks": [{
            "name": "BAM",
            "url": "https://s3.amazonaws.com/igv.org.demo/gstt1_sample.bam",
            "indexURL": "https://s3.amazonaws.com/igv.org.demo/gstt1_sample.bam.bai",
            "format": "bam",
            "type": "alignment"
        }],
        "roi": [
            {
                "name": "ROI set 1",
                "url": "https://s3.amazonaws.com/igv.org.test/data/roi/roi_bed_1.bed",
                "indexed": False,
                "color": "rgba(94,255,1,0.25)"
            }
        ]
    }
)