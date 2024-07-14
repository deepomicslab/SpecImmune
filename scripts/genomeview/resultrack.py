import pysam


from genomeview.utilities import match_chrom_format
from genomeview.track import Track



class ResultTrack(Track):
    def __init__(self, res):
        super().__init__([])
        
        self.res = res
        self.name = "test"
            
    # def layout(self, scale):
    #     self.scale = scale
    #     margin_left = 50
    #     margin_top = 30
    #     line_height = 40
    #     rect_height = 30
    #     rect_width = 180
    #     text_margin = 10
    #     max_alleles_per_row = 4
    #     gap_between_alleles = 20  # Adding a gap between alleles
    #     # calculate the number of rows needed and height of the track
    #     pdf_width = 900
    #     pdf_x_margin = 50
    #     self.rows = []
    #     current_row = 1
    #     current_row_height = 0
    #     for i, allele in enumerate(self.res["alleles"]):
    #         if i % max_alleles_per_row == 0:
    #             current_row += 1
    #             current_row_height = 0
    #         current_row_height = max(current_row_height, rect_height)
    #         self.rows.append(current_row)
    #     self.height = current_row * (current_row_height + line_height) + margin_top +100
    #     print(f"Height: {self.height}")
    
    def layout(self, scale):
        self.scale = scale
        margin_left = 50
        margin_top = 30
        line_height = 40
        rect_height = 30
        text_margin = 10
        gap_between_alleles = 20  # Adding a gap between alleles
        pdf_width = 900
        pdf_x_margin = 50
        max_columns = 3  # Set the maximum number of columns to 3
        self.rows = []
        current_row = 0
        y_pos = margin_top
        x_pos = margin_left
        column_count = 0

        for allele in self.res["alleles"]:
            text_width = len(allele) * 8  # Adjust 8 to match your font's character width
            cell_width = text_width + text_margin * 2 + gap_between_alleles

            if column_count >= max_columns or x_pos + cell_width > pdf_width - pdf_x_margin:
                current_row += 1
                x_pos = margin_left
                y_pos += line_height
                column_count = 0

            self.rows.append(current_row)
            x_pos += cell_width
            column_count += 1

        self.height = (current_row + 1) * line_height + margin_top + 200
        print(f"Height: {self.height}")
        
        
    def render(self, renderer):
        # print(f"Rendering track: {self.name}")
        # yield from renderer.text(10, 10, self.name, **{"font-color": "black"})
        yield from renderer.hla_typing_results(10, 10,self.res)