import itertools
import numpy

class GraphicsBackend:
    pass
    
def _addOptions(kwdargs, defaults=None):
    if defaults is None: defaults = {}
    options = []
    defaults.update(kwdargs)
    for key, arg in defaults.items():
        if arg is not None and arg != "":
            options.append("""{key}="{arg}" """.format(key=key, arg=arg))
    return "".join(options)

class SVG(GraphicsBackend):
    _filter_id = 0

    def text(self, x, y, text, size=10, anchor="middle", family="Helvetica", **kwdargs):
        defaults = {}
        assert anchor in ["start", "middle", "end"]
        # print(f"render text: {x}, {y}, {text}, {size}, {anchor}, {family}")
        yield """<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" font-family="{family}" text-anchor="{anchor}" {more}>{text}</text>""".format(
            x=x, y=y, size=size, family=family, anchor=anchor, more=_addOptions(kwdargs, defaults), text=text)

    def text_with_background(self, x, y, text, size=10, anchor="middle", text_color="black", bg="white", bg_opacity=0.8, **kwdargs):
        self._filter_id += 1

        text_filter = [
            """<defs>""",
            """    <filter x="0" y="0" width="1" height="1" id="solid{}">""".format(self._filter_id),
            """        <feFlood flood-opacity="{}" flood-color="{}"/>""".format(bg_opacity, bg),
            """        <feComposite in="SourceGraphic"/>""",
            """    </filter>""",
            """</defs>"""]
        for line in text_filter:
            yield line

        # this is a stoopid hack to get the filter to be fully behind the text, without making it blurry
        kwdargs["fill"] = bg
        kwdargs["filter"] = "url(#solid{})".format(self._filter_id)
        yield from self.text(x, y, text, size, anchor, **kwdargs)

        kwdargs["fill"] = text_color
        del kwdargs["filter"]
        yield from self.text(x, y, text, size, anchor, **kwdargs)

    def text_with_background2(self, x, y, text, size=10, anchor="middle", text_color="black", bg="white", bg_opacity=0.8, bg_padding=2, **kwdargs):
        self._filter_id += 1

        # Estimate text dimensions
        char_width = size * 0.6  # Approximate width of a character, adjust if needed
        text_width = char_width * len(text)
        text_height = size  # Height roughly equal to font size
        text_ascent = size * 0.8  # Ascent is roughly 80% of the font size

        # Background rectangle dimensions
        rect_width = text_width + 2 * bg_padding+20
        rect_height = text_height + 2 * bg_padding

        # Adjust x and y based on the anchor
        if anchor == "middle":
            rect_x = x - rect_width / 2
            rect_y = y - rect_height / 2
            text_x = x
            text_y = y + text_ascent / 2
        elif anchor == "start":
            rect_x = x
            rect_y = y - rect_height / 2
            text_x = x + bg_padding
            text_y = y + text_ascent / 2
        elif anchor == "end":
            rect_x = x - rect_width
            rect_y = y - rect_height / 2
            text_x = x - text_width - bg_padding
            text_y = y + text_ascent / 2

        # Draw the background rectangle
        yield from self.rect(rect_x, rect_y, rect_width, rect_height, fill=bg, opacity=bg_opacity)

        # Draw the text over the background
        kwdargs["fill"] = text_color
        yield from self.text(text_x, text_y, text, size, anchor=anchor, **kwdargs)

    def rect(self, x, y, width, height, **kwdargs):
        defaults = {"fill":"white", "stroke":"black"}
        tag = """<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" {more}/>""".format(
            x=x, y=y, w=width, h=height, more=_addOptions(kwdargs, defaults))
        yield tag

    def line(self, x1, y1, x2, y2, **kwdargs):
        defaults = {"stroke":"black"}
        yield """<line x1="{x1:.2f}" x2="{x2:.2f}" y1="{y1:.2f}" y2="{y2:.2f}" {more} />""".format(
            x1=x1, x2=x2,  y1=y1, y2=y2, more=_addOptions(kwdargs, defaults))
        
    def polygon(self, points, **kwdargs):
        defaults = {"fill":"blue", "stroke":"black"}
        points = ["{:.2f},{:.2f}".format(x, y) for x, y in points]
        yield """<polygon points="{}" {more} />""".format(" ".join(points), more=_addOptions(kwdargs, defaults))

    def line_with_arrows(self, x1, y1, x2, y2, n=None, arrows=None, direction="right",
                         color="black", filled=True,
                         arrow_scale=None, arrowKwdArgs=None, **kwdargs):
        
        defaults = {"stroke":color}
        defaults.update(kwdargs)
        
        yield from self.line(x1, y1, x2, y2, **defaults)
        if arrowKwdArgs is None: arrowKwdArgs = {}

        if arrow_scale is None:
            arrow_scale = kwdargs.get("stroke-width", 1)

        if n is not None:
            arrows = numpy.arange(n) / n

        for arrow in arrows:
            x_arrow = x1+float(x2-x1)*arrow
            y_arrow = y1+float(y2-y1)*arrow
            yield from self.arrow(x_arrow, y_arrow, direction, filled=filled,
                color=color, scale=arrow_scale, **arrowKwdArgs)

    def arrow(self, x, y, direction, color="black", scale=1.0, filled=True, **kwdargs):
        more = _addOptions(kwdargs)

        if filled:
            fill = "fill=\"{color}\"".format(color=color)
            close = " z" # closes the path
        else:
            fill = "fill=\"transparent\""
            close = ""

        if direction == "right":
            path = """<path d="M {x0} {y0} L {x1} {y1} L {x2} {y2}{close}" stroke="{color}" """ \
                   """{fill} xcenter="{xcenter}" {more}/>"""
            a = path.format(
                x0=(x-2.5*scale), y0=(y-5*scale), 
                x1=(x+2.5*scale), y1=y, 
                x2=(x-2.5*scale), y2=(y+5*scale),
                close=close,
                fill=fill,
                color=color,
                xcenter=x,
                more=more)
        elif direction == "left":
            path = """<path d="M {x0} {y0} L {x1} {y1} L {x2} {y2}{close}" stroke="{color}" """ \
                   """{fill} xcenter="{xcenter}" {more}/>"""
            a = path.format(
                x0=(x+2.5*scale), y0=(y-5*scale), 
                x1=(x-2.5*scale), y1=y, 
                x2=(x+2.5*scale), y2=(y+5*scale),
                close=close,
                fill=fill,
                color=color,
                xcenter=x,
                more=more)
        yield a
        
    def block_arrow(self, left, top, width, height, arrow_width, direction, **kwdargs):
        coords = {"stroke": kwdargs.pop("stroke", "none"), "fill":kwdargs.pop("fill", "black")}
        coords["more"] = _addOptions(kwdargs)

        if direction == "right":
            path = """<path d="M {x0} {y0} L {x1} {y1} L {x2} {y2} L {x3} {y3} """ \
                   """L {x4} {y4} z" stroke="{stroke}" fill="{fill}" {more}/>"""
            coords["x0"], coords["y0"] = left, top,
            coords["x1"], coords["y1"] = left+width, top
            coords["x2"], coords["y2"] = left+width+arrow_width, top+height/2
            coords["x3"], coords["y3"] = left+width, top+height
            coords["x4"], coords["y4"] = left, top+height
        else:
            path = """<path d="M {x0} {y0} L {x1} {y1} L {x2} {y2} L {x3} {y3} """ \
                   """L {x4} {y4} z" stroke="{stroke}" fill="{fill}" {more}/>"""

            coords["x0"], coords["y0"] = left, top
            coords["x1"], coords["y1"] = left+width, top
            coords["x2"], coords["y2"] = left+width, top+height
            coords["x3"], coords["y3"] = left, top+height
            coords["x4"], coords["y4"] = left-arrow_width, top+height/2

        path = path.format(**coords)
        yield path


    def start_clipped_group(self, x, y, width, height, name):
        yield """<clipPath id="clip_path_{}"><rect x="{}" y="{}" width="{}" height="{}" /></clipPath>""".format(
            name, x, y, width, height)
        yield """<g clip-path="url(#clip_path_{})">""".format(name)
        
    def stop_clipped_group(self):
        yield "</g>"

    # def hla_typing_results(self, x, y, sample_info):
    #     sample = sample_info["Sample"]
    #     locus = sample_info["Locus"]
    #     alleles = sample_info["alleles"]
    #     resolution = sample_info["resolution"]
    #     pdf_width = 900
    #     pdf_x_margin = 50

    #     # Constants for layout
    #     margin_left = 50
    #     margin_top = 30
    #     line_height = 40
    #     rect_height = 30
    #     rect_width = 200
    #     text_margin = 10
    #     max_alleles_per_row = 4
    #     gap_between_alleles = 20  # Adding a gap between alleles

    #     # Calculate dynamic positions
    #     y_pos = margin_top

    #     # Drawing the title
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Sample:</tspan> <tspan font-weight='bold'>{sample}</tspan>", size=20, anchor="start", family="Arial", fill="#213271")
    #     y_pos += line_height
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Locus:</tspan> <tspan font-weight='bold'>{locus}</tspan>", size=16, anchor="start", family="Arial", fill="#E86349")
    #     y_pos += line_height
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Resolution:</tspan> <tspan font-weight='bold'>4th field</tspan>", size=16, anchor="start", family="Arial", fill="black")
            
    #     y_pos += 60
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Full typing result:</tspan>", size=16, anchor="start", family="Arial", fill="black")
    #     y_pos += line_height
        
    #     # Drawing the alleles horizontally
    #     row_count = 0
    #     x_pos = margin_left
        
    #     for idx, allele in enumerate(alleles):
    #         if idx > 0 and idx % max_alleles_per_row == 0:
    #             y_pos += line_height
    #             x_pos = margin_left
    #             row_count += 1

    #         # Calculate text width for centering (assuming fixed-width font for simplicity)
    #         text_width = len(allele) * 8  # Adjust 8 to match your font's character width

    #         # Center text horizontally in the cell
    #         centered_x_pos = x_pos + (rect_width - text_width) / 2

    #         yield from self.rect(x_pos - text_margin, y_pos - rect_height / 2, rect_width, rect_height, fill="#DDFAFB", stroke="none")
    #         yield from self.text(centered_x_pos, y_pos, allele, size=14, anchor="start", family="Arial", fill="black", **{"font-weight":'bold'})
    #         x_pos += rect_width + gap_between_alleles
        
    #     y_pos += line_height * (row_count + 1)
    
    # def hla_typing_results(self, x, y, sample_info):
    #     sample = sample_info["Sample"]
    #     locus = sample_info["Locus"]
    #     alleles = sample_info["alleles"]
    #     resolution = sample_info["resolution"]
    #     pdf_width = 900
    #     pdf_x_margin = 50

    #     margin_left = 50
    #     margin_top = 30
    #     line_height = 40
    #     rect_height = 30
    #     rect_width = 200
    #     text_margin = 10
    #     gap_between_alleles = 20

    #     y_pos = margin_top

    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Sample:</tspan> <tspan font-weight='bold'>{sample}</tspan>", size=20, anchor="start", family="Arial", fill="#213271")
    #     y_pos += line_height
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Locus:</tspan> <tspan font-weight='bold'>{locus}</tspan>", size=16, anchor="start", family="Arial", fill="#E86349")
    #     y_pos += line_height
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Resolution:</tspan> <tspan font-weight='bold'>4th field</tspan>", size=16, anchor="start", family="Arial", fill="black")
    #     y_pos += 60
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Full typing result:</tspan>", size=16, anchor="start", family="Arial", fill="black")
    #     y_pos += line_height

    #     row_count = 0
    #     x_pos = margin_left
        
    #     for idx, allele in enumerate(alleles):
    #         text_width = len(allele) * 8  # Adjust 8 to match your font's character width
    #         cell_width = rect_width + text_margin * 2 + text_width + gap_between_alleles

    #         if x_pos + cell_width > pdf_width - pdf_x_margin:
    #             y_pos += line_height
    #             x_pos = margin_left
    #             row_count += 1

    #         centered_x_pos = x_pos + (rect_width - text_width) / 2

    #         yield from self.rect(x_pos - text_margin, y_pos - rect_height / 2, rect_width, rect_height, fill="#DDFAFB", stroke="none")
    #         yield from self.text(centered_x_pos, y_pos, allele, size=14, anchor="start", family="Arial", fill="black", **{"font-weight": 'bold'})
    #         x_pos += cell_width
        
    #     y_pos += line_height * (row_count + 1)
    # def hla_typing_results(self, x, y, sample_info):
    #     sample = sample_info["Sample"]
    #     locus = sample_info["Locus"]
    #     alleles = sample_info["alleles"]
    #     resolution = sample_info["resolution"]
    #     tag= "Homozygous" if sample_info["tag"]=="hom" else "Heterozygous"
    #     pdf_width = 900
    #     pdf_x_margin = 50

    #     margin_left = 50
    #     margin_top = 30
    #     line_height = 40
    #     rect_height = 30
    #     text_margin = 10
    #     gap_between_alleles = 20
    #     max_columns = 3  # Set the maximum number of columns to 3

    #     y_pos = margin_top


    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Sample:</tspan> <tspan font-weight='bold'>{sample}</tspan>", size=20, anchor="start", family="Arial", fill="#213271")
    #     y_pos += line_height
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Locus:</tspan> <tspan font-weight='bold'>{locus}</tspan>", size=16, anchor="start", family="Arial", fill="#E86349")
    #     y_pos += line_height
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Resolution:</tspan> <tspan font-weight='bold'>4th field</tspan>", size=16, anchor="start", family="Arial", fill="black")
    #     y_pos += 60
    #     yield from self.text(margin_left, y_pos, f"<tspan font-weight='bold' font-style='italic'>Full typing result {tag} :</tspan>", size=16, anchor="start", family="Arial", fill="black")
    #     y_pos += line_height

    #     x_pos = margin_left
    #     column_count = 0
        
    #     for idx, allele in enumerate(alleles):
    #         text_width = len(allele) * 8  # Adjust 8 to match your font's character width
    #         cell_width = text_width + text_margin * 2 + gap_between_alleles

    #         if column_count >= max_columns or x_pos + cell_width > pdf_width - pdf_x_margin:
    #             y_pos += line_height
    #             x_pos = margin_left
    #             column_count = 0

    #         centered_x_pos = x_pos + text_margin

    #         yield from self.text_with_background2(centered_x_pos, y_pos, allele, size=14, anchor="start", family="Arial", fill="black", bg="#DDFAFB", bg_opacity=1, **{"font-weight": 'bold'})
    #         x_pos += cell_width
    #         column_count += 1
        
    #     y_pos += line_height
    def hla_typing_results(self, x, y, sample_info):
        """
        Draw HLA typing results at position (x, y).

        - Display at most 5 alleles per haplotype; if there are more, append a "+N more" line.
        - If sample_info["tag"] == "hom", draw two identical boxes (homozygous).
          Otherwise draw a single box (heterozygous).
        """
        # 1) Extract sample data
        sample        = sample_info["Sample"]
        locus         = sample_info["Locus"]
        alleles       = sample_info["alleles"]
        is_homozygous = (sample_info.get("tag") == "hom")

        # 2) Limit number of displayed alleles
        max_display      = 5
        display_alleles  = alleles[:max_display]
        remaining        = len(alleles) - max_display
        if remaining > 0:
            display_alleles.append(f"+{remaining} more")

        # 3) Styling & layout parameters
        font_family         = "Arial"
        title_font_size     = 20
        subtitle_font_size  = 16
        allele_font_size    = 14
        line_height         = title_font_size * 1.2    # consistent line spacing
        padding             = 6
        box_gap             = 20

        # 4) Starting coordinates (shift down by one line to avoid clipping)
        current_x = x
        current_y = y + line_height

        # 5) Draw the header lines
        yield from self.text(
            current_x, current_y,
            f"<tspan font-weight='bold' font-style='italic'>Sample:</tspan> "
            f"<tspan font-weight='bold'>{sample}</tspan>",
            size=title_font_size, anchor="start",
            family=font_family, fill="#213271"
        )
        current_y += line_height

        yield from self.text(
            current_x, current_y,
            f"<tspan font-weight='bold' font-style='italic'>Locus:</tspan> "
            f"<tspan font-weight='bold'>{locus}</tspan>",
            size=subtitle_font_size, anchor="start",
            family=font_family, fill="#E86349"
        )
        current_y += line_height

        yield from self.text(
            current_x, current_y,
            f"<tspan font-weight='bold' font-style='italic'>Resolution:</tspan> "
            f"<tspan font-weight='bold'>4th field</tspan>",
            size=subtitle_font_size, anchor="start",
            family=font_family, fill="black"
        )
        current_y += line_height

        tag_str = "Homozygous" if is_homozygous else "Heterozygous"
        yield from self.text(
            current_x, current_y,
            f"<tspan font-weight='bold' font-style='italic'>Full typing result ({tag_str}):</tspan>",
            size=subtitle_font_size, anchor="start",
            family=font_family, fill="black"
        )
        current_y += line_height

        # 6) Compute dimensions for the allele box
        char_width_est = allele_font_size * 0.6
        max_chars      = max(len(a) for a in display_alleles)
        box_width      = max_chars * char_width_est + 2 * padding
        box_height     = len(display_alleles) * line_height + 2 * padding

        # 7) Draw one or two boxes
        box_count = 2 if is_homozygous else 1
        for _ in range(box_count):
            # 7a) Background rectangle
            yield from self.rect(
                current_x, current_y,
                box_width, box_height,
                fill="#DDFAFB",
                stroke="none",
                opacity=1.0
            )

            # 7b) Draw each allele line
            text_x = current_x + padding
            text_y = current_y + padding + (allele_font_size * 0.8)  # account for ascent
            for allele in display_alleles:
                yield from self.text(
                    text_x, text_y,
                    allele,
                    size=allele_font_size,
                    anchor="start",
                    family=font_family,
                    fill="black",
                    **{"font-weight": "bold"}
                )
                text_y += line_height

            # Move right for next box
            current_x += box_width + box_gap


    

class Renderer:
    newid = itertools.count()
    
    def __init__(self, backend, x, y, width, height):
        self.backend = backend
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.id = next(Renderer.newid)
    
    def text(self, x, y, *args, **kwdargs):
        yield from self.backend.text(x+self.x, y+self.y, *args, **kwdargs)
        
    def text_with_background(self, x, y, *args, **kwdargs):
        yield from self.backend.text_with_background(x+self.x, y+self.y, *args, **kwdargs)
        
    def rect(self, x, y, *args, **kwdargs):
        yield from self.backend.rect(x+self.x, y+self.y, *args, **kwdargs)
    
    def line(self, x1, y1, x2, y2, *args, **kwdargs):
        yield from self.backend.line(x1+self.x, y1+self.y, x2+self.x, y2+self.y, *args, **kwdargs)
    
    def line_with_arrows(self, x1, y1, x2, y2, *args, **kwdargs):
        yield from self.backend.line_with_arrows(x1+self.x, y1+self.y, x2+self.x, y2+self.y, *args, **kwdargs)

    def arrow(self, x, y, *args, **kwdargs):
        yield from self.backend.arrow(x+self.x, y+self.y, *args, **kwdargs)

    def block_arrow(self, left, top, *args, **kwdargs):
        yield from self.backend.block_arrow(left+self.x, top+self.y, *args, **kwdargs)

    def polygon(self, points, *args, **kwdargs):
        points = [(x+self.x, y+self.y) for x, y in points]
        yield from self.backend.polygon(points, *args, **kwdargs)

    def render(self, element):
        yield "<!-- {} -->".format(element.name)
        yield from self.backend.start_clipped_group(self.x, self.y, self.width, self.height, self.id)
        # yield self.backend.rect(self.x, self.y, self.width, self.height, fill="blue")

        if hasattr(element, "prerenderers"):
            for prerenderer in element.prerenderers:
                for subelement in prerenderer(self, element):
                    yield "   " + subelement

        
        for subelement in element.render(self):
            yield "   " + subelement

        if hasattr(element, "postrenderers"):
            for postrenderer in element.postrenderers:
                for subelement in postrenderer(self, element):
                    yield "   " + subelement

        yield from self.backend.stop_clipped_group()
        
    def subrenderer(self, x=0, y=0, width=None, height=None):
        if width is None: width = self.width
        if height is None: height = self.height

        assert width <= self.width
        assert height <= self.height, "{} {}".format(height, self.height)

        x += self.x
        y += self.y
        
        renderer = Renderer(self.backend, x, y, width, height)

        return renderer
    
    def hla_typing_results(self,x,y,sample_info):
        yield from self.backend.hla_typing_results(x+self.x, y+self.y, sample_info)

