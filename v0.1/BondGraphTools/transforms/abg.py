import logging

logger = logging.getLogger(__name__)

def load(filename):
    """
    ABG file structure in .xfig format

    See Also:
        http://mcj.sourceforge.net/fig-format.html

    Header is given by:
	string	orientation		("Landscape" or "Portrait")
	string	justification		("Center" or "Flush Left")
	string	units			("Metric" or "Inches")
	string	papersize		("Letter", "Legal", "Ledger", "Tabloid",
					 "A", "B", "C", "D", "E",
					 "A4",   "A3", "A2", "A1", "A0" and "B5")
	float	magnification		(export and print magnification, %)
	string	multiple-page		("Single" or "Multiple" pages)
	int	transparent color	(color number for transparent color for GIF
					 export. -3=background, -2=None, -1=Default,
					 0-31 for standard colors or 32- for user colors)
	# optional comment		(An optional set of comments may be here,
					 which are associated with the whole figure)
	int	resolution coord_system	(Fig units/inch and coordinate system:
					   1: origin at lower left corner (NOT USED)
					   2: upper left)
    Args:
        filename:
    Returns:
    """
    print_options = {"orientation": str,
                     "justification": str,
                     "units": str,
                     "papersize": str,
                     "magnification": float,
                     "multiple-page": str,
                     "transparent color": int}

    with open(filename) as abg:
        line = abg.readline()
        header_comment = []
        raw_bonds = []
        raw_nodes = []
        while line[0] == '#':
            header_comment.append(line)
            line = abg.readline()

        for key, cls in print_options.items():
            print_options[key] = cls(line)
            line = abg.readline()

        figure_comment = []
        while line[0] == "#":
            figure_comment.append(line)
            line = abg.readline()

        while line.isspace():
            line = abg.readline()

        resolution, coords = (int(x) for x in line.split())

        while line.isspace():
            line = abg.readline()
        while line:
            if line.isspace():
                pass
            elif line[0] == '2':
                # polyline, we need to know where all the nodes are before we
                # build the bonds
                raw_bonds.append(line)
            elif line[0] == '4':
                # got ourselves a bond
                raw_nodes.append(_build_node_recipe(line))

def _build_node_recipe(string):
    return {}

class MalformedABG(Exception):
    pass
