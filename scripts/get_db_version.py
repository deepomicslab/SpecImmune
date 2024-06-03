import re
import os

def get_IMGT_version(args):
    g_file = "%s/whole/release_version.txt"%(args['db'])

    ## if the release_version.txt does not exist, return "N/A"
    if not os.path.exists(g_file):
        return "N/A"

    version_info = "N/A"
    for line in open(g_file):
        if re.search("# version:", line):
            version_info = line.strip()
    return version_info