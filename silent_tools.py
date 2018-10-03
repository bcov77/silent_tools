#!/usr/bin/env python

# a collection of python routines to deal with silent files that don't require pyrosetta

# Add the silent_tools folder to your path, and then do this to import silent tools

#import distutils
#import os
#import sys
#sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
#import silent_tools


# If dont_overwrite_duplicates true, then duplicate entries are simply listed one
# after another in tags
def silent_slurp(fname, clear_scores=False, dont_overwite_duplicates=False):
    f = open(fname)
    lines = f.readlines()
    f.close()

    assert( len(lines) >= 2 )
    assert( lines[0].startswith("SEQUENCE:" ) )
    assert( lines[1].startswith("SCORE:" ) )

    scoreline = lines[1].strip()
    if ( clear_scores ):
        scoreline = "SCORE: score description"

    silent_dict = {}
    silent_tags = []

    current_tag = ""
    skip_next = False
    for line in lines:
        line = line.strip()
        if (len(line) == 0):
            continue

        if (skip_next):
            skip_next = False
            continue

        if (line.startswith("SEQUENCE")):
            skip_next = True
            continue

        sp = line.split()
        if (line.startswith("SCORE")):
            current_tag = sp[-1]

            silent_tags.append(current_tag)

            if ( not dont_overwite_duplicates or current_tag not in silent_dict ):
                silent_dict[current_tag] = []

            if ( clear_scores ):
                sp = ["SCORE:", "0", current_tag]

        if ( current_tag == ""):
            continue

        if (sp[-1] == current_tag or line.startswith("NONCANONICAL")):
            silent_dict[current_tag].append(line)
            continue


    return silent_dict, silent_tags, scoreline
