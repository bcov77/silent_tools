#!/usr/bin/env python

# a collection of python routines to deal with silent files that don't require pyrosetta

# Add the silent_tools folder to your path, and then do this to import silent tools

#import distutils
#import os
#import sys
#sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
#import silent_tools


#Duplicates get @number appended to their name
def silent_slurp(fname, clear_scores=False):
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
    match_tag = ""

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

        if (line.startswith("SCORE")):
            sp = line.split()
            current_tag = sp[-1]
            match_tag = current_tag

            number = 0
            while ( current_tag in silent_dict ):
                number += 1
                current_tag = match_tag + "@%i"%number

            silent_dict[current_tag] = []
            silent_tags.append(current_tag)

            if ( clear_scores ):
                line = "SCORE: 0 %s"%match_tag

        if ( match_tag == ""):
            continue

        if (line.endswith(match_tag) or line.startswith("NONCANONICAL")):
            silent_dict[current_tag].append(line)
            continue


    return silent_dict, silent_tags, scoreline



def silent_header(scoreline):
    return "SEQUENCE: A\n%s\n"%scoreline.strip()


def silent_body(silent_dict, tags=None):
    if ( tags is None ):
        tags = list(silent_dict.keys())









