#!/usr/bin/env python

# a collection of python routines to deal with silent files that don't require pyrosetta

# Add the silent_tools folder to your path, and then do this to import silent tools

#import distutils
#import os
#import sys
#sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
#import silent_tools


#Quickly parse a silent file without pulling it into memory
def silent_parse(file, clear_scores=False):


    line1 = next(file)
    line2 = next(file)

    assert( line1.startswith("SEQUENCE:" ) )
    assert( line2.startswith("SCORE:" ) )

    scoreline = line2.strip()
    if ( clear_scores ):
        scoreline = "SCORE: score description"




    silent_offsets = {}
    silent_tags = []

    current_tag = ""
    match_tag = ""

    pos = len(line1) + len(line2)
    while (True):
        try:
            line = next(file)
        except:
            break

        line_len = len(line)

        last_pos = pos
        pos += line_len

        if ( line_len == 0 ):
            continue
        if ( line[0] != "S" ):
            continue
        if ( line.startswith("SEQUENCE") ):
            try:
                next(file)
            except:
                pass
            continue
        if ( line.startswith("SCORE") ):
            sp = line.split()

            tag = sp[-1]

            number = 0
            while ( tag in silent_offsets ):
                number += 1
                current_tag = match_tag + "@%i"%number

            silent_offsets[tag] = last_pos
            silent_tags.append(tag)


    return silent_offsets, silent_tags, scoreline



def silent_header(scoreline):
    return "SEQUENCE: A\n%s\n"%scoreline.strip()


def silent_body(silent_dict, tags=None):
    if ( tags is None ):
        tags = list(silent_dict.keys())









