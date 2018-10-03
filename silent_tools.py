#!/usr/bin/env python
from __future__ import print_function

# a collection of python routines to deal with silent files that don't require pyrosetta

# Add the silent_tools folder to your path, and then do this to import silent tools

#import distutils
#import os
#import sys
#sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
#import silent_tools

import os
import sys
import subprocess
import json
from collections import defaultdict


# Returns the silent index which allows rapid
#  parsing of the silent file
def get_silent_index(file):

    index_name = get_index_name(file)

    if ( not os.path.exists( index_name ) ):
        return build_silent_index(file)

    with open(index_name) as f:
        silent_index = json.loads(f.read())

    if ( validate_silent_index(file, silent_index) ):
        return silent_index

    eprint("Warning!! Silent file changed size. Rebuilding index!")
    return build_silent_index(file)


def cmd(command, wait=True):
    # print ""
    # print command
    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    return str(the_stuff[0]) + str(the_stuff[1])

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def get_index_name(file):
    return file + ".idx"


def build_silent_index(file):
    assert(os.path.exists(file))

    with open(file) as f:
        line1 = next(f)
        line2 = next(f)

        assert( line1.startswith("SEQUENCE:" ) )
        assert( line2.startswith("SCORE:" ) )

        scoreline = line2.strip()

    # I'm sorry. If you put description in the name of your pose, well..
    lines = cmd("command grep --byte-offset SCORE: %s | grep -v description | awk '{print $1,$NF}'"%file).strip().split("\n")

    index = defaultdict(lambda : {}, {})
    order = []

    for line in lines:
        # print(line)
        sp = line.strip().split()
        index[sp[1]]["seek"] = int(sp[0][:-7])
        order.append(sp[1])

    size = file_size(file)

    silent_index = {"index":index, "tags":order, "scoreline":scoreline, "size":size}

    try:
        f = open(file + ".idx", "w")
        f.write(json.dumps(silent_index))
        f.close()
    except:
        eprint("Warning!!! Unable to save index file. Must reindex every time!")

    return silent_index

def validate_silent_index(file, silent_index):
    size = file_size(file)
    return size == silent_index["size"]

def file_size(file):
    return int(cmd("du -b %s | awk '{print $1}'"%file).strip())

def silent_header(scoreline):
    return "SEQUENCE: A\n%s\n"%scoreline.strip()

def silent_body(silent_dict, tags=None):
    if ( tags is None ):
        tags = list(silent_dict.keys())









