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

SILENT_INDEX_VERSION = "3"

# Returns the silent index which allows rapid
#  parsing of the silent file
def get_silent_index(file):

    index_name = get_index_name(file)

    if ( not os.path.exists( index_name ) ):
        return build_silent_index(file)

    if ( os.path.getmtime(file) > os.path.getmtime(index_name) ):
        eprint("Silent file newer than index. Rebuilding index!")
        return build_silent_index(file)

    with open(index_name) as f:
        silent_index = json.loads(f.read())

    if ( validate_silent_index(file, silent_index) ):
        return silent_index

    eprint("Silent file changed size. Rebuilding index!")
    return build_silent_index(file)


def get_silent_structures(file, silent_index, tags):
    with open(file) as f:
        return get_silent_structures_file_open(f, silent_index, tags)

def get_silent_structure(file, silent_index, tag):
    with open(file) as f:
        return get_silent_structure_file_open(f, silent_index, tag)

def get_silent_structures_file_open( f, silent_index, tags ):
    structures = []
    for tag in tags:
        structures.append(get_silent_structure_file_open(f, silent_index, tag))

    return structures


def get_silent_structure_file_open( f, silent_index, tag ):
    assert( tag in silent_index['index'] )
    entry = silent_index['index'][tag]

    f.seek( entry['seek'] )

    first_line = next(f)
    assert(first_line.startswith("SCORE"))

    structure = [first_line]

    while (True):
        try:
            line = next(f)
        except:
            break

        if ( len(line) == 0 ):
            continue
        if ( line[0] == "S" ):  # score or sequence, either way we're done
            break

        structure.append(line)

    return structure



def write_silent_file( file, silent_index, structures ):
    with open(file, "w") as f:
        f.write(silent_header(silent_index))

        for structure in structures:
            f.write("".join(structure))



def cmd(command, wait=True):
    # print ""
    # print command
    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    return str(the_stuff[0]) + str(the_stuff[1])

    
def cmd2(command, wait=True):
    # print ""
    # print command
    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    return str(the_stuff[0]),  str(the_stuff[1]), the_command.returncode

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def get_index_name(file):
    return file + ".idx"


def build_silent_index(file):
    assert(os.path.exists(file))

    with open(file) as f:
        line1 = next(f)
        if ( line1.startswith("SEQUENCE:" ) ):
            line1 = next(f)

        assert( line1.startswith("SCORE:" ) )

        scoreline = line1

    # I'm sorry. If you put description in the name of your pose, it will disappear
    lines = cmd("command grep --byte-offset SCORE: %s | grep -v description | awk '{print $1,$NF}'"%file).strip().split("\n")

    index = defaultdict(lambda : {}, {})
    order = []
    orig_order = []
    unique_tags = True

    for line in lines:
        # print(line)
        sp = line.strip().split()
        name = sp[1]
        orig_order.append(name)
        if ( name in index ):
            number = 1
            while (name + "_%i"%number in index):
                number += 1
            new_name = name + "_%i"%number
            index[new_name]["orig"] = name
            name = new_name
            unique_tags = False

        index[name]["seek"] = int(sp[0][:-7])
        order.append(name)

    size = file_size(file)

    silent_index = {"index":index, "tags":order, "orig_tags":orig_order, "scoreline":scoreline, "size":size, 
                    "unique_tags":unique_tags, "version":SILENT_INDEX_VERSION}

    try:
        f = open(get_index_name(file), "w")
        f.write(json.dumps(silent_index))
        f.close()
    except:
        eprint("Warning!!! Unable to save index file. Must reindex every time!")

    return silent_index

def validate_silent_index(file, silent_index):
    if ( "version" not in silent_index ):
        return False
    if ( silent_index['version'] != SILENT_INDEX_VERSION ):
        eprint("Silentindex from older version of silent_tools")
        return False
    size = file_size(file)
    return size == silent_index["size"]

def file_size(file):
    return int(cmd("du -b %s | awk '{print $1}'"%file).strip())

def silent_header(silent_index):
    return "SEQUENCE: A\n%s\n"%silent_index['scoreline'].strip()










