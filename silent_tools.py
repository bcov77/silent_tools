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
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np
import re
import struct


SILENT_INDEX_VERSION = "4"

# Returns the silent index which allows rapid
#  parsing of the silent file
def get_silent_index(file):

    index_name = get_index_name(file)

    if ( not os.path.exists( index_name ) ):
        return build_silent_index(file)

    if ( os.path.getmtime(get_real_file(file)) > os.path.getmtime(index_name) ):
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


def get_real_file(file):
    real_file = cmd("realpath %s"%file).strip()
    if ( not os.path.exists(file) or not os.path.exists(real_file) ):
        eprint("silent_tools: Error file doesn't exist: file")
        assert(False)
    return real_file


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
    lines = cmd("command grep -a --byte-offset '^SCORE:' %s | grep -v description | awk '{print $1,$NF}'"%file).strip().split("\n")

    index = defaultdict(lambda : {}, {})
    order = []
    orig_order = []
    unique_tags = True

    dup_index = {}

    for line in lines:
        # print(line)
        sp = line.strip().split()
        name = sp[1]
        orig_order.append(name)
        if ( name in index ):
            # speedup
            if ( name in dup_index ):
                number = dup_index[name]
            else:
                number = 1
            # /speedup

            while (name + "_%i"%number in index):
                number += 1

            # speedup
            dup_index[name] = number
            # /speedup
            new_name = name + "_%i"%number
            index[new_name]["orig"] = name
            name = new_name
            unique_tags = False

        index[name]["seek"] = int(sp[0][:-7])
        order.append(name)

    size = file_size(file)

    silent_index = {"index":index, "tags":order, "orig_tags":orig_order, "scoreline":scoreline, "size":size, 
                    "unique_tags":unique_tags, "version":SILENT_INDEX_VERSION}

    sequence = "A"
    if ( len(order) > 0 ):
        try:
            structure = get_silent_structure(file, silent_index, order[0])
            sequence = "".join(get_sequence_chunks(structure))
        except:
            eprint("Failed to get sequence. Please tell Brian")

    silent_index['sequence'] = sequence


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
    file = get_real_file(file)
    return int(cmd("du -b %s | awk '{print $1}'"%file).strip())

def silent_header(silent_index):
    return "SEQUENCE: %s\n%s\n"%(silent_index['sequence'], silent_index['scoreline'].strip())



def get_sequence_chunks(structure):

    full_sequence = None
    chain_endings = None

    for line in structure:
        if ( line.startswith("ANNOTATED_SEQUENCE") ):
            tmp = line
            tmp = tmp.strip()
            tmp = tmp.split()[1]
            full_sequence = re.sub(r"\[[^]]*\]", "", tmp)
        if ( line.startswith("CHAIN_ENDINGS") ):
            tmp = line
            tmp = tmp.strip()
            tmp = tmp.split()
            chain_endings = [int(x) for x in tmp[1:len(tmp)-1] ]

    bad = False
    if ( full_sequence is None ):
        eprint("silentsequence: no ANNOTATED_SEQUENCE for tag %s"%tag)
        bad = True
    if ( chain_endings is None ):
        #eprint("silentsequence: no CHAIN_ENDINGS for tag %s"%tag)
        #bad = True
        chain_endings=[]

    if (bad):
        return None

    sequence_chunks = []
    last_end = 0
    for end in chain_endings + [len(full_sequence)]:
        sequence_chunks.append( full_sequence[last_end:end] )
        last_end = end

    return sequence_chunks


########
# Everything below this point is a little sketchy



_atom_record_format = (
    "ATOM  {atomi:5d} {atomn:^4}{idx:^1}{resn:3s} {chain:1}{resi:4d}{insert:1s}   "
    "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}\n"
)
def format_atom(
        atomi=0,
        atomn='ATOM',
        idx=' ',
        resn='RES',
        chain='A',
        resi=0,
        insert=' ',
        x=0,
        y=0,
        z=0,
        occ=1,
        b=0
):
    return _atom_record_format.format(**locals())



name1_to_name3 = {
    "R":"ARG",
    "K":"LYS",
    "N":"ASN",
    "D":"ASP",
    "E":"GLU",
    "Q":"GLN",
    "H":"HIS",
    "P":"PRO",
    "Y":"TYR",
    "W":"TRP",
    "S":"SER",
    "T":"THR",
    "G":"GLY",
    "A":"ALA",
    "M":"MET",
    "C":"CYS",
    "F":"PHE",
    "L":"LEU",
    "V":"VAL",
    "I":"ILE",
}

def write_pdb_atoms(atoms, sequence, atom_names):
    lines = []
    assert(len(atoms) / len(sequence) == len(atom_names))

    for i in range(len(sequence)):
        name3 = name1_to_name3[sequence[i]]

        for iatom, atom in enumerate(atom_names):
            atom_offset = i*len(atom_names)+iatom
            a = atoms[atom_offset]

            lines.append( format_atom( 
                    atomi=(atom_offset)%100000,
                    resn=name3,
                    resi=(i+1)%10000,
                    atomn=atom_names[iatom],
                    x=a[0],
                    y=a[1],
                    z=a[2]
                    ))

    return lines



#########
# Everything below this point is sketchy





silent_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"



def code_from_6bit(_8bit):
    _8bit = ord(_8bit[0])
    if ( ( _8bit >= ord('A')) and (_8bit <= ord('Z')) ): return _8bit - ord('A')
    if ( ( _8bit >= ord('a')) and (_8bit <= ord('z')) ): return _8bit - ord('a') + 26
    if ( ( _8bit >= ord('0')) and (_8bit <= ord('9')) ): return _8bit - ord('0') + 52
    if (   _8bit == ord('+')  ): return 62
    return 63


def decode_32_to_24( i0, i1, i2, i3 ):
    i0 = code_from_6bit( i0 )
    i1 = code_from_6bit( i1 )
    i2 = code_from_6bit( i2 )
    i3 = code_from_6bit( i3 )

    o0 = 0xFF & (i0 | (i1 << 6))
    o1 = 0xFF & ((i1 >> 2) | (i2 << 4))
    o2 = 0xFF & ((i3 << 2) | (i2 >> 4))

    return o0, o1, o2

def decode6bit( jar ):

    ba = bytearray()

    valid_bits = 0
    i = 0
    while ( i < len(jar) ):

        this_str = ["!", "!", "!", "!"]

        j = 0
        while ( i < len(jar) and j < 4 ):
            this_str[j] = jar[i]
            i += 1
            j += 1
            valid_bits += 6

        # print(this_str)
        bytess = decode_32_to_24(*this_str)
        # print(bytess)

        ba.append(bytess[0])
        ba.append(bytess[1])
        ba.append(bytess[2])
    valid_bytes = int( valid_bits / 8 )
    ba = ba[:valid_bytes]
    assert(len(ba) % 4 == 0)
    return ba


def silent_line_to_atoms(line):
    ba = decode6bit( line )

    float_packer = struct.Struct("f"*(len(ba)//4))

    floats = float_packer.unpack(ba)

    assert(len(floats) % 3 == 0)

    return np.array(floats).reshape(-1, 3)


def sketch_get_atoms_by_residue(structure):

    sequence = "".join(get_sequence_chunks(structure))
    if ( sequence is None ):
        return None

    residues = []

    for line in structure:

        # Ok, so we're going to use some really crappy detection here
        sp = line.split()
        if ( len(sp) != 2 ):
            continue

        binary = sp[0]
        # Ensure there are lowercase letters
        if ( binary.upper() == binary ):
            continue

        if ( not binary.startswith("L") ):
            continue
        binary = binary[1:]

        residues.append( silent_line_to_atoms( binary ) )

    print(len(sequence), len(residues))
    assert(len(sequence) == len(residues))
    return residues

def sketch_get_atoms(structure, atom_nums):

    atoms_by_res = sketch_get_atoms_by_residue(structure)

    final = []
    for residue in atoms_by_res:
        final.append(residue[atom_nums])

    final = np.array(final).reshape(-1, 3)

    return final


def sketch_get_cas_protein_struct(structure):

    sequence = "".join(get_sequence_chunks(structure))

    cas = []

    for line in structure:
        line = line.strip()
        if (len(line) == 0):
            continue
        sp = line.split()


        if (len(sp) != 13):
            continue

        try:
            seqpos = int(sp[0])
            if ( not sp[1] in "HEL" ):
                raise Exception()
            x = float(sp[5])
            y = float(sp[6])
            z = float(sp[7])
        except:
            continue
        cas.append([x, y, z])

        assert(seqpos == len(cas))

    assert(len(cas) == len(sequence))

    return np.array(cas)





















