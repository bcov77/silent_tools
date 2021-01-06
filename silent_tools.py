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
import bz2


SILENT_INDEX_VERSION = "5"

# Returns the silent index which allows rapid
#  parsing of the silent file
def get_silent_index(file, accept_garbage=False):

    index_name = get_index_name(file)

    if ( not os.path.exists( index_name ) ):
        return build_silent_index(file, accept_garbage=accept_garbage)

    if ( os.path.getmtime(get_real_file(file)) > os.path.getmtime(index_name) ):
        eprint("Silent file newer than index. Rebuilding index!")
        return build_silent_index(file)

    try:
        with open(index_name) as f:
            silent_index = json.loads(f.read())
    except:
        eprint("Silent index is corrupt. Rebuilding index!")
        return build_silent_index(file)

    if ( validate_silent_index(file, silent_index) ):
        return silent_index

    eprint("Silent file changed size. Rebuilding index!")
    return build_silent_index(file)


def get_silent_structures(file, silent_index, tags):
    with open(file, errors='ignore') as f:
        return get_silent_structures_file_open(f, silent_index, tags)

def get_silent_structure(file, silent_index, tag):
    with open(file, errors='ignore') as f:
        return get_silent_structure_file_open(f, silent_index, tag)

def get_silent_structures_file_open( f, silent_index, tags ):
    structures = []
    for tag in tags:
        structures.append(get_silent_structure_file_open(f, silent_index, tag))

    return structures


def get_silent_structure_file_open( f, silent_index, tag, return_first_line=False ):
    assert( tag in silent_index['index'] )
    entry = silent_index['index'][tag]

    f.seek( entry['seek'] )

    first_line = next(f)
    structure, first_line = rip_structure_by_lines(f, first_line)
    
    if ( return_first_line ):
        return structure, first_line
    else:
        return structure


# can throw
def rip_structure_by_lines_arbitrary_start(f, first_line, save_structure=True):
    while ( not first_line.startswith("SCORE") ):
        first_line = next(f) # throw

    return rip_structure_by_lines(f, first_line, save_structure=save_structure)

# can throw
def rip_structures_till(f, first_line, till_structure):

    while True:
        while ( not first_line.startswith("SCORE") ):
            first_line = next(f) # throw

        cur_tag = first_line.strip().split()[-1]

        if ( cur_tag == till_structure ):
            break

        _, first_line = rip_structure_by_lines(f, first_line, save_structure=False)


    return rip_structure_by_lines(f, first_line, save_structure=True)



def rip_structure_by_lines(f, first_line, save_structure=True):

    assert(first_line.startswith("SCORE"))

    structure = [first_line] if save_structure else None

    while (True):
        try:
            line = next(f)
        except:
            line = None
            break

        if ( len(line) == 0 ):
            continue
        if ( line[0] == "S" and (line.startswith("SCORE") or line.startswith("SEQUENCE"))):  # score or sequence, either way we're done
            break

        if ( save_structure ):
            structure.append(line)

    first_non_structure_line = line
    return structure, first_non_structure_line


def get_silent_structures_true_slice( f, silent_index, idx_start, idx_stop_py, oneline=False, raw_string=False ):
    assert( idx_start >= 0 and idx_stop_py <= len(silent_index['index']) )

    start_seek = silent_index['index'][silent_index['tags'][idx_start]]['seek']

    if ( idx_stop_py == len(silent_index['tags']) ):
        stop_seek = None
    else:
        stop_seek = silent_index['index'][silent_index['tags'][idx_stop_py]]['seek']

    f.seek( start_seek )

    if ( stop_seek is None ):
        data = f.read()
    else:
        data = f.read(stop_seek - start_seek)

    if ( raw_string ):
        return data

    structures = []
    for idx in range(idx_start, idx_stop_py):
        start = silent_index['index'][silent_index['tags'][idx]]['seek']

        if ( idx + 1 == idx_stop_py ):
            stop = None
        else:
            stop = silent_index['index'][silent_index['tags'][idx+1]]['seek']
            assert( stop - start_seek <= len(data) + 1 )

        if ( stop is None ):
            structure_dat = data[start-start_seek:]
        else:
            structure_dat = data[start-start_seek:stop-start_seek]

        if ( not oneline ):
            structure_dat = [ x + "\n" for x in structure_dat.split("\n") if len(x) > 0 ]

        structures.append(structure_dat)

    return structures


def get_real_file(file):
    real_file, error, code = cmd2("realpath %s"%file)
    if ( code != 0 ):
        real_file = cmd("readlink -f %s"%file)
    real_file = real_file.strip()
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

def detect_silent_type(structure):
    is_binary = False
    is_protein = False
    for line in structure:
        if ( len(line) == 0 ):
            continue
        if ( line[0] in "HEL" ):
            is_binary = True
        if ( len(line) < 6 ):
            continue
        if ( line[5] in "HEL" ):
            is_protein = True

    if ( is_binary and is_protein ):
        eprint("silent_tools: Silent file is both BINARY and PROTEIN? Using UNKNOWN")
        return "UNKNOWN"

    if ( is_binary ):
        return "BINARY"
    if ( is_protein ):
        return "PROTEIN"

    eprint("silent_tools: Can't determine silent type. Using UNKNOWN")
    return "UNKNOWN"


def assert_is_silent_and_get_scoreline(file, return_f=False, accept_garbage=False):
    if ( not os.path.exists(file) ):
        sys.exit("silent_tools: Error! Silent file doesn't exist: " + file)

    try:
        if ( file.endswith(".bz2") ):
            f = bz2.open(file, "rt")
        else:
            f = open(file, errors='ignore')
    except:
        sys.exit("silent_tools: Error! Can't open silent file: " + file)

    try:
        line1 = next(f)
    except:
        sys.exit("silent_tools: Error! Silent file is empty: " + file)

    if ( line1.startswith("SEQUENCE:" ) ):
        try:
            line1 = next(f)
        except:
            sys.exit("silent_tools: Error! Truncated silent file: " + file)
    else:
        eprint("silent_tools: Warning! Silent file doesn't have SEQUENCE line")

    if ( not line1.startswith("SCORE:" ) ):
        if ( accept_garbage ):
            eprint("silent_tools: Error! Silent file doesn't have SCORE: header")
        else:
            sys.exit("silent_tools: Error! Silent file doesn't have SCORE: header")

    scoreline = line1

    sp = scoreline.split()
    if ( len(sp) < 2 or sp[1] != "score" and sp[1] != "total_score" ):
        eprint("silent_tools: Warning! First score is not \"score\"! Rosetta won't like this!")

    if ( return_f ):
        return scoreline, f

    f.close()

    return scoreline


def build_silent_index(file, accept_garbage=False):

    scoreline = assert_is_silent_and_get_scoreline(file, accept_garbage=accept_garbage)


    # I'm sorry. If you put description in the name of your pose, it will disappear
    lines = cmd2("command grep -a --byte-offset '^SCORE:' %s | grep -va description | tr -d '\r' | awk '{print $1,$NF}'"%file)[0].strip().split("\n")


    # with open("tmp", "w") as f:
    #     f.write("\n".join(lines))
    # with open("tmp") as f:
    #     lines = f.read().split("\n")

    index = defaultdict(lambda : {}, {})
    order = []
    orig_order = []
    unique_tags = True

    dup_index = {}

    for line in lines:
        try:
            # eprint(line)
            sp = line.strip().split()

            # this might seem like a weird test, but it catches when awk only gets 1 field
            if ( sp[0] == sp[1] ):
                offset = 0 if len(order) == 0 else index[order[-1]]['seek']
                eprint("silent_tools: corruption: file_offset: %i"%(offset))
                continue

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
        except:
            offset = 0 if len(order) == 0 else index[order[-1]]['seek']
            eprint("silent_tools: corruption: file_offset: %i -- %s"%(offset, line))

    size = file_size(file)

    silent_index = {"index":index, "tags":order, "orig_tags":orig_order, "scoreline":scoreline, "size":size, 
                    "unique_tags":unique_tags, "version":SILENT_INDEX_VERSION}

    sequence = "A"
    silent_type = "UNKNOWN"
    if ( len(order) > 0 ):
        try:
            structure = get_silent_structure(file, silent_index, order[0])
            sequence = "".join(get_sequence_chunks(structure))
            silent_type = detect_silent_type(structure)
        except:
            eprint("Failed to get sequence. Please tell Brian")

    silent_index['sequence'] = sequence
    silent_index['silent_type'] = silent_type


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

def silent_header_fix_corrupt(silent_index):
    use_scoreline = silent_index['scoreline']

    sp = use_scoreline.split()
    if ( len(sp) < 2 or sp[1] != "score" and sp[1] != "total_score" ):
        use_scoreline = "SCORE: score description"

    return silent_header_slim(silent_index['sequence'], use_scoreline, silent_index['silent_type'])

def silent_header(silent_index):
    return silent_header_slim(silent_index['sequence'], silent_index['scoreline'], silent_index['silent_type'])

def silent_header_slim(sequence, scoreline, silent_type):
    header = "SEQUENCE: %s\n%s\n"%(sequence, scoreline.strip())
    if ( silent_type != "UNKNOWN" ):
        header += "REMARK %s SILENTFILE\n"%silent_type
    return header



def get_sequence_chunks(structure, tag="FIXME"):

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

def get_chain_ids(structure, tag="FIXME"):

    resnum_line = None
    for line in structure:
        if ( line.startswith("RES_NUM") ):
            resnum_line = line
            break

    if ( resnum_line is None ):
        eprint("silent_tools: no RES_NUM for tag %s"%tag)
        return ""


    parts = line.split()
    usable_parts = [x for x in parts if ":" in x]

    chain_ids = ""
    for part in usable_parts:
        idd, rangee = part.split(":")
        assert(len(idd) == 1)

        start, end = [int(x) for x in rangee.split('-')]

        chain_ids += idd*(end-start + 1)

    return chain_ids


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


import importlib.util
package_name = 'numba'
spec = importlib.util.find_spec(package_name)
if not spec is None:
    
    from numba import njit


    @njit(fastmath=True)
    def code_from_6bit(_8bit):

        if ( ( _8bit >= 65) and (_8bit <= 90) ): return _8bit - 65
        if ( ( _8bit >= 97) and (_8bit <= 122) ): return _8bit - 97 + 26
        if ( ( _8bit >= 48) and (_8bit <= 57) ): return _8bit - 48 + 52
        if (   _8bit == 43  ): return 62
        return 63


    @njit(fastmath=True)
    def decode_32_to_24( i0, i1, i2, i3 ):
        i0 = code_from_6bit( i0 )
        i1 = code_from_6bit( i1 )
        i2 = code_from_6bit( i2 )
        i3 = code_from_6bit( i3 )

        o0 = 0xFF & (i0 | (i1 << 6))
        o1 = 0xFF & ((i1 >> 2) | (i2 << 4))
        o2 = 0xFF & ((i3 << 2) | (i2 >> 4))

        return o0, o1, o2

    scr = np.zeros(1000, np.byte)
    def decode6bit( jar ):
        return numba_decode6bit( jar.encode(), scr )

    @njit(fastmath=True)
    def numba_decode6bit( jar, ba ):
        ba_len = 0

        this_str = np.zeros(4, np.byte)

        valid_bits = 0
        i = 0
        while ( i < len(jar) ):

            this_str[0] = 0
            this_str[1] = 0
            this_str[2] = 0
            this_str[3] = 0

            j = 0
            while ( i < len(jar) and j < 4 ):
                this_str[j] = jar[i]
                i += 1
                j += 1
                valid_bits += 6

            # print(this_str)
            o0, o1, o2 = decode_32_to_24(this_str[0], this_str[1], this_str[2], this_str[3])
            # print(bytess)

            ba[ba_len] = o0
            ba[ba_len+1] = o1
            ba[ba_len+2] = o2
            ba_len += 3
        valid_bytes = int( valid_bits / 8 )
        ba = ba[:valid_bytes]
        assert(len(ba) % 4 == 0)
        return ba




_float_packer_by_len = None
def silent_line_to_atoms(line):
    global _float_packer_by_len
    if ( _float_packer_by_len is None ):
        _float_packer_by_len = []
        for i in range(1000):
            _float_packer_by_len.append(struct.Struct("f"*(i)))


    ba = decode6bit( line )

    float_packer = _float_packer_by_len[len(ba)//4] #struct.Struct("f"*(len(ba)//4))

    floats = float_packer.unpack(ba)

    assert(len(floats) % 3 == 0)

    return np.array(floats).reshape(-1, 3)


def get_chains_mask(chunks, chains):
    sequence = "".join(chunks)
    if ( chains is None ):
        mask = np.ones(len(sequence))
    else:
        mask = np.zeros(len(sequence))
        for chain in chains:
            lb = np.sum([len(chunk) for chunk in chunks[:chain]]).astype(int)
            ub = np.sum([len(chunk) for chunk in chunks[:chain+1]]).astype(int)
            mask[lb:ub] = True
    return mask


def sketch_get_atoms_by_residue(structure, chains=None):

    chunks = get_sequence_chunks(structure)
    sequence = "".join(chunks)
    if ( sequence is None ):
        return None

    mask = get_chains_mask(chunks, chains)

    residues = []

    ires = -1
    for line in structure:

        if ( len(line) == 0 ):
            continue

        if ( line[0] not in "EHL" ):
            continue

        # Ok, so we're going to use some really crappy detection here
        sp = line.split()
        if ( len(sp) != 2 ):
            continue

        ires += 1
        if ( not mask[ires] ):
            continue

        binary = sp[0][1:]

        residues.append( silent_line_to_atoms( binary ) )

    # print(np.sum(mask), len(residues))
    assert(np.sum(mask) == len(residues))
    return residues

def sketch_get_atoms(structure, atom_nums, chains=None):

    atoms_by_res = sketch_get_atoms_by_residue(structure, chains)

    final = []
    for residue in atoms_by_res:
        try:
            final.append(residue[atom_nums])
        except:
            arr = []
            for atom_num in atom_nums:
                try:
                    arr.append(residue[atom_num])
                except:
                    arr.append(np.array([np.nan, np.nan, np.nan]))
            final.append(np.array(arr))

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


def sketch_get_ncac_protein_struct(structure):

    sequence = "".join(get_sequence_chunks(structure))

    ncac = []

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
            nx = float(sp[2])
            ny = float(sp[3])
            nz = float(sp[4])
            cax = float(sp[5])
            cay = float(sp[6])
            caz = float(sp[7])
            cx = float(sp[8])
            cy = float(sp[9])
            cz = float(sp[10])
        except:
            continue
        ncac.append([nx, ny, nz])
        ncac.append([cax, cay, caz])
        ncac.append([cx, cy, cz])

        assert(seqpos*3 == len(ncac))

    assert(len(ncac) == len(sequence)*3)

    return np.array(ncac)



















