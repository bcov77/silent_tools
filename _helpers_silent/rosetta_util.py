import sys
import os
import re
import struct

from collections import defaultdict

import numpy as np

#########################
# These functions all used to be in silent_tools.py
#########################


def write_pdb_atoms(atoms, sequence, atom_names):
    lines = []
    assert(len(atoms) / len(sequence) == len(atom_names))

    for i in range(len(sequence)):
        try:
            name3 = name1_to_name3[sequence[i]]
        except:
            name3 = "UNK"

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

def inner_get_silent_res_data(ba):
    line = "L"

    iters = int(np.ceil(len(ba) / 3))

    len_ba = len(ba)
    for i in range(iters):
        i0 = 0
        i1 = 0
        i2 = 0
        i0 = ba[i*3+0]
        if ( i*3 + 1 < len_ba ):
            i1 = ba[i*3+1]
        if ( i*3+2 < len_ba ):
            i2 = ba[i*3+2]

        line += encode_24_to_32(i0, i1, i2)

    return line

def code_to_6bit(byte):
    return silent_chars[byte]

def encode_24_to_32(i0, i1, i2):
    return code_to_6bit( i0 & 63 ) + \
            code_to_6bit( ((i1 << 2) | (i0 >> 6)) & 63 ) + \
            code_to_6bit( ((i1 >> 4) | ((i2 << 4) & 63)) & 63 ) + \
            code_to_6bit( i2 >> 2 )


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

    @njit(fastmath=True, cache=True)
    def inner_get_silent_res_data(ba):
        line = "L"

        iters = int(np.ceil(len(ba) / 3))

        len_ba = len(ba)
        for i in range(iters):
            i0 = 0
            i1 = 0
            i2 = 0
            i0 = ba[i*3+0]
            if ( i*3 + 1 < len_ba ):
                i1 = ba[i*3+1]
            if ( i*3+2 < len_ba ):
                i2 = ba[i*3+2]

            line += encode_24_to_32(i0, i1, i2)

        return line

    @njit(fastmath=True, cache=True)
    def code_to_6bit(byte):
        return silent_chars[byte]

    @njit(fastmath=True, cache=True)
    def encode_24_to_32(i0, i1, i2):
        return code_to_6bit( i0 & 63 ) + \
                code_to_6bit( ((i1 << 2) | (i0 >> 6)) & 63 ) + \
                code_to_6bit( ((i1 >> 4) | ((i2 << 4) & 63)) & 63 ) + \
                code_to_6bit( i2 >> 2 )


_float_packer_0 = struct.Struct("f")
def get_silent_res_data(coords):

    ba = bytearray()
    for coord in coords:
        ba += _float_packer_0.pack(coord)

    return inner_get_silent_res_data(ba)


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












#########################################################################################################################################
#                                               Loading PDBs from silent files without Rosetta
#########################################################################################################################################

rosetta_aa_order = {
'ALA': [' N  ',' CA ',' C  ',' O  ',' CB ',' H  ',' HA ','1HB ','2HB ','3HB '],
'CYS': [' N  ',' CA ',' C  ',' O  ',' CB ',' SG ',' H  ',' HA ','1HB ','2HB ',' HG '],
'ASP': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' OD1',' OD2',' H  ',' HA ','1HB ','2HB '],
'GLU': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ',' OE1',' OE2',' H  ',' HA ','1HB ','2HB ','1HG ','2HG '],
'PHE': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD1',' CD2',' CE1',' CE2',' CZ ',' H  ',' HA ','1HB ','2HB ',' HD1',' HD2',' HE1',' HE2',' HZ '],
'GLY': [' N  ',' CA ',' C  ',' O  ',' H  ','1HA ','2HA '],
'HIS': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' ND1',' CD2',' CE1',' NE2',' H  ',' HA ','1HB ','2HB ',' HD2',' HE1',' HE2'],
'HIS_D': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' ND1',' CD2',' CE1',' NE2',' H  ',' HA ','1HB ','2HB ',' HD1',' HD2',' HE1'],
'ILE': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG1',' CG2',' CD1',' H  ',' HA ',' HB ','1HG1','2HG1','1HG2','2HG2','3HG2','1HD1','2HD1','3HD1'],
'LYS': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ',' CE ',' NZ ',' H  ',' HA ','1HB ','2HB ','1HG ','2HG ','1HD ','2HD ','1HE ','2HE ','1HZ ','2HZ ','3HZ '],
'LEU': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD1',' CD2',' H  ',' HA ','1HB ','2HB ',' HG ','1HD1','2HD1','3HD1','1HD2','2HD2','3HD2'],
'MET': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' SD ',' CE ',' H  ',' HA ','1HB ','2HB ','1HG ','2HG ','1HE ','2HE ','3HE '],
'ASN': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' OD1',' ND2',' H  ',' HA ','1HB ','2HB ','1HD2','2HD2'],
'PRO': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ',' NV ',' HA ','1HB ','2HB ','1HG ','2HG ','1HD ','2HD '],
'GLN': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ',' OE1',' NE2',' H  ',' HA ','1HB ','2HB ','1HG ','2HG ','1HE2','2HE2'],
'ARG': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD ',' NE ',' CZ ',' NH1',' NH2',' H  ',' HA ','1HB ','2HB ','1HG ','2HG ','1HD ','2HD ',' HE ','1HH1','2HH1','1HH2','2HH2'],
'SER': [' N  ',' CA ',' C  ',' O  ',' CB ',' OG ',' H  ',' HA ','1HB ','2HB ',' HG '],
'THR': [' N  ',' CA ',' C  ',' O  ',' CB ',' OG1',' CG2',' H  ',' HA ',' HB ',' HG1','1HG2','2HG2','3HG2'],
'VAL': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG1',' CG2',' H  ',' HA ',' HB ','1HG1','2HG1','3HG1','1HG2','2HG2','3HG2'],
'TRP': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD1',' CD2',' NE1',' CE2',' CE3',' CZ2',' CZ3',' CH2',' H  ',' HA ','1HB ','2HB ',' HD1',' HE1',' HE3',' HZ2',' HZ3',' HH2'],
'TYR': [' N  ',' CA ',' C  ',' O  ',' CB ',' CG ',' CD1',' CD2',' CE1',' CE2',' CZ ',' OH ',' H  ',' HA ','1HB ','2HB ',' HD1',' HD2',' HE1',' HE2',' HH '],

# DNA
' DT': [' P  ', ' OP2', ' OP1', " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " C1'", ' N1 ', ' C2 ', ' O2 ', ' N3 ', ' C4 ', ' O4 ', ' C5 ', ' C7 ', ' C6 ', "H5''", " H5'", " H4'", " H3'", "H2''", " H2'", " H1'", ' H3 ', ' H71', ' H72', ' H73', ' H6 '],
' DG': [' P  ', ' OP2', ' OP1', " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " C1'", ' N9 ', ' C4 ', ' N3 ', ' C2 ', ' N1 ', ' C6 ', ' C5 ', ' N7 ', ' C8 ', ' N2 ', ' O6 ', "H5''", " H5'", " H4'", " H3'", "H2''", " H2'", " H1'", ' H1 ', ' H8 ', ' H22', ' H21'],
' DA': [' P  ', ' OP2', ' OP1', " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " C1'", ' N9 ', ' C4 ', ' N3 ', ' C2 ', ' N1 ', ' C6 ', ' C5 ', ' N7 ', ' C8 ', ' N6 ', "H5''", " H5'", " H4'", " H3'", "H2''", " H2'", " H1'", ' H2 ', ' H8 ', ' H61', ' H62'],
' DC': [' P  ', ' OP2', ' OP1', " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " C1'", ' N1 ', ' C2 ', ' O2 ', ' N3 ', ' C4 ', ' N4 ', ' C5 ', ' C6 ', "H5''", " H5'", " H4'", " H3'", "H2''", " H2'", " H1'", ' H42', ' H41', ' H5 ', ' H6 '],

# RNA
'  C': [' P  ', ' OP2', ' OP1', " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C1'", " C2'", " O2'", ' N1 ', ' C2 ', ' O2 ', ' N3 ', ' C4 ', ' N4 ', ' C5 ', ' C6 ', " H5'", "H5''", " H4'", " H3'", " H1'", " H2'", "HO2'", ' H42', ' H41', ' H5 ', ' H6 '],
'  G': [' P  ', ' OP2', ' OP1', " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C1'", " C2'", " O2'", ' N1 ', ' C2 ', ' N2 ', ' N3 ', ' C4 ', ' C5 ', ' C6 ', ' O6 ', ' N7 ', ' C8 ', ' N9 ', " H5'", "H5''", " H4'", " H3'", " H1'", " H2'", "HO2'", ' H1 ', ' H22', ' H21', ' H8 '],
'  U': [' P  ', ' OP2', ' OP1', " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C1'", " C2'", " O2'", ' N1 ', ' C2 ', ' O2 ', ' N3 ', ' C4 ', ' O4 ', ' C5 ', ' C6 ', " H5'", "H5''", " H4'", " H3'", " H1'", " H2'", "HO2'", ' H3 ', ' H5 ', ' H6 '],
'  A': [' P  ', ' OP2', ' OP1', " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C1'", " C2'", " O2'", ' N1 ', ' C2 ', ' N3 ', ' C4 ', ' C5 ', ' C6 ', ' N6 ', ' N7 ', ' C8 ', ' N9 ', " H5'", "H5''", " H4'", " H3'", " H1'", " H2'", "HO2'", ' H2 ', ' H61', ' H62', ' H8 ']

}

synonyms = {
    " NV ":" N  ", # PRO
    "CAV ":" CA ", # PRO
}

# Rosetta, why did you do this???
#   I think it's because it needs to actually be 3 letters
silent_name3_to_pdb_name3 = {
    'GUA':' DG',
    'ADE':' DA',
    'CYT':' DC',
    'THY':' DT',
    'RAD':'  A',
    'RCY':'  C',
    'RGU':'  G',
    'URA':'  U'
}

pdb_name3_to_silent_name3 = {}
for key, val in silent_name3_to_pdb_name3.items():
    pdb_name3_to_silent_name3[val] = key



needs_annotated_name = {
    ' DG',
    ' DA',
    ' DC',
    'HIS_D'
}

name3_to_name1 = {
    "ARG":"R",
    "LYS":"K",
    "ASN":"N",
    "ASP":"D",
    "GLU":"E",
    "GLN":"Q",
    "HIS":"H",
    "PRO":"P",
    "TYR":"Y",
    "TRP":"W",
    "SER":"S",
    "THR":"T",
    "GLY":"G",
    "ALA":"A",
    "MET":"M",
    "CYS":"C",
    "PHE":"F",
    "LEU":"L",
    "VAL":"V",
    "ILE":"I",

# DNA
    " DA":"a",
    " DT":"t",
    " DG":"g",
    " DC":"c",
# RNA -- must come second because RNA technically has "dibs" on the name1 letters
    "  A":"a",
    "  U":"u",
    "  G":"g",
    "  C":"c"
}

name1_to_name3 = {}
for key, val in name3_to_name1.items():
    name1_to_name3[val] = key

name3_to_name1['HIS_D'] = 'H'


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

_atom_record_format = (
    "ATOM  {atomi:5d} {atomn:^4}{idx:^1}{resn:3s} {chain:1}{resi:4d}{insert:1s}   "
    "{x:8s}{y:8s}{z:8s}{occ:6.2f}{b:6.2f}           {elem:3s}\n"
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
        b=0,
        elem=''
):
    x = f'{x:8.3f}'[:8] # this isn't exactly how rosetta handles huge coords but what are you going to do anyways
    y = f'{y:8.3f}'[:8]
    z = f'{z:8.3f}'[:8]
    return _atom_record_format.format(**locals())

_disulfide_record_format = (
    "SSBOND     CYS {chain_a:1s} {resi_a:4d}    CYS {chain_b:1s} {resi_b:4d}                                           \n"
    )
def format_disulfide(
    chain_a='A',
    resi_a=0,
    chain_b='A',
    resi_b=0
):
    return _disulfide_record_format.format(**locals())


def get_atom_order_from_variants(name3, variants):
    atom_order = list(rosetta_aa_order[name3])

    for variant in variants:
        if variant == "NtermProteinFull":
            if name3 == "PRO":
                NV_idx = atom_order.index(" NV ")
                atom_order.insert(NV_idx+1, "CAV ")

                for new_atom in ["1H  ", "2H  "]:
                    atom_order.append(new_atom)
            else:
                H_idx = atom_order.index(" H  ")
                atom_order.pop(H_idx)
                for new_atom in reversed(["1H  ", "2H  ", "3H  "]):
                    atom_order.insert(H_idx, new_atom)
            continue

        if variant == "CtermProteinFull":
            O_idx = atom_order.index(" O  ")
            atom_order.insert(O_idx+1," OXT")
            continue

        if variant == "protein_cutpoint_lower":
            O_idx = atom_order.index(" O  ")
            atom_order.insert(O_idx+1, "OVL1")
            atom_order.insert(O_idx+2, "OVL1")
            continue

        if variant == "protein_cutpoint_upper":
            O_idx = atom_order.index(" O  ")
            atom_order.insert(O_idx+1, "OVU1")
            continue

        if variant == "disulfide":
            atom_order.remove(" HG ")
            continue

        # RNA
        if variant == 'UpperRNA':
            # no change
            continue
        if variant == 'LowerRNA':
            # no change
            continue

        # DNA
        if variant == "LowerDNA":
            # There's no change
            continue

        if variant == "UpperDNA":
            idx = atom_order.index("H2''")
            atom_order.insert(idx, "HO3'")
            continue


        assert False, "Warning! Unknown variant:" + variant

    return atom_order



def structure_to_pdb(structure, renumber_resno_increasing=False):

    annotated_sequence = None
    res_num = None
    chain_endings = None
    atom_lines = []
    pdb_infos = defaultdict(list)
    disulfides = []

    for line in structure:
        line = line.strip()
        if len(line) == 0:
            continue
        sp = line.split()
        if line.startswith("ANNOTATED_SEQUENCE:"):
            annotated_sequence = sp[1]
            continue
        if line.startswith("CHAIN_ENDINGS"):
            chain_endings = sp[1:-1]
            continue
        if line.startswith("RES_NUM"):
            res_num = sp[1:-1]
            continue
        if line.startswith("REMARK PDBinfo-LABEL:"):
            seqpos = int(sp[2])
            pdb_infos[seqpos] += sp[3:]
            continue
        if line.startswith("NONCANONICAL_CONNECTION"):
            if ((sp[2] == "N" and sp[4] == "C") or (sp[2] == "C" and sp[4] == "N")) and abs(int(sp[3]) - int(sp[1])) == 1:
                continue
            assert(sp[2] == "SG")
            assert(sp[4] == "SG")
            disulfides.append([int(sp[1]), int(sp[3])])
            continue
        if line.startswith("FOLD_TREE"):
            continue
        if line.startswith("SCORE:"):
            continue
        if line.startswith("REMARK BINARY SILENTFILE"):
            continue
        if line.startswith("RT"):
            continue
        if line[0] in "ELH" and len(sp) == 2:
            atom_lines.append(sp[0])
            continue
        if line.startswith("REMARK"):
            continue
        eprint("Warning! Unknown silentfile line:", line)


    assert(annotated_sequence is not None)

    annotated_sequence = [x.groups()[0] for x in re.finditer(r"([A-Za-z]([\[][^\]]+[\]])?)", annotated_sequence)]
    assert(len(annotated_sequence) == len(atom_lines))

    L = len(annotated_sequence)

    if ( res_num is None ):
        res_num = ["A:1-%s"%L]
    if ( chain_endings is None ):
        chain_endings = []

    # make it A:1-10 B:11-20 C:21-30 etc
    if renumber_resno_increasing:
        res_num = []
        last_ending = 0
        for iending, ending in enumerate(chain_endings + [L]):
            ending = int(ending)
            letter = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[iending]
            res_num.append(f"{letter}:{last_ending+1}-{ending}")
            last_ending = ending


    chain_endings = set([int(x) for x in chain_endings])


    expanded_res_num = [None]*L
    expanded_chain = [None]*L
    seqpos0 = -1
    for elem in res_num:
        sp = elem.split(":") # we're dropping segment_id if it exists
        chain = sp[0]
        rang = sp[-1]
        if '-' in rang:
            lb, ub = [int(x) for x in rang.split("-")]
        else:
            lb = ub = int(rang)
        for num in range(lb, ub+1):
            seqpos0 += 1
            expanded_res_num[seqpos0] = num
            expanded_chain[seqpos0] = chain
    assert expanded_chain[-1] is not None




    pdb_lines = []

    disulfides = [list(sorted(x)) for x in disulfides]
    disulfides = list(sorted(disulfides, key=lambda x: x[0]))

    for disulfide in disulfides:
        seqpos_a, seqpos_b = disulfide
        pdb_lines.append(format_disulfide(
            chain_a=expanded_chain[seqpos_a-1],
            resi_a=expanded_res_num[seqpos_a-1],
            chain_b=expanded_chain[seqpos_b-1],
            resi_b=expanded_res_num[seqpos_b-1]
            ))



    atomno = 0
    for seqpos0 in range(L):
        annotated_name = annotated_sequence[seqpos0]
        name1 = annotated_name[0]
        name3 = name1_to_name3[name1]


        variants = []
        if ( len(annotated_name) > 1 ):
            var_split = annotated_name[2:-1].split(":")
            variants = var_split[1:]

            # Because the silent file uses a different naming convention for DNA/RNA...
            if var_split[0] in silent_name3_to_pdb_name3:
                name3 = silent_name3_to_pdb_name3[var_split[0]]
            else:
                name3 = var_split[0]

        atom_order = get_atom_order_from_variants(name3, variants)

        atom_coords = silent_line_to_atoms(atom_lines[seqpos0][1:])

        assert(len(atom_coords) == len(atom_order))


        for iatom in range(len(atom_order)):
            atom_name = atom_order[iatom]
            if "V" in atom_name:
                continue
            atomno += 1

            element = atom_name[1]
            if atom_name[0] not in ' 0123456789':
                element = atom_name[0]

            pdb_lines.append(format_atom(
                atomi=atomno, 
                atomn=atom_name,
                resn=name3[:3],
                chain=expanded_chain[seqpos0],
                resi=expanded_res_num[seqpos0],
                x=atom_coords[iatom][0],
                y=atom_coords[iatom][1],
                z=atom_coords[iatom][2],
                elem=element
                ))


        if seqpos0 + 1 in chain_endings:
            atomno += 1
            pdb_lines.append("TER                                                                             \n")


    if pdb_lines[-1] != "TER                                                                             \n":
        pdb_lines.append("TER                                                                             \n")

    for seqpos in sorted(pdb_infos):
        pdb_lines.append(f'REMARK PDBinfo-LABEL: {seqpos:4d} ' + " ".join(pdb_infos[seqpos]) + "\n")

    return pdb_lines





#################################################################################################
#                                Writing Functions
#################################################################################################


# converts AAAAAAAABBBBBBB to 1-7:A 8-14:B
def chain_ids_to_silent_format(chain_ids, res_nums=None):
    if res_nums is None:
        res_nums = list(range(1, len(chain_ids)+1))
    res_nums = res_nums + [None]
    parts = []
    cur_letter = None
    cur_start = None
    start_no = None
    for i, letter in enumerate(chain_ids + "\n"): # chain id can never be \n
        if letter != cur_letter or (start_no + i - cur_start != res_nums[i]):
            if cur_letter != None:
                parts.append("%s:%i-%i"%(cur_letter, start_no, start_no + i - cur_start - 1))
            cur_letter = letter
            cur_start = i
            start_no = res_nums[i]
    return " ".join(parts)



def get_stub_from_n_ca_c(n, ca, c):
    e1 = ca - n
    e1 /= np.linalg.norm(e1)

    e3 = np.cross( e1, c - n )
    e3 /= np.linalg.norm(e3)

    e2 = np.cross( e3, e1 )

    stub = np.identity(4)
    stub[:3,0] = e1
    stub[:3,1] = e2
    stub[:3,2] = e3
    stub[:3,3] = ca

    return stub

# this function is crazy
# this is just how it works
# /home/bcov/pyrosetta_experiments/write_silent_no_rosetta/determine_jump_atoms.py
def get_jump_stub(atom_dict, is_pro=False):
    if ' CA ' in atom_dict:
        N = np.array(atom_dict[' N  '])
        CA = np.array(atom_dict[' CA '])
        if ( ' H  ' in atom_dict ):
            H = np.array(atom_dict[' H  '])
        else:
            if is_pro:
                H = np.array(atom_dict[' C  '])
            else:
                H = np.array(atom_dict['1H  '])

        return get_stub_from_n_ca_c(N, CA, H)

    if " C5'" in atom_dict:
        return get_stub_from_n_ca_c( np.array(atom_dict[" O5'"]), 
                                     np.array(atom_dict[" C5'"]),
                                     np.array(atom_dict[" P  "]))


def parse_pdb_into_needed_format(pdb_lines):
    # PDBs must be parsed into this format
    chain_ends = []      # 1-indexed last residue number of each chain
    aa_atom_dicts = []   # dictionary of xyz coords for each res.
    name3s = []          # name3s
    chain_letters = []   # string or list of chain letters for each aa
    res_nums = []        # the text residue number from the pdb
    in_a_disulfides = [] # is this in a disulfide?
    disulfides = []      # disulfides in 1-indexed continuous numerical space [a, b]
    pdb_info_labels = {} # int -> list(string)

    raw_disulfides = []     # temporary chain_no -- chain_no that we parse from SSBOND
    chain_nos = []          # temporary chain-number used for disulfides
    in_a_disulfide_set = set() # temporary which res_nos are in disulfides

    last_ident = ''

    for line in pdb_lines:
        line = line.strip()

        if line.startswith('SSBOND'):
            sp = line.split()
            assert len(sp) >= 7, 'Weird SSBOND line: ' + line 
            _, _, chain_a, resi_a, _, chain_b, resi_b = sp[:7]
            raw_disulfides.append([chain_a + resi_a, chain_b + resi_b])
            disulfides.append([None, None])
            in_a_disulfide_set.add(chain_a + resi_a)
            in_a_disulfide_set.add(chain_b + resi_b)
            continue

        if line.startswith("REMARK PDBinfo-LABEL:"):
            sp = line.split()
            seqpos = int(sp[2])
            if seqpos in pdb_info_labels:
                pdb_info_labels[seqpos] += sp[3:]
            else:
                pdb_info_labels[seqpos] = sp[3:]

        if not line.startswith('ATOM') or len(line) < 54:
            if len(aa_atom_dicts) > 0:
                if len(chain_ends) == 0 or chain_ends[-1] != len(aa_atom_dicts):
                    chain_ends.append(len(aa_atom_dicts))
            last_ident = ''

            continue

        ident = line[17:27]

        x = float(line[30:30+8])
        y = float(line[38:38+8])
        z = float(line[46:46+8])
        name3 = line[17:20]
        chain = line[21]
        atom = line[12:16]
        res_num = int(line[22:26].strip())

        if ident != last_ident:
            aa_atom_dicts.append({})
            name3s.append(name3)
            chain_letters.append(chain)
            chain_no = ident[3:].replace(' ', '')
            chain_nos.append(chain_no)
            last_ident = ident
            res_nums.append(res_num)

            seqpos_ros = len(name3s)
            if chain_no in in_a_disulfide_set:
                for i_dis, (a, b) in enumerate(raw_disulfides):
                    if a == chain_no:
                        disulfides[i_dis][0] = seqpos_ros
                    if b == chain_no:
                        disulfides[i_dis][1] = seqpos_ros

            in_a_disulfides.append( chain_no in in_a_disulfide_set )


        aa_atom_dicts[-1][atom] = [x, y, z]

    for a, b in disulfides:
        assert not a is None and not b is None


    return chain_ends, aa_atom_dicts, name3s, chain_letters, res_nums, in_a_disulfides, disulfides, pdb_info_labels


def fmt(value):
    if isinstance(value, str):
        return value
    else:
        return '%.3f'%value

def parsed_pdb_to_silent( chain_ends, aa_atom_dicts, name3s, chain_letters, res_nums, in_a_disulfides, disulfides, tag,
                                                                            write_header=False, score_dict={}, pdb_info_labels={},
                                                                            allow_missing_atoms=False,
                                                                            allow_extra_atoms=False):
    assert len(aa_atom_dicts) > 0, 'Error! There are no atoms in your pdb?'

    # annotated_sequence has to come really early which is why we split this
    silent = ''
    post_silent = ''

    if write_header:
        silent += 'SEQUENCE: A\n'
        silent += f'SCORE:     score {" ".join(score_dict.keys())} description\n'
        silent += 'REMARK BINARY SILENTFILE\n'

    # SCORE:
    scores_string = ' '.join([fmt(score) for score in score_dict.values()])
    silent += f'SCORE:     0.000 {scores_string}        {tag}\n'

    # ANNOTATED_SEQUENCE: goes here

    # PDBinfo-LABEL
    for seqpos in sorted(list(pdb_info_labels)):
        labels = pdb_info_labels[seqpos]
        post_silent += f'REMARK PDBinfo-LABEL:{seqpos:5d} {" ".join(labels)}\n'

    # FOLD_TREE
    ft_parts = ['FOLD_TREE']
    last_end = 0
    for iend, end in enumerate(chain_ends):
        if iend > 0:
            ft_parts.append(f'EDGE 1 {last_end+1} {iend}')
        ft_parts.append(f'EDGE {last_end+1} {end} -1')
        last_end = end
    ft_parts.append(tag)

    post_silent += '  '.join(ft_parts) + '\n'

    # CHAIN_ENDINGS
    if len(chain_ends) > 1:
        post_silent += f'CHAIN_ENDINGS {" ".join(str(x) for x in chain_ends[:-1])} {tag}\n'

    # RES_NUM
    post_silent += f'RES_NUM {chain_ids_to_silent_format("".join(chain_letters), res_nums)} {tag}\n'

    # NONCANONICAL_CONNECTION
    for a, b in disulfides:
        post_silent += f'NONCANONICAL_CONNECTION: {a} SG {b} SG\n'

    # RT
    stub1 = get_jump_stub(aa_atom_dicts[0], is_pro=name3s[0] == 'PRO')
    for end in chain_ends[:-1]:
        stub2 = get_jump_stub(aa_atom_dicts[end], is_pro=name3s[end] == 'PRO')
        jump_rt = np.linalg.inv(stub1) @ stub2

        post_silent += f'RT {" ".join("%.8f"%x for x in list(jump_rt[:3,:3].flat) + list(jump_rt[:3,3].flat))} {tag}\n'


    # residue lines
    missing_atom = [0, 0, 0]
    annotated_seq = ''

    for ires, (d, name3, in_a_disulfide) in enumerate(zip(aa_atom_dicts, name3s, in_a_disulfides)):

        assert name3 in rosetta_aa_order, 'Error! Unknown name3: ' + name3
        if name3 == 'HIS':
            if ' HD1' in d:
                name3 = 'HIS_D'

        variants = []
        if '1H  ' in d:
            variants.append('NtermProteinFull')
        if ' OXT' in d:
            variants.append("CtermProteinFull")
        if in_a_disulfide:
            assert name3 == 'CYS'
            variants.append('disulfide')

        # there's no way to detect LowerDNA, but it isn't necessary
        if "HO3'" in d:
            variants.append('UpperDNA')

        # RNA doesn't even change anything with LowerDNA and UpperDNA

        atom_order = get_atom_order_from_variants(name3, variants)

        atom_used = {}
        for atom_name in d:
            atom_used[atom_name] = False

        messages = []
        no_atom_errors = True
        res_coords = []
        for atom_name in atom_order:
            og_atom_name = atom_name

            not_found = atom_name not in d
            if not_found and atom_name in synonyms:
                atom_name = synonyms[atom_name]
                not_found = atom_name not in d

            if not_found:
                messages.append(f'Residue {ires+1} {name3} missing atom "{og_atom_name}"')
                eprint(messages[-1])
                if not allow_missing_atoms:
                    no_atom_errors = False
                res_coords += missing_atom
            else:
                res_coords += d[atom_name]
                atom_used[atom_name] = True

        for atom_name in atom_used:
            if not atom_used[atom_name]:
                messages.append(f'Residue {ires+1} {name3} extra atom "{atom_name}"')
                eprint(messages[-1])
                if not allow_extra_atoms:
                    no_atom_errors = False

        assert no_atom_errors, '\n'.join(messages)


        post_silent += get_silent_res_data(res_coords) + f' {tag}\n'

        annotated_seq += name3_to_name1[name3]
        if name3 in needs_annotated_name or len(variants) > 0:

            write_name3 = name3
            if name3 in pdb_name3_to_silent_name3:
                write_name3 = pdb_name3_to_silent_name3[name3]
            to_add = '[' + ':'.join( [write_name3] + variants ) + ']'
            annotated_seq += to_add


    silent += "ANNOTATED_SEQUENCE: %s %s\n"%(annotated_seq, tag)
    silent += post_silent

    return silent

# score_dict is string -> float or string
# pdb_info_labels is int -> list(string)
def pdb_to_structure(pdb_lines, tag, write_header=True, score_dict={}, pdb_info_labels={}):

    chain_ends, aa_atom_dicts, name3s, chain_letters, res_nums, in_a_disulfides, disulfides, more_info_labels = parse_pdb_into_needed_format(pdb_lines)

    # if you add into the pdb_info_labels object you modify global state
    for key in pdb_info_labels:
        if key in more_info_labels:
            more_info_labels[key] += pdb_info_labels[key]
        else:
            more_info_labels[key] = pdb_info_labels[key]


    silent = parsed_pdb_to_silent( chain_ends, aa_atom_dicts, name3s, chain_letters, res_nums, in_a_disulfides, disulfides, tag,
                                                write_header=write_header, score_dict=score_dict, pdb_info_labels=more_info_labels)

    return silent








