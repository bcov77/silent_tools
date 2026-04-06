import sys
import os
import re
import struct
import base64

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

atom_type_to_element = {
    'CNH2':'C',
    'COO':'C',
    'CH0':'C',
    'CH1':'C',
    'CH2':'C',
    'CH3':'C',
    'aroC':'C',
    'Ntrp':'N',
    'Nhis':'N',
    'NtrR':'N',
    'NH2O':'N',
    'Nlys':'N',
    'Narg':'N',
    'Npro':'N',
    'OH':'O',
    'OW':'O',
    'ONH2':'O',
    'OOC':'O',
    'Oaro':'O',
    'Oet2':'O',
    'Oet3':'O',
    'S':'S',
    'SH1':'S',
    'Nbb':'N',
    'CAbb':'C',
    'CObb':'C',
    'OCbb':'O',
    'Phos':'P',
    'Pbb':'P',
    'Hpol':'H',
    'HS':'H',
    'Hapo':'H',
    'Haro':'H',
    'HNbb':'H',
    'Hwat':'H',
    'Owat':'O',
    'Opoint':'O',
    'HOH':'O',
    'Bsp2':'B',
    'F':'F',
    'Cl':'CL',
    'Br':'BR',
    'I':'I',
    'Zn2p':'ZN',
    'Co2p':'CO',
    'Cu2p':'CU',
    'Fe2p':'FE',
    'Fe3p':'FE',
    'Mg2p':'MG',
    'Ca2p':'CA',
    'Pha':'P',
    'OPha':'O',
    'OHha':'O',
    'Hha':'H',
    'CO3':'C',
    'OC3':'O',
    'Si':'Si',
    'OSi':'O',
    'Oice':'O',
    'Hice':'H',
    'Na1p':'NA',
    'K1p':'K',
    'He':'HE',
    'Li':'LI',
    'Be':'BE',
    'Ne':'NE',
    'Al':'AL',
    'Ar':'AR',
    'Sc':'SC',
    'Ti':'TI',
    'V':'V',
    'Cr':'CR',
    'Mn':'MN',
    'Ni':'NI',
    'Ga':'GA',
    'Ge':'GE',
    'As':'AS',
    'Se':'SE',
    'Kr':'KR',
    'Rb':'RB',
    'Sr':'SR',
    'Y':'Y',
    'Zr':'ZR',
    'Nb':'NB',
    'Mo':'MO',
    'Tc':'TC',
    'Ru':'RU',
    'Rh':'RH',
    'Pd':'PD',
    'Ag':'AG',
    'Cd':'CD',
    'In':'IN',
    'Sn':'SN',
    'Sb':'SB',
    'Te':'TE',
    'Xe':'XE',
    'Cs':'CS',
    'Ba':'BA',
    'La':'LA',
    'Ce':'CE',
    'Pr':'PR',
    'Nd':'ND',
    'Pm':'PM',
    'Sm':'SM',
    'Eu':'EU',
    'Gd':'GD',
    'Tb':'TB',
    'Dy':'DY',
    'Ho':'HO',
    'Er':'ER',
    'Tm':'TM',
    'Yb':'YB',
    'Lu':'LU',
    'Hf':'HF',
    'Ta':'TA',
    'W':'W',
    'Re':'RE',
    'Os':'OS',
    'Ir':'IR',
    'Pt':'PT',
    'Au':'AU',
    'Hg':'HG',
    'Tl':'TL',
    'Pb':'PB',
    'Bi':'BI',
    'Po':'PO',
    'At':'AT',
    'Rn':'RN',
    'Fr':'FR',
    'Ra':'RA',
    'Ac':'AC',
    'Th':'TH',
    'Pa':'PA',
    'U':'U',
    'Np':'NP',
    'Pu':'PU',
    'Am':'AM',
    'Cm':'CM',
    'Bk':'BK',
    'Cf':'CF',
    'Es':'ES',
    'Fm':'FM',
    'Md':'MD',
    'No':'NO',
    'Lr':'LR',
    'SUCK':'Z',
    'REPL':'Z',
    'REPLS':'Z',
    'HREPS':'Z',
    'VIRT':'X',
    'MPct':'X',
    'MPnm':'X',
    'MPdp':'X',
    'MPtk':'X',
}

def fix_element_name(element):
    if len(element) == 1:
        return element
    return element[0] + element[1].lower()
element_to_atom_type = {}
for key, value in atom_type_to_element.items():
    value = fix_element_name(value)
    atom_type_to_element[key] = value
    element_to_atom_type[value] = key


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

_atom_record_format = (
    "{ATOM:6s}{atomi:5d} {atomn:^4}{idx:^1}{resn:3s} {chain:1}{resi:4d}{insert:1s}   "
    "{x:8s}{y:8s}{z:8s}{occ:6.2f}{b:6.2f}           {elem:3s}\n"
)
def format_atom(
        ATOM='ATOM',
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


def get_atom_order_from_variants(name3, variants, params_dict={}):
    if name3 in params_dict:
        atom_order = params_dict[name3]['atom_order']
        assert len(variants) == 0, f'Variants on ligands no supported {name3} {variants}'
        return atom_order

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

# Rosetta is my favorite program
# It reorders params file on the fly
def reorder_params_rosetta_order(atom_order, elements, bonds, name_to_4char, first_sidechain_atom):
    # Rosetta's internal atom ordering:
    #   1. Backbone heavy atoms (params-file order, i.e. before FIRST_SIDECHAIN_ATOM)
    #   2. Sidechain heavy atoms (params-file order, i.e. at/after FIRST_SIDECHAIN_ATOM)
    #   3. Hydrogens grouped by parent heavy atom (following the heavy atom order above),
    #      with each group ordered by BOND-declaration order in the params file.

    elem_map = dict(zip(atom_order, elements))
    is_H = {name: (elem == 'H') for name, elem in zip(atom_order, elements)}

    if first_sidechain_atom is None:
        fsa_idx = 0          # no FIRST_SIDECHAIN_ATOM: all heavy atoms are sidechain
    elif first_sidechain_atom == 'NONE':
        fsa_idx = len(atom_order)   # all heavy atoms are backbone
    else:
        fsa_4char = name_to_4char.get(first_sidechain_atom)
        fsa_idx = atom_order.index(fsa_4char) if fsa_4char in atom_order else 0

    backbone_heavy = []
    sidechain_heavy = []
    for i, name in enumerate(atom_order):
        if not is_H[name]:
            if i < fsa_idx:
                backbone_heavy.append(name)
            else:
                sidechain_heavy.append(name)

    heavy_order = backbone_heavy + sidechain_heavy
    heavy_set = set(heavy_order)

    h_by_parent = {a: [] for a in heavy_order}
    h_seen = set()
    for a1_str, a2_str in bonds:
        a1 = name_to_4char.get(a1_str)
        a2 = name_to_4char.get(a2_str)
        if a1 is None or a2 is None:
            continue
        for heavy, h in [(a1, a2), (a2, a1)]:
            if heavy in heavy_set and is_H.get(h, False) and h not in h_seen:
                h_by_parent[heavy].append(h)
                h_seen.add(h)

    new_order = []
    new_elements = []
    for a in heavy_order:
        new_order.append(a)
        new_elements.append(elem_map[a])
    for a in heavy_order:
        for h in h_by_parent[a]:
            new_order.append(h)
            new_elements.append(elem_map[h])
    # Safety net: include any H atoms with no BOND entry
    for name in atom_order:
        if is_H[name] and name not in h_seen:
            new_order.append(name)
            new_elements.append(elem_map[name])

    return new_order, new_elements


def parse_raw_params(params_files):
    params_dict = {}
    for params_file in params_files:
        atom_order = []
        elements = []
        bonds = []
        name3 = None
        name_to_4char = {}
        first_sidechain_atom = None
        for line in params_file.split('\n'):
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith('IO_STRING') and len(line) >= 13:
                name3 = line[10:13]
            elif line.startswith('ATOM') and len(line) >= 9:
                atom_name_4 = line[5:9]
                atom_order.append(atom_name_4)
                elements.append(atom_type_to_element[line.split()[2]])
                name_to_4char[atom_name_4.strip()] = atom_name_4
            elif line.startswith('BOND'):
                sp = line.split()
                if len(sp) >= 3:
                    bonds.append((sp[1], sp[2]))
            elif line.startswith('FIRST_SIDECHAIN_ATOM'):
                sp = line.split()
                if len(sp) >= 2:
                    first_sidechain_atom = sp[1]
        if name3 is not None:
            atom_order, elements = reorder_params_rosetta_order(atom_order, elements, bonds, name_to_4char, first_sidechain_atom)
            atom_index = {name: i for i, name in enumerate(atom_order)}
            conect = [[] for _ in atom_order]
            for a1_str, a2_str in bonds:
                a1 = name_to_4char.get(a1_str)
                a2 = name_to_4char.get(a2_str)
                if a1 is not None and a2 is not None:
                    if a1 in atom_index:
                        conect[atom_index[a1]].append(a2)
                    if a2 in atom_index:
                        conect[atom_index[a2]].append(a1)
            d = {}
            d['atom_order'] = atom_order
            d['conect'] = conect
            d['elements'] = elements

            params_dict[name3] = d
    return params_dict



def structure_to_pdb(structure, renumber_resno_increasing=False, raw_params_files=[]):

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

    for infos in pdb_infos.values():
        for info in infos:
            if info.startswith("PARAMS:"):
                raw_params_files.append(base_64_decode_params(info.split(':')[-1]))

    params_dict = parse_raw_params(raw_params_files)

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
    conect_lines = []

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
        name3 = name1_to_name3.get(name1, None)


        variants = []
        if ( len(annotated_name) > 1 ):
            var_split = annotated_name[2:-1].split(":")
            variants = var_split[1:]

            # Because the silent file uses a different naming convention for DNA/RNA...
            if var_split[0] in silent_name3_to_pdb_name3:
                name3 = silent_name3_to_pdb_name3[var_split[0]]
            else:
                name3 = var_split[0]
            name3 = ' '*(3 - len(name3)) + name3

        assert name3 is not None, f"Couldn't figure out name3 for {name1}"

        atom_order = get_atom_order_from_variants(name3, variants, params_dict)

        atom_coords = silent_line_to_atoms(atom_lines[seqpos0][1:])

        assert(len(atom_coords) == len(atom_order))


        atom_name_to_atomno = {}
        ATOM = 'HETATM' if name3 in params_dict else 'ATOM'
        for iatom in range(len(atom_order)):
            atom_name = atom_order[iatom]
            if "V" in atom_name:
                continue
            atomno += 1
            atom_name_to_atomno[atom_name] = atomno

            if name3 in params_dict and params_dict[name3]['elements'][iatom]:
                element = params_dict[name3]['elements'][iatom]
            else:
                element = atom_name[1]
                if atom_name[0] not in ' 0123456789':
                    element = atom_name[0]

            pdb_lines.append(format_atom(
                ATOM=ATOM,
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

        if name3 in params_dict:
            for iatom, bonded_names in enumerate(params_dict[name3]['conect']):
                src_name = atom_order[iatom]
                if src_name not in atom_name_to_atomno:
                    continue
                src_no = atom_name_to_atomno[src_name]
                bonded_nos = [atom_name_to_atomno[n] for n in bonded_names if n in atom_name_to_atomno]
                for chunk_start in range(0, len(bonded_nos), 4):
                    chunk = bonded_nos[chunk_start:chunk_start+4]
                    conect_lines.append(f"CONECT{src_no:5d}" + "".join(f"{n:5d}" for n in chunk) + "\n")

        if seqpos0 + 1 in chain_endings:
            atomno += 1
            pdb_lines.append("TER                                                                             \n")


    if pdb_lines[-1] != "TER                                                                             \n":
        pdb_lines.append("TER                                                                             \n")

    for seqpos in sorted(pdb_infos):
        pdb_lines.append(f'REMARK PDBinfo-LABEL: {seqpos:4d} ' + " ".join(pdb_infos[seqpos]) + "\n")

    pdb_lines += conect_lines

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
def get_jump_stub(atom_dict, is_pro=False, atom_order=None):
    if atom_order is not None:
        a0, a1, a2 = atom_order[0], atom_order[1], atom_order[2]
        return get_stub_from_n_ca_c(np.array(atom_dict[a0]),
                                    np.array(atom_dict[a1]),
                                    np.array(atom_dict[a2]))

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


def build_conect_dict(raw_conect, atomno_to_info, hetatm_idents, atomno_to_xyz=None):
    # Build ident -> [atom_name_4, ...] in serial order
    ident_to_atom_order = {}
    for serial in sorted(atomno_to_info.keys()):
        atom_name, ident = atomno_to_info[serial]
        if ident not in ident_to_atom_order:
            ident_to_atom_order[ident] = []
        ident_to_atom_order[ident].append(atom_name)

    # Build bonds from raw_conect: ident -> {atom_name -> set of bonded atom_names}
    # Only track intra-residue bonds for HETATM idents
    bonds = {}
    for serials in raw_conect:
        src_serial = serials[0]
        if src_serial not in atomno_to_info:
            continue
        src_name, src_ident = atomno_to_info[src_serial]
        if src_ident not in hetatm_idents:
            continue
        if src_ident not in bonds:
            bonds[src_ident] = {}
        if src_name not in bonds[src_ident]:
            bonds[src_ident][src_name] = set()
        for s in serials[1:]:
            if s in atomno_to_info:
                dst_name, dst_ident = atomno_to_info[s]
                if dst_ident == src_ident:
                    bonds[src_ident][src_name].add(dst_name)

    # Make bidirectional: if A->B exists, ensure B->A exists
    for ident in list(bonds.keys()):
        for src, bonded_set in list(bonds[ident].items()):
            for dst in list(bonded_set):
                if dst not in bonds[ident]:
                    bonds[ident][dst] = set()
                bonds[ident][dst].add(src)

    # Invent bonds for hetatm atoms that have no connections using distance cutoffs
    if atomno_to_xyz is not None:
        # Build ident -> list of (serial, atom_name, xyz) for hetatm idents
        ident_to_serials = {}
        for serial in sorted(atomno_to_info.keys()):
            atom_name, ident = atomno_to_info[serial]
            if ident not in hetatm_idents:
                continue
            if serial not in atomno_to_xyz:
                continue
            if ident not in ident_to_serials:
                ident_to_serials[ident] = []
            ident_to_serials[ident].append((serial, atom_name, atomno_to_xyz[serial]))

        for ident, atoms in ident_to_serials.items():
            if len(atoms) < 2:
                continue

            # Find atoms with no bonds
            unconnected = []
            for serial, atom_name, xyz in atoms:
                has_bond = (ident in bonds and
                            atom_name in bonds[ident] and
                            len(bonds[ident][atom_name]) > 0)
                if not has_bond:
                    unconnected.append((serial, atom_name, xyz))

            if not unconnected:
                continue

            if ident not in bonds:
                bonds[ident] = {}

            for serial, atom_name, xyz in unconnected:
                # Determine hydrogen vs heavy atom by first non-space char of atom_name_4
                first_non_space = atom_name.lstrip(' ')
                is_hydrogen = len(first_non_space) > 0 and first_non_space[0] == 'H'
                cutoff = 1.3 if is_hydrogen else 2.0

                xyz_arr = np.array(xyz)

                # Compute distances to all other atoms in same ident
                candidates = [(other_name, np.linalg.norm(xyz_arr - np.array(other_xyz)))
                              for other_serial, other_name, other_xyz in atoms
                              if other_serial != serial]

                bonded = [name for name, dist in candidates if dist <= cutoff]

                if not bonded:
                    # No atoms within cutoff: pick the closest
                    bonded = [min(candidates, key=lambda x: x[1])[0]]

                if atom_name not in bonds[ident]:
                    bonds[ident][atom_name] = set()
                for b in bonded:
                    bonds[ident][atom_name].add(b)
                    # Keep bidirectional
                    if b not in bonds[ident]:
                        bonds[ident][b] = set()
                    bonds[ident][b].add(atom_name)

    # Build final conect_dict in serial order for all HETATM idents
    conect_dict = {}
    for ident in bonds:
        if ident not in ident_to_atom_order:
            continue
        conect_dict[ident] = []
        for atom_name in ident_to_atom_order[ident]:
            if atom_name in bonds[ident]:
                conect_dict[ident].append([atom_name] + list(bonds[ident][atom_name]))

    # Catch single-atom ligands like ZN and make sure they end up in conect_dict
    for ident in hetatm_idents:
        if ident not in conect_dict and ident in ident_to_atom_order:
            first_atom = ident_to_atom_order[ident][0]
            conect_dict[ident] = [[first_atom]]

    return conect_dict


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

    atomno_to_info = {}  # serial int -> (atom_name_4, ident)
    atomno_to_xyz = {}   # serial int -> [x, y, z]
    raw_conect = []      # list of lists of int serials
    hetatm_idents = set()
    hetatm_elements = {}  # ident -> {atom_name_4 -> element}
    hetatm_coords = {}   # ident -> {atom_name_4 -> [x, y, z]}

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

        if line.startswith('CONECT'):
            serials = []
            for col_start in [6, 11, 16, 21, 26]:
                s = line[col_start:col_start+5].strip() if len(line) >= col_start+5 else ''
                if s:
                    serials.append(int(s))
            if serials:
                raw_conect.append(serials)
            continue

        if not (line.startswith('ATOM') or line.startswith('HETATM')) or len(line) < 54:
            if len(aa_atom_dicts) > 0:
                if len(chain_ends) == 0 or chain_ends[-1] != len(aa_atom_dicts):
                    chain_ends.append(len(aa_atom_dicts))
            last_ident = ''

            continue

        ident = line[17:27]
        serial = int(line[6:11].strip())
        atom = line[12:16]

        # Stupid biotite
        if atom[0] == ' ' and atom[1].isnumeric():
            atom = atom[1:] + ' '

        if line.startswith('HETATM'):
            hetatm_idents.add(ident)
            element = line[76:78].strip() if len(line) >= 78 else ''
            hetatm_elements.setdefault(ident, {})[atom] = element

        atomno_to_info[serial] = (atom, ident)

        x = float(line[30:30+8])
        y = float(line[38:38+8])
        z = float(line[46:46+8])
        atomno_to_xyz[serial] = [x, y, z]

        if line.startswith('HETATM'):
            hetatm_coords.setdefault(ident, {})[atom] = [x, y, z]

        name3 = line[17:20]
        chain = line[21]
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

    conect_dict = build_conect_dict(raw_conect, atomno_to_info, hetatm_idents, atomno_to_xyz)

    return chain_ends, aa_atom_dicts, name3s, chain_letters, res_nums, in_a_disulfides, disulfides, pdb_info_labels, conect_dict, hetatm_elements, hetatm_coords


def fmt(value):
    if isinstance(value, str):
        return value
    else:
        return '%.3f'%value


def _icoor_values(child_xyz, parent_xyz, angle_xyz, torsion_xyz):
    """Compute Rosetta ICOOR_INTERNAL (phi_deg, theta_deg, d) from Cartesian coords.

    Mirrors Stub::from_four_points(center=parent, a=parent, b=angle, c=torsion):
      e1 (x) = normalized(parent - angle)
      e3 (z) = normalized(cross(e1, torsion - angle))
      e2 (y) = cross(e3, e1)
    Then inverts Stub::spherical(phi, theta, d).
    """
    c  = np.array(child_xyz,   dtype=float)
    p  = np.array(parent_xyz,  dtype=float)
    a  = np.array(angle_xyz,   dtype=float)
    t  = np.array(torsion_xyz, dtype=float)

    diff = c - p
    d = float(np.linalg.norm(diff))
    if d < 1e-9:
        return 0.0, 0.0, 0.0

    e1 = p - a
    e1_norm = np.linalg.norm(e1)
    if e1_norm < 1e-9:
        return 0.0, 0.0, d
    e1 = e1 / e1_norm

    cb = t - a  # c - b in from_four_points notation
    e3 = np.cross(e1, cb)
    e3_norm = np.linalg.norm(e3)
    if e3_norm < 1e-9:  # collinear — nudge as Rosetta does
        cb = cb + np.array([1.4e-7, 6.7e-8, 2.3e-7])
        e3 = np.cross(e1, cb)
        e3_norm = np.linalg.norm(e3)
    e3 = e3 / e3_norm

    e2 = np.cross(e3, e1)

    x_comp = float(np.dot(diff, e1))
    y_comp = float(np.dot(diff, e2))
    z_comp = float(np.dot(diff, e3))

    theta = float(np.degrees(np.arccos(np.clip(x_comp / d, -1.0, 1.0))))
    phi   = float(np.degrees(np.arctan2(z_comp, y_comp)))
    return phi, theta, d


def generate_raw_params(name3, d):
    atom_order = d['atom_order']
    conect = d['conect']
    elements = d['elements']
    n = len(atom_order)
    stripped = [a.strip() for a in atom_order]

    lines = []
    lines.append(f'NAME {name3.strip()}')
    lines.append(f'IO_STRING {name3} X')
    lines.append('TYPE LIGAND')
    lines.append('AA UNK')

    for atom_name, element in zip(atom_order, elements):
        if atom_name in (' V1 ', ' V2 '):
            lines.append(f'ATOM {atom_name} VIRT VIRT 0.00')
        else:
            atom_type = element_to_atom_type[fix_element_name(element)]
            lines.append(f'ATOM {atom_name} {atom_type:4s}  CT1  0.00')

    written_bonds = set()
    for i, bonded_names in enumerate(conect):
        src = stripped[i]
        for dst_4 in bonded_names:
            dst = dst_4.strip()
            bond = tuple(sorted([src, dst]))
            if bond not in written_bonds:
                written_bonds.add(bond)
                lines.append(f'BOND  {src}  {dst}')

    lines.append(f'NBR_ATOM {stripped[0]}')
    lines.append('NBR_RADIUS 3.0')

    if n == 1:
        a0 = stripped[0]
        lines.append(f'ICOOR_INTERNAL  {a0}   0   0   0  {a0}  {a0}  {a0}')
    elif n == 2:
        a0, a1 = stripped[0], stripped[1]
        lines.append(f'ICOOR_INTERNAL  {a0}   0   0   0  {a0}  {a1}  {a1}')
        lines.append(f'ICOOR_INTERNAL  {a1}   1   1   1  {a0}  {a1}  {a1}')
    else:
        # Build bidirectional adjacency list using stripped atom names
        stripped_set = set(stripped)
        adj = defaultdict(list)
        for i, bonded_list in enumerate(conect):
            src = stripped[i]
            for dst_4 in bonded_list:
                dst = dst_4.strip()
                if dst not in stripped_set:
                    continue
                if dst not in adj[src]:
                    adj[src].append(dst)
                if src not in adj[dst]:
                    adj[dst].append(src)

        # Choose initial chain a0->a1->a2 along bonded atoms.
        # a2 is preferably bonded to a1 (extending the chain), or falls back to
        # a second neighbor of a0 if a1 is a leaf.
        a0 = stripped[0]
        a1 = adj[a0][0] if adj[a0] else stripped[1]
        a2_via_a1 = [x for x in adj[a1] if x != a0]
        if a2_via_a1:
            a2, a2_parent = a2_via_a1[0], a1
        else:
            a2_via_a0 = [x for x in adj[a0] if x != a1]
            a2, a2_parent = (a2_via_a0[0], a0) if a2_via_a0 else (stripped[2], a1)

        # Build BFS spanning tree, honoring a1 as first child of a0 and a2
        # under a2_parent, so all parent references are bonded.
        parent_map = {a0: None, a1: a0, a2: a2_parent}
        seen = {a0, a1, a2}
        bfs = [a0, a1, a2]
        bi = 0
        while bi < len(bfs):
            atom = bfs[bi]; bi += 1
            for nbr in adj[atom]:
                if nbr not in seen:
                    seen.add(nbr)
                    parent_map[nbr] = atom
                    bfs.append(nbr)

        # Build stripped-key coords lookup for ICOOR computation
        raw_coords = d.get('coords', {})
        coords_stripped = {k.strip(): v for k, v in raw_coords.items()}

        def icoor_line(child, parent, angle, torsion, i):
            if coords_stripped and all(x in coords_stripped for x in [child, parent, angle, torsion]):
                phi, theta, dist = _icoor_values(
                    coords_stripped[child], coords_stripped[parent],
                    coords_stripped[angle], coords_stripped[torsion])
                return f'ICOOR_INTERNAL  {child}   {phi:.6f}   {theta:.6f}   {dist:.6f}  {parent}  {angle}  {torsion}'
            return f'ICOOR_INTERNAL  {child}   {i}.000000   {i}.000000   {i}.000000  {parent}  {angle}  {torsion}'

        # Write the three stub ICOOR lines (special Rosetta convention).
        # Line 3: parent is a2_parent; angle is the other first-chain atom.
        ang3 = a0 if a2_parent == a1 else a1
        lines.append(icoor_line(a0, a0, a1, a2, 0))
        lines.append(icoor_line(a1, a0, a1, a2, 1))
        lines.append(icoor_line(a2, a2_parent, ang3, a2, 2))

        written = {a0, a1, a2}

        # Write remaining atoms in BFS order.  For each atom we need:
        #   parent  – tree parent (bonded to atom by construction)
        #   angle   – tree grandparent, or any placed bonded neighbor of parent
        #   torsion – tree great-grandparent, or any placed bonded neighbor of angle
        for atom in bfs:
            if atom in written:
                continue
            p = parent_map[atom]

            # angle: grandparent in tree, else any placed neighbor of p (not atom)
            g = parent_map.get(p)
            if g is None or g not in written:
                g = next((x for x in adj[p] if x != atom and x in written), p)

            # torsion: great-grandparent in tree, else any placed neighbor of g
            # (excluding p and atom; also reject t==p which makes the dihedral degenerate)
            t = parent_map.get(g)
            if t is None or t not in written or t == atom or t == p:
                t = next((x for x in adj[g] if x != p and x in written and x != atom), None)
                if t is None:
                    t = next((x for x in adj[g] if x in written and x != atom), p)

            lines.append(icoor_line(atom, p, g, t, len(written)))
            written.add(atom)

    return '\n'.join(lines) + '\n'


def base_64_encode_params(params_string):
    return base64.b64encode(params_string.encode()).decode()


def base_64_decode_params(encoded_string):
    return base64.b64decode(encoded_string.encode()).decode()


def conect_dict_to_params_dict(conect_dict, hetatm_elements=None, hetatm_coords=None):
    params_dict = {}
    for ident, conect_records in conect_dict.items():
        name3 = ident[0:3]

        atom_order = []
        seen = set()
        bonded_map = {}
        for record in conect_records:
            src = record[0]
            if src not in seen:
                atom_order.append(src)
                seen.add(src)
                bonded_map[src] = []
            bonded_map[src].extend(record[1:])

        conect = [bonded_map[name] for name in atom_order]

        if len(atom_order) < 3:
            first_atom = atom_order[0]
            conect[0].append(' V1 ')
            atom_order.append(' V1 ')
            conect.append([first_atom, ' V2 '])
            atom_order.append(' V2 ')
            conect.append([' V1 '])

        elem_map = hetatm_elements.get(ident, {}) if hetatm_elements is not None else {}
        element = [elem_map.get(name, '') for name in atom_order]

        # Build coords dict, adding synthetic coords for virtual atoms
        ident_coords = hetatm_coords.get(ident, {}) if hetatm_coords is not None else {}
        coords = dict(ident_coords)
        if coords:
            first_xyz = coords.get(atom_order[0])
            if first_xyz is not None:
                if ' V1 ' in atom_order and ' V1 ' not in coords:
                    coords[' V1 '] = [first_xyz[0] + 1.0, first_xyz[1], first_xyz[2]]
                if ' V2 ' in atom_order and ' V2 ' not in coords:
                    v1_xyz = coords.get(' V1 ', first_xyz)
                    coords[' V2 '] = [v1_xyz[0], v1_xyz[1] + 1.0, v1_xyz[2]]

        d = {'atom_order': atom_order, 'conect': conect, 'elements': element, 'coords': coords}

        if name3 in params_dict:
            existing = params_dict[name3]
            assert existing['atom_order'] == atom_order and existing['conect'] == conect, \
                f'Conflicting params for {name3.strip()!r}: atom_order or connectivity differs between residues'
        else:
            params_dict[name3] = d

    return params_dict


def add_virtual_coords(d, atom_order):
    d = dict(d)
    first_xyz = d[atom_order[0]]
    if ' V1 ' in atom_order and ' V1 ' not in d:
        d[' V1 '] = [first_xyz[0] + 1.0, first_xyz[1], first_xyz[2]]
    if ' V2 ' in atom_order and ' V2 ' not in d:
        v1_xyz = d[' V1 ']
        d[' V2 '] = [v1_xyz[0], v1_xyz[1] + 1.0, v1_xyz[2]]
    return d


def parsed_pdb_to_silent( chain_ends, aa_atom_dicts, name3s, chain_letters, res_nums, in_a_disulfides, disulfides, tag,
                                                                            write_header=False, score_dict={}, pdb_info_labels={},
                                                                            params_dict={},
                                                                            allow_missing_atoms=False,
                                                                            allow_extra_atoms=False,
                                                                            full_params_files={}):
    assert len(aa_atom_dicts) > 0, 'Error! There are no atoms in your pdb?'

    # Add virtual atoms to metal ions
    for i, name3 in enumerate(name3s):
        if name3 in params_dict:
            atom_order = params_dict[name3]['atom_order']
            if ' V1 ' in atom_order or ' V2 ' in atom_order:
                aa_atom_dicts[i] = add_virtual_coords(aa_atom_dicts[i], atom_order)


    local_name3_to_name1 = {k:v for k, v in name3_to_name1.items()}
    local_needs_annotated_name = {k for k in needs_annotated_name}


    pdb_info_labels = {k:v for k,v in pdb_info_labels.items()}

    # Store params into first pdb info label and update our name3 vectors
    for params_name3, d in params_dict.items():
        if params_name3 in full_params_files:
            raw_params = full_params_files[params_name3]
        else:
            raw_params = generate_raw_params(params_name3, d)
        params_header = f'PARAMS:{params_name3.replace(" ", "")}:'
        duplicate = False
        for prev_info in pdb_info_labels.setdefault(1, []):
            if prev_info.startswith(params_header):
                duplicate = True
        if not duplicate:
            pdb_info_labels[1].append(params_header + base_64_encode_params(raw_params))
        local_name3_to_name1[params_name3] = 'X'
        local_needs_annotated_name.add(params_name3)


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
    stub1 = get_jump_stub(aa_atom_dicts[0], is_pro=name3s[0] == 'PRO',
                          atom_order=params_dict.get(name3s[0], {}).get('atom_order'))
    for end in chain_ends[:-1]:
        stub2 = get_jump_stub(aa_atom_dicts[end], is_pro=name3s[end] == 'PRO',
                              atom_order=params_dict.get(name3s[end], {}).get('atom_order'))
        jump_rt = np.linalg.inv(stub1) @ stub2

        post_silent += f'RT {" ".join("%.8f"%x for x in list(jump_rt[:3,:3].flat) + list(jump_rt[:3,3].flat))} {tag}\n'

    # residue lines
    missing_atom = [0, 0, 0]
    annotated_seq = ''

    for ires, (d, name3, in_a_disulfide) in enumerate(zip(aa_atom_dicts, name3s, in_a_disulfides)):

        assert name3 in rosetta_aa_order or name3 in params_dict, 'Error! Unknown name3: ' + name3
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
        elif name3 == 'CYS' and ' HG ' not in d:
            variants.append('disulfide')

        # there's no way to detect LowerDNA, but it isn't necessary
        if "HO3'" in d:
            variants.append('UpperDNA')

        # RNA doesn't even change anything with LowerDNA and UpperDNA

        atom_order = get_atom_order_from_variants(name3, variants, params_dict)

        atom_used = {}
        for atom_name in d:
            atom_used[atom_name] = False

        messages = []
        no_atom_errors = True
        res_coords = []
        renamed_h_used = set()
        for atom_name in atom_order:
            og_atom_name = atom_name

            not_found = atom_name not in d
            if not_found and atom_name in synonyms:
                atom_name = synonyms[atom_name]
                not_found = atom_name not in d

            # We're missing a hydrogen, check for biotite naming scheme
            if not_found and atom_name[0].isnumeric() and atom_name[1] == 'H':
                biotite_format = atom_name[1:].strip() + '%i'
                if len(biotite_format) < 5:
                    biotite_format = ' ' + biotite_format
                biotite_format += ' '*(5 - len(biotite_format))
                for i_h in range(1, 4):
                    try_name = biotite_format%i_h
                    if try_name in renamed_h_used:
                        continue
                    if try_name in d:
                        atom_name = try_name
                        renamed_h_used.add(try_name)
                        not_found = False
                        break


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

        annotated_seq += local_name3_to_name1[name3]
        if name3 in local_needs_annotated_name or len(variants) > 0:

            write_name3 = name3.strip()
            if name3 in pdb_name3_to_silent_name3:
                write_name3 = pdb_name3_to_silent_name3[name3]
            to_add = '[' + ':'.join( [write_name3] + variants ) + ']'
            annotated_seq += to_add


    silent += "ANNOTATED_SEQUENCE: %s %s\n"%(annotated_seq, tag)
    silent += post_silent

    return silent

# score_dict is string -> float or string
# pdb_info_labels is int -> list(string)
def pdb_to_structure(pdb_lines, tag, write_header=True, score_dict={}, pdb_info_labels={}, allow_extra_atoms=False, full_params_files={}):

    chain_ends, aa_atom_dicts, name3s, chain_letters, res_nums, in_a_disulfides, disulfides, more_info_labels, conect_dict, hetatm_elements, hetatm_coords = parse_pdb_into_needed_format(pdb_lines)

    params_dict = conect_dict_to_params_dict(conect_dict, hetatm_elements, hetatm_coords)

    for key in full_params_files:
        params_dict[key] = parse_raw_params([full_params_files[key]])[key]

    # if you add into the pdb_info_labels object you modify global state
    for key in pdb_info_labels:
        if key in more_info_labels:
            more_info_labels[key] += pdb_info_labels[key]
        else:
            more_info_labels[key] = pdb_info_labels[key]


    silent = parsed_pdb_to_silent( chain_ends, aa_atom_dicts, name3s, chain_letters, res_nums, in_a_disulfides, disulfides, tag,
                                                write_header=write_header, score_dict=score_dict, pdb_info_labels=more_info_labels,
                                                params_dict=params_dict, allow_extra_atoms=allow_extra_atoms, full_params_files=full_params_files)

    return silent



def generate_params_for_pdb(pdb_lines):

    chain_ends, aa_atom_dicts, name3s, chain_letters, res_nums, in_a_disulfides, disulfides, more_info_labels, conect_dict, hetatm_elements, hetatm_coords = parse_pdb_into_needed_format(pdb_lines)

    params_dict = conect_dict_to_params_dict(conect_dict, hetatm_elements, hetatm_coords)

    raw_params_dict = {}
    for name3, d in params_dict.items():
        raw_params_dict[name3] = generate_raw_params(name3, d)

    return raw_params_dict


def strip_to_ligands(pose):
    from pyrosetta.rosetta.core.select.residue_selector import ResiduePropertySelector
    from pyrosetta.rosetta.core.chemical import PROTEIN
    # 1. Create a selector for all PROTEIN residues
    protein_selector = ResiduePropertySelector(PROTEIN)
    
    # 2. Get a boolean mask (vector1_bool) of protein residues
    protein_mask = protein_selector.apply(pose)
    
    # 3. Convert mask to a list of indices (in reverse to avoid index shifting)
    # Rosetta indices are 1-based
    protein_indices = [i for i, is_protein in enumerate(protein_mask, 1) if is_protein]
    
    if len(protein_indices) == pose.size():
        return None

    # 4. Delete the protein residues
    # Deleting in reverse order is safer when using delete_residue()
    for res_idx in sorted(protein_indices, reverse=True):
        pose.delete_residue_slow(res_idx)
    
    return pose


def inject_params_into_pyrosetta_pose(in_pose):
    from pyrosetta.rosetta.std import stringstream

    pose = strip_to_ligands(in_pose.clone())
    if pose is None:
        return

    ss = stringstream()
    pose.dump_pdb(ss)
    pdb_lines = ss.str().split("\n")

    raw_params_dict = generate_param_for_pdb(pdb_lines)

    prev_labels = list(in_pose.pdb_info().get_reslabels(1))
    prev_ligands = set([item.split(':')[1] for item in prev_labels if item.startswith("PARAMS:")])
    for name3, raw_params in raw_params_dict.items():
        if name3 in prev_ligands:
            continue
        params_header = f'PARAMS:{name3.replace(" ", "")}:'
        in_pose.pdb_info().add_reslabel(1, params_header + base_64_encode_params(raw_params))



        

