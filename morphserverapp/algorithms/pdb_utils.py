import subprocess
#from hashlib import md5
from Bio import PDB as pdb
from . import amino_acids_dict


id_struct = 0
def load_from_file(filename):
    parser = pdb.PDBParser(QUIET=True)
    global id_struct
    if filename[-4:] != '.pdb':
        name = filename + '.pdb'
    else:
        name = filename
    struct = parser.get_structure('id_'+str(id_struct), name)
    id_struct += 1
    return struct


def __to_broken_line_struct(struct):
    ans = []
    for residue in struct.get_chains():
        for atom in residue.get_atoms():
            if atom.name == 'CA':
                coords = atom.get_coord()
                ans.append((coords[0], coords[1], coords[2]))
    return ans

def to_broken_line(filename):
    struct = load_from_file(filename)
    return __to_broken_line_struct(struct)

def global_align(pdb1,pdb2):
    s1 = load_from_file(pdb1)
    s2 = load_from_file(pdb2)
    rs1 = [residue.get_resname() for residue in s1.get_residues()]
    rs2 = [residue.get_resname() for residue in s2.get_residues()]
    seq1, seq2 = '', ''
    for r in rs1:
        if r.lower() in amino_acids_dict.aminos:
            seq1 += amino_acids_dict.c3to1(r.lower())
    for r in rs2:
        if r.lower() in amino_acids_dict.aminos:
            seq2 += amino_acids_dict.c3to1(r.lower())

    n, m = len(seq1), len(seq2)
    if n == m:
        start = __to_broken_line_struct(s1)
        finish = __to_broken_line_struct(s2)
        idx1 = [i for i in range(n)]
        return start, finish, idx1, idx1

    script = subprocess.Popen(['morphserverapp/algorithms/extern/global-alignment/run',seq1,seq2], stdout=subprocess.PIPE)
    output = script.communicate()
    res, error = output
    if error:
        print("Error", error)
        return error
    res = res.decode()
    tokens = res.split('\r\n')
    nn = int(tokens[0])

    bl1, idx1 = [], []
    ids = tokens[1].split()
    for i in range(nn):
        residue = rs1[int(ids[i])]
        for atom in residue.get_atoms():
            if atom.name == 'CA':
                coords = atom.get_coord()
                bl1.append((coords[0], coords[1], coords[2]))
                idx1.append(ids[i])

    bl2, idx2 = [], []
    ids = tokens[2].split()
    for i in range(nn):
        residue = rs2[int(ids[i])]
        for atom in residue.get_atoms():
            if atom.name == 'CA':
                coords = atom.get_coord()
                bl2.append((coords[0], coords[1], coords[2]))
                idx2.append(ids[i])
    return bl1, bl2, idx1, idx2