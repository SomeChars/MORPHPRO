import subprocess
from hashlib import md5
from Bio import PDB as pdb
from . import amino_acids_dict

directory = '../public/data/'


# saves pdb to a file with a name equal to its md5 hash
# returns this hash as a result
def save_to_file(text):
    '''saves pdb to a file with a name equal to its md5 hash
    returns this hash as a result'''
    hasher = md5()
    hasher.update(text.encode())
    name = hasher.hexdigest()
    path = directory + name + '.pdb'
    try:
        with open(path, mode='r') as f:
            pass
            # file already exists
    except IOError:
        f = open(path, mode='w')
        f.write(text)
        f.flush()
        f.close()
    return name


id_struct = 0


def load_from_file(file_name):
    parser = pdb.PDBParser(QUIET=True)
    global id_struct
    if file_name[-4:] != '.pdb':
        name = file_name + '.pdb'
    else:
        name = file_name
    struct = parser.get_structure('id_' + str(id_struct), name)
    id_struct += 1
    return struct


def __to_broken_line_struct(struct):
    ans = []
    for residue in struct.get_chains():
        for atom in residue.get_atoms():
            if atom.name == 'CA':
                coords = atom.get_coord()
                ans.append([coords[0], coords[1], coords[2]])
    return ans


def to_broken_line(text):
    struct = load_from_file(text)
    return __to_broken_line_struct(struct)


def global_align(pdb1, pdb2):
    s1 = load_from_file(pdb1)
    s2 = load_from_file(pdb2)
    rs1 = [residue for residue in s1.get_residues()]
    rs2 = [residue for residue in s2.get_residues()]
    seq1, seq2 = '', ''
    for r in rs1:
        rname = r.get_resname()
        if rname.lower() in amino_acids_dict.aminos:
            seq1 += amino_acids_dict.c3to1(rname.lower())
    for r in rs2:
        rname = r.get_resname()
        if rname.lower() in amino_acids_dict.aminos:
            seq2 += amino_acids_dict.c3to1(rname.lower())

    n, m = len(seq1), len(seq2)
    if n == m:
        start = __to_broken_line_struct(s1)
        finish = __to_broken_line_struct(s2)
        idx1 = [i for i in range(n)]
        return (start, finish, idx1, idx1)

    script = subprocess.Popen(['morphserverapp/algorithms/extern/global-alignment/run', seq1, seq2], stdout=subprocess.PIPE)
    output = script.communicate()
    res, error = output
    if error:
        print("Error", error)
        return error
    res = res.decode()
    tokens = res.split('\n')
    # print("\n\nTOKENS:", tokens)
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
    return (bl1, bl2, idx1, idx2)


if __name__ == '__main__':
    ans = global_align('2J53', '6l4i')
    print(ans)
    pass
