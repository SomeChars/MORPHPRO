from .basic_operations import *
from .statistics import *
import time
import subprocess
from .pdb_utils import *
from .broken_line import morph as bl_m
from ..models import MorphRequest, Pdb


def do_test():
    start = [
        [3.8, 0, 0],
        [3.8 * 2, 0, 0],
        [3.8 * 3, 0, 0],
        [3.8 * 4, 0, 0],
        [3.8 * 6, 0, 0],
        [3.8 * 7, 0, 0],
        [3.8 * 8, 0, 0],
        [3.8 * 9, 0, 0]
    ]

    finish = [
        [3.8 * 3, 1, 0],
        [3.8 * 2, 0, 0],
        [3.8 * 3, 0, 0],
        [3.8 * 4, 0, 0],
        [3.8 * 6, 0, 0],
        [3.8 * 7, 0, 0],
        [3.8 * 8, 0, 0],
        [3.8 * 7, 0, 1]
    ]
    idx_start = [0, 1, 2, 3, 5, 6, 7, 8]
    idx_finish = [0, 1, 2, 3, 5, 6, 7, 8]
    return start, finish, idx_start, idx_finish


def naive_morph(morph, test=False):
    if test:
        start, finish, idx_start, idx_finish = do_test()
    else:
        start, finish, idx_start, idx_finish = __read_lines_and_align(morph)

    if len(start) != len(finish):
        print(len(start), ' is len of start')
        return None
    start, finish = qcp_align(start, finish)
    # what are those?
    return bl_m(start, finish, morph.morphing_count)


def analyze_morph(morph):
    start, finish = __read_lines(morph)
    if len(start) == len(finish):
        start, finish = qcp_align(start, finish)
    analyze(start, finish)


def castellana_pevzner_oe_morph(morph, test=False):
    return castellana_pevzner_morph(morph, "oe", test)


def castellana_pevzner_oea_morph(morph, test=False):
    castellana_pevzner_morph(morph, "oea", test)


def castellana_pevzner_oeac_morph(morph, test=False):
    castellana_pevzner_morph(morph, "oeac", test)


def __max_dist(start, finish):
    d = dist(start[0], finish[0])
    for i in range(len(start)):
        x, y = start[i], finish[i]
        d = max(d, dist(x, y))
    return d

# если вы чините этот сервер - простите меня за эту функцию, мне жаль
def __answer_parser(line):
    a = line[0].decode()
    a = a.split('\n')[1:-1]
    b = []
    c = []
    for i in a:
        for j in i.split():
            b.append(float(j))
        c.append(b)
        b = []

    n = len(c[0])
    d = []
    e = []
    for i in c:
        for j in range(n // 3):
            e.append([i[j], i[j + 1], i[j + 2]])
        d.append(e)
        e = []
    return d


def castellana_pevzner_morph(morph_request, filename, test=False):
    if test:
        start, finish, idx_start, idx_finish = do_test()
    else:
        start, finish, idx_start, idx_finish = __read_lines_and_align(morph_request)

    if len(start) != len(finish):
        return None

    start, finish = qcp_align(start, finish)
    n = len(start)
    t1 = time.time()

    script = subprocess.Popen('morphserverapp/algorithms/extern/castellana-pevzner/' + filename, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    script_input = str(n) + '\n'
    for i in start:
        script_input += '{} {} {}\n'.format(i[0], i[1], i[2])
    for i in idx_start:
        script_input += '{}\n'.format(i)
    for i in finish:
        script_input += '{} {} {}\n'.format(i[0], i[1], i[2])
    for i in idx_finish:
        script_input += '{}\n'.format(i)
    if morph_request.auto_interpolation:
        interpolation = '{}\n'.format(3+int(__max_dist(start, finish)))
    else:
        interpolation = '{}\n'.format(morph_request.morphing_count)
    script_input += interpolation
    script_input = script_input.encode()

    line = script.communicate(input=script_input)
    t2 = time.time()

    parsed_line = __answer_parser(line)
    return parsed_line


def qcp_align(start, finish, test=False):
    if test:
        start, finish = do_test()[:2]

    n = len(start)
    script = subprocess.Popen('morphserverapp/algorithms/extern/qcp/run', stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    script_input = str(n) + '\n'
    for i in start:
        script_input += '{} {} {}\n'.format(i[0], i[1], i[2])
    for i in finish:
        script_input += '{} {} {}\n'.format(i[0], i[1], i[2])
    script_input = script_input.encode()
    line = script.communicate(input=script_input)
    # парсим ответ
    tokens = line[0]
    tokens = tokens.decode()
    tokens = tokens.split()
    rms = tokens[0]
    # TODO
    #  Rails.logger.info("Theoretical RMSD: #{rms}")
    #  translate this to python

    m = [tokens[3 * i + 1:3 * i + 4] for i in range(3)]
    for i in range(3):
        for j in range(3):
            m[i][j] = float(m[i][j])

    # центрируем начало и конец
    sc = centroid(start)
    fc = centroid(finish)
    s, f = [], []
    for i in range(n):
        s += [minus(start[i], sc)]
        f += [minus(finish[i], fc)]

    for i in range(1):
        f[i] = multiply(m, f[i])
    return s, f


def __read_lines(mr):
    start = to_broken_line(Pdb.objects.get(name=mr.protein_a_name).file.path)
    finish = to_broken_line(Pdb.objects.get(name=mr.protein_b_name).file.path)
    return start, finish


def __read_lines_and_align(mr):
    return global_align(Pdb.objects.get(name=mr.protein_a_name).file.path, Pdb.objects.get(name=mr.protein_b_name).file.path)


if __name__ == '__main__':
    # castellana_pevzner_morph(0, 'oe', True)
    ans = qcp_align(0, 'oe', test=True)
    print(ans)
    pass
