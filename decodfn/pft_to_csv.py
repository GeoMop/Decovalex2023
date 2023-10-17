import os.path
import sys
import pandas
import re
import os

def pass_none(fn, none_val=None):
    return lambda x: none_val if (x is None) else fn(x)
def parse_col(col_str_q):
    col_str = col_str_q.strip(' "\n')
    pattern="^(?P<var>.+?)( (?P<obs_name>[a-zA-Z_]+) \((?P<cell_id>[0-9]+)\) \((?P<point>[0-9 -.]+)\))?$"

    match = re.match(pattern, col_str)
    cell_id = pass_none(int, 0)(match.group('cell_id'))
    split_point = lambda s: [float(coord) for coord in s.split(' ')]
    point = pass_none(split_point, [])(match.group('point'))
    return match.group('var'), pass_none(str, '')(match.group('obs_name')), cell_id, point

    # items = col_str.split(' obs_slope ') + ['() (  )',]
    # var, cell_point = items[:2]
    # cell, point = cell_point.split(") (")
    # return var, cell.lstrip('('), point.rstrip(')').split(' ')

file_in = sys.argv[1]
with open(file_in, 'r') as fin:
    header = fin.readline()
    cols = [parse_col(c) for c in header.split(',')]
    col_names = [v+':'+o for v, o, c, p in cols]
    cell_ids = [c for v, o, c, p in cols]
    points = [str(p) for v, o, c, p in cols]
    df = pandas.read_csv(fin, delim_whitespace=True, names=col_names, index_col=False)

csv_file_name = os.path.splitext(file_in)[0] + '.csv'
with open(csv_file_name, 'w') as fout:
    col_ln = '; '.join(col_names)
    cell_ln = '; '.join([str(c) for c in cell_ids])
    point_ln = '; '.join([str(p) for p in points])
    lines = f"{col_ln}\n{cell_ln}\n{point_ln}\n"
    fout.write(lines)
    df.to_csv(fout, header=False, index=False, sep=';')