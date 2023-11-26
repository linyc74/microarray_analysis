import os
import pandas as pd
from os.path import join
from typing import List


def get_temp_path(
        prefix: str = 'temp',
        suffix: str = '') -> str:

    i = 1
    while True:
        fpath = f'{prefix}{i:03}{suffix}'
        if not os.path.exists(fpath):
            return fpath
        i += 1


def get_files(
        source: str = '.',
        startswith: str = '',
        endswith: str = '',
        isfullpath: bool = False) -> List[str]:

    ret = []
    s, e = startswith, endswith
    for path, dirs, files in os.walk(source):
        if path == source:
            ret = [f for f in files if (f.startswith(s) and f.endswith(e))]

    if isfullpath:
        ret = [join(source, f) for f in ret]

    if ret:
        ret.sort()  # make the order consistent across OS platforms
    return ret


def left_join(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    n = len(left)
    merged = left.merge(
        right=right,
        right_index=True,
        left_index=True,
        how='left'
    )
    assert len(merged) == n
    return merged
