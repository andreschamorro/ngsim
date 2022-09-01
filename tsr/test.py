# Test for the extension.
#
# Eli Bendersky [http://eli.thegreenplace.net]
# This code is in the public domain.
from collections import defaultdict
import ngsim
DEFAULT_READ_CONFIG = {'x_fold': 2, 'len_r': 150, 'len_l': 150, 'std_dev': 50, 'dist': 500}

if __name__ == '__main__':

    count = defaultdict(int)
    for k, (i, r1, r2) in enumerate(ngsim.readgen("A.fa.gz", **DEFAULT_READ_CONFIG)):
        count[i.split('|')[0]] += 1
        print(k, end=' ')

    print(count)
