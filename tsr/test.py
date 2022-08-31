# Test for the extension.
#
# Eli Bendersky [http://eli.thegreenplace.net]
# This code is in the public domain.
from collections import defaultdict
import ngsim

if __name__ == '__main__':

    count = defaultdict(int)
    for i, r1, r2 in ngsim.readgen("tsr/A.fa.gz", x_fold=2, dist=25):
        count[i] += 1
        print(i, r1, r2)

    print(count)
