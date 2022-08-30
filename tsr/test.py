# Test for the extension.
#
# Eli Bendersky [http://eli.thegreenplace.net]
# This code is in the public domain.
import ngsim

if __name__ == '__main__':

    for elem in ngsim.readgen("A.fa.gz", 100, dist=25):
        print(elem)
