from __future__ import print_function, division
from plumbum import local, cli, TEE
import re
import pandas as pd
import numpy as np

TIMES = re.compile(r"Time = ([-+]?\d*\.\d+|\d+) ms")
SIZE = 6

class TimeIt(cli.Application):

    def main(self):
        PhSp = local['../../../build/PhSp']
        nprocs = int(local['nproc']())

        ranges = range(1,12) + range(12,24,2) + range(24,64,4) + range(64,128,16) + range(128,256,64) + [256]
        ranges = [r for r in ranges if r <= nprocs]

        results = pd.DataFrame(index=ranges, columns='Generate GenStdDev Copy CopyStdDev'.split())

        for j in ranges:
            with local.env(OMP_NUM_THREADS=j):
                times = np.empty([SIZE,2])
                for i in range(SIZE):
                    _, output, _ = PhSp & TEE
                    times[i] = map(float, TIMES.findall(output))
                generate = times[:,0].mean()
                gstd = times[:,0].std()
                copy = times[:,1].mean()
                cstd = times[:,1].std()

            results.loc[j] = [generate, gstd, copy,  cstd]

        print(results)

if __name__ == "__main__":
    TimeIt()
