import sys
import numpy as np
import subprocess as sp
import json
import time
import datetime
from pathlib import Path

rowcnts = 2**np.arange(10,15)
colcnts = [128, 256, 512]
# colcnts = 2**np.arange(10,14)
pvals = [2, 5, 10]
bvals = [4, 8]
# bvals = [4, 8, 16]

for m in rowcnts:
    for n in colcnts:
        for p in pvals:
            for nprocs in bvals:
                failed = False
                payload = dict()
                payload["timestamp"] = str(datetime.datetime.now().timestamp())
                payload["nprocs"] = int(nprocs)
                payload["rows"] = int(m)
                payload["cols"] = int(n)
                payload["ncomps"] = int(p)
                payload["seed"] = np.random.randint(1<<31)
                t1 = -time.perf_counter()
                # tcpunch = sp.Popen(["fmi/extern/TCPunch/server/build/tcpunchd"], stdout=sp.PIPE, stderr=sp.PIPE)
                procs = [None]*nprocs
                for rank in range(nprocs):
                    payload["myrank"] = int(rank)
                    payload_str = json.dumps(payload)
                    # cmd = f"./build/fmi_svd_bench {rank} {nprocs} {m} {n} {p} {np.random.randint(1<<31)}"
                    cmd = f"aws lambda invoke --cli-read-timeout 600 --function-name aws_svd_bench --cli-binary-format raw-in-base64-out "
                    cmd += f"--payload '{payload_str}' out/rows{m}_cols{n}_comps{p}_nprocs{nprocs}.json"
                    sys.stdout.write(cmd + "\n")
                    sys.stdout.flush()
                    # procs[rank] = sp.Popen(cmd.split(), stdout=sp.PIPE, stderr=sp.PIPE)
                # for rank in range(nprocs):
                    # procs[rank].wait()
                    # if procs[rank].returncode != 0:
                        # failed = True
                        # for i in range(rank+1, nprocs):
                            # procs[rank].terminate()
                        # break
                time.sleep(1)
                # tcpunch.terminate()
                t1 += time.perf_counter()
                if not failed:
                    sys.stdout.write(f"[nprocs={nprocs},m={m},n={n},p={p}] finished in {t1:.5f} seconds\n")
                    sys.stdout.flush()
                else:
                    sys.stdout.write(f"[nprocs={nprocs},m={m},n={n},p={p}] failed\n")
                    sys.stdout.flush()
