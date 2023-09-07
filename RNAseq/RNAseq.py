# %%
import re
import subprocess
from pathlib import Path
from multiprocessing import Process, Manager

# cwd = Path.cwd()
th = 6
wkdir = "."
abswd = Path(wkdir).absolute()
rawdata = f"{abswd}/rawdata"
outdir = f"{abswd}/clean"

if not Path(outdir).exists(): Path(outdir).mkdir() 
# %%
if Path(rawdata).exists():
    fq = list(Path(rawdata).iterdir())

sample = list(
    {re.split("\W|_", str(Path(f).name))[0] for f in fq}
)

def gg(s: str, return_dict):
    global th
    global fq
    fin = [str(f) for f in fq
           if str(Path(f).name)
           .startswith(s)]
    fq1 = [i for i in fin
           if re.split("\.", i)[0]
           .endswith("1")][0]
    fq2 = [i for i in fin
           if re.split("\.", i)[0]
           .endswith("2")][0]
    ii = f"fastp --thread {th} --length_required 21 -i {fq1} -I {fq2}"
    oo = f"-o {outdir}/{s}_1.cl.fq.gz -O {outdir}/{s}_2.cl.fq.gz"
    jj = f"-j {outdir}/{s}.fastp.json -h {outdir}/{s}.fastp.html"
    cmd = f"{ii} {oo} {jj}".split(" ")
    res = subprocess.Popen(cmd,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE, text=True).communicate()
    return_dict[res] = res

if __name__ == '__main__':
    manager = Manager()
    return_dict = manager.dict()
    jobs = []
    for _ in sample:
        p = Process(target=gg, args=(_, return_dict))  # 实例化进程对象
        p.start()
        jobs.append(p)

    for proc in jobs:
        proc.join()
    print(return_dict.values())

# %%
