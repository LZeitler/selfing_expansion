import numpy as np
import pandas as pd
from itertools import product
from datetime import datetime
import os
import subprocess

config={}
config["workpath"]="~/temp/slim/"
here = os.path.expanduser(config["workpath"])
wdir = here  # + now


def expand_grid(data_dict):
    rows = product(*data_dict.values())
    return pd.DataFrame.from_records(rows, columns=data_dict.keys())


# set up replicates
reps = 10
# set up parameters
allparams = {
    "hdel": [.3],
    "hben": [.3],
    "mdel": [-0.001,-0.0001],
    "mben": [0.01,0.000],
    "nneutral": [0.25],
    "ndel": [0.75],
    "nben": [0.01],
    "Na": [5000],
    "Ns": [200],
    "subpops": [200],
    "subpopsA": [25],
    "selfRate": [0,.5,1],
    "sh": [0],
    "mig": [0.05],
    "maxAge": [1,50],
    "rCoeff": [1.2],
}



nreps = np.arange(reps)
dparams = expand_grid(allparams)
dparams.index = range(1000,len(dparams) + 1000)
nrep=7; dparam=1006
d = dparams[dparams.index == dparam]
vals = [d.iloc[0,i] for i in range(d.shape[1])]
rowl = dparams.columns.tolist()
pars = [str(i) + "=" + str(j) for i,j in zip(rowl, vals)]
pars = " -d ".join(pars)
pars = "-d " + pars + " -d parcomb=" + str(d.index[0]) + " -d rep=" + str(nrep)
"slim" + " " + pars + " " + "slimscript"


## check if file exist
os.chdir(wdir)
os.path.isfile("parcomb.csv")


## stages
stages = ["edge-at90","edge-ap50","mids-at90","core-finl"]


## burn in
burnvcf = []
for r in dparamsB.index:
    for n in range(breps):
        s = ""
        for c in bpar.columns:
            s += "_" + str(c) + str(dparamsB.loc[r,c])
        burnvcf.append(str(n) + s + ".slimout")

burnfind = pd.DataFrame({'vcfname':burnvcf,"parcomb":np.repeat(dparamsB.index,breps)})
burnvcf = pd.unique(burnvcf)


# burnin
output = '36_hdel0.3_hben0.3_mdel-0.0001_mben0.0_nneutral0.499_ndel0.5_nben0.01_Na5000_rCoeff1.2.slimout'
pars = [str(i) + "=" + str(j) for i,j in zip(rowl, vals)]
pars = " -d ".join(pars)
pars = "-d " + pars + " -d bfile=\\\"" + output + "\\\" -d rep=" + output.split("_")[0]
pars

# these are the wildcards
wildc = {"nrep":nreps,
         "dparam":dparams.index,
         "stage":stages,
         "burnin":burnvcf,
}


# sim
os.chdir("/home/leo/temp/")
bfile = "\\\"" + burnfind[burnfind.parcomb == int("1006")].sample(1).vcfname.values[0] + "\\\""
bfile
bfile = '\\"16_hdel0.3_hben0.3_mdel-0.0001_mben0.0_nneutral0.499_ndel0.5_nben0.01_Na5000_rCoeff1.2.slimout\\"'
bfile

pars = [str(i) + "=" + str(j) for i,j in zip(rowl, vals)]
pars = " -d ".join(pars)
pars = "-d " + pars + " -d parcomb=" + str(d.index[0]) + " -d bfile=" + bfile + " -d rep=" + str(1)
pars

subprocess.Popen("/home/leo/programs/qtslim/build/slim" + " " + pars + " " + "~/pro/selfing-sims/slim/test_selfing_fecundity_w_burnin.slim", shell=True).wait()
