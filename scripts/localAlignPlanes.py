import os
import awkward as ak
import uproot
import numpy as np
import sys
import io
from alignmentTable import *

filename = sys.argv[1]
oldalignfn = sys.argv[2]
newalignfn = sys.argv[3]

file = uproot.open(filename)
t = file["EventNtuple/ntuple"].arrays(filter_name="trkhits")

oldalign = AlignmentSet()
oldalign.fromFile(oldalignfn)
newalign = AlignmentSet()
newalign.fromFile(oldalignfn)

planederiv = ak.flatten(ak.flatten(t['trkhits']['planeDeriv'])).to_numpy()
udresid = ak.flatten(ak.flatten(t['trkhits']['udresid'])).to_numpy()
udresidmvar = ak.flatten(ak.flatten(t['trkhits']['udresidmvar'])).to_numpy()
udresidpvar = ak.flatten(ak.flatten(t['trkhits']['udresidpvar'])).to_numpy()
var = udresidmvar+udresidpvar
plane = ak.flatten(ak.flatten(t['trkhits']['plane'])).to_numpy()
state = ak.flatten(ak.flatten(t['trkhits']['state'])).to_numpy()

planederiv *= -1

nparams = 6
dp_num = np.array([np.zeros((nparams,1)) for i in range(36)])
dp_den = np.array([np.zeros((nparams,nparams)) for i in range(36)])

chi2 = 0
counts = 0

for i in range(len(plane)):
    if np.abs(state[i]) == 1: 
        jacob = planederiv[i][np.newaxis,:]
        dp_den[plane[i]] += np.matmul(jacob.T,jacob)/var[i]
        dp_num[plane[i]] += jacob.T*udresid[i]/var[i]
        counts += 1
        chi2 += udresid[i]**2/var[i]
print(chi2/counts,chi2,counts)
for i in range(36):
    dp = np.matmul(np.linalg.inv(dp_den[i]),dp_num[i])[:,0]    
    for j in range(nparams):
        newalign.tables["TrkAlignPlane"].constants[i][j] -= dp[j]
fout = open(newalignfn,"w")
fout.write(newalign.tables["TrkAlignPlane"].to_proditions_table())
fout.close()

oldalign.getStats()
newalign.getStats()
print("OLD ALIGNMENT:")
oldalign.printStats()
print("NEW ALIGNMENT:")
newalign.printStats()
