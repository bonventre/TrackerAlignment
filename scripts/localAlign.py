import os
import uproot
import awkward as ak
import numpy as np
import sys
import io
from alignmentTable import *

def process(filename,oldalignfn,newalignfn,minplanes,minpanels,minhits,mincounts,alignplane=True,rotations=True):
    file = uproot.open(filename)
    t = file["EventNtuple/ntuple"].arrays(filter_name=["trkhits","trkhitcalibs","trk.*"])
    cut = (t["trk.nactive"] == t["trk.nhits"]) & (t["trk.nplanes"] > minplanes) & (t["trk.npanels"] > minpanels) & (t["trk.nactive"] > minhits)
    t = t[cut]
    
    oldalign = AlignmentSet()
    oldalign.fromFile(oldalignfn)
    newalign = AlignmentSet()
    newalign.fromFile(oldalignfn)
    
    planederiv = ak.flatten(ak.flatten(t['trkhitcalibs']['dDdPlane'])).to_numpy()
    panelderiv = ak.flatten(ak.flatten(t['trkhitcalibs']['dDdPanel'])).to_numpy()
    udresid = ak.flatten(ak.flatten(t['trkhits']['udresid'])).to_numpy()
    udresidmvar = ak.flatten(ak.flatten(t['trkhits']['udresidmvar'])).to_numpy()
    udresidpvar = ak.flatten(ak.flatten(t['trkhits']['udresidpvar'])).to_numpy()
    plane = ak.flatten(ak.flatten(t['trkhits']['plane'])).to_numpy()
    panel = ak.flatten(ak.flatten(t['trkhits']['panel'])).to_numpy()
    state = ak.flatten(ak.flatten(t['trkhits']['state'])).to_numpy()

    planederiv *= -1
    panelderiv *= -1
    
    if alignplane:
        if rotations:
            nparams = 6
        else:
            nparams = 3
        npstart = 0
        nchannels = 36
        deriv = planederiv
        channel = plane
        table = "TrkAlignPlane"
    else:
        if rotations:
            nparams = 5
        else:
            nparams = 2
        npstart = 1
        nchannels = 216
        deriv = panelderiv[:,1:]
        channel= plane*6+panel
        table = "TrkAlignPanel"

    dp_num = np.array([np.zeros((nparams,1)) for i in range(nchannels)])
    dp_den = np.array([np.zeros((nparams,nparams)) for i in range(nchannels)])
    dp = np.array([np.zeros(nparams) for i in range(nchannels)])
    variance = np.array([np.zeros((nparams,nparams)) for i in range(nchannels)])
    
    for ifit in range(1):
        chi2 = np.zeros(nchannels)
        counts = np.zeros(nchannels)
        
        for i in range(len(channel)):
          if state[i] > -2:
            jacob = np.array([[deriv[i][k] for k in range(nparams)]])
            v = udresidmvar[i]+udresidpvar[i]
            r = udresid[i] - np.dot(dp[channel[i]],jacob[0])
            dp_den[channel[i]] += np.matmul(jacob.T,jacob)/v
            dp_num[channel[i]] += jacob.T*r/v
            counts[channel[i]] += 1
            chi2[channel[i]] += r*r/v
        for i in range(nchannels):
          if counts[i] < mincounts:
              continue
        #  dp[i] = dp_num[i]/dp_den[i]
          dp[i] += np.matmul(np.linalg.inv(dp_den[i]),dp_num[i])[:,0]
          variance[i] = np.linalg.inv(dp_den[i])
      
        print(ifit,chi2/counts)
    
    for i in range(nchannels):
        for j in range(nparams):
            newalign.tables[table].constants[i][j+npstart] -= dp[i][j]
    
    fout = open(newalignfn,"w")
    fout.write(newalign.tables["TrkAlignPlane"].to_proditions_table())
    fout.write(newalign.tables["TrkAlignPanel"].to_proditions_table())
    fout.close()
    
    oldalign.getStats()
    newalign.getStats()
    oldalign.printStats()
    newalign.printStats()

if __name__ == "__main__":
    artfcl = sys.argv[1]
    ntuplefcl = sys.argv[2]
    source = sys.argv[3]
    iterations = int(sys.argv[4])
    minplanes = 2
    minpanels = 2
    minhits = 10
    mincounts = 100
    for i in range(iterations):
        os.system("cp iter%d.txt current.txt" % (i))
        os.system(f"mu2e -c {artfcl} -S {source} --output iter{i}.art")
        os.system(f"mu2e -c {ntuplefcl} -s iter{i}.art -T iter{i}.root")
        process(f"iter{i}.root","current.txt",f"iter{i+1}_pl.txt",minplanes,minpanels,minhits,mincounts,alignplane=True)
        process(f"iter{i}.root",f"iter{i+1}_pl.txt",f"iter{i+1}.txt",minplanes,minpanels,minhits,mincounts,alignplane=False)
