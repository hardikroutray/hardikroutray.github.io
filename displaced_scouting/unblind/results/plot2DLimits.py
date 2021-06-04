import json, glob, sys, os, gzip, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from scipy.interpolate import interp2d
from scipy.interpolate import griddata

import ROOT
ROOT.gROOT.SetBatch(True)

from ROOT import TTree, TFile
from array import array

model = str(sys.argv[-1])

print("-----Please run in python 3-----")

if model == "BphiX":
    df_obs = pd.read_csv("https://hardikroutray.github.io/displaced_scouting/unblind/results/bphilimits_alldata_allctaus_wstheta_analysisDP_v1.csv")
elif model == "HZdZd_2mu4mu" or model=="HZdZd_2mu4mu_eps":
    df_obs = pd.read_csv("https://hardikroutray.github.io/displaced_scouting/unblind/results/hzdlimits_alldata_allctaus_weps_analysisDP_twomufourmu_v1.csv")
else:
    print("-----Please input last argument as BphiX or HZdZd_2mu4mu or HZdZd_2mu4mu_eps-----")
    exit()

print("Running on model", model)

if model == "HZdZd_2mu4mu" or model == "HZdZd_2mu4mu_eps":
    df_obs = df_obs[df_obs["BRbrULobs"] <= 1]

mask_ranges = np.array([
    [0.43,0.49],
    [0.52,0.58],
    [0.73,0.84],
    [0.96,1.08],
    [2.91,3.27],
    [3.47,3.89],
    [8.99,9.87],
    [9.61,10.77],
    ])
mask_ranges[:,0] = mask_ranges[:,0]/(1+0.05)
mask_ranges[:,1] = mask_ranges[:,1]/(1-0.05)

mask_ranges=np.array([[ 0.40952381,  0.51578947],
 [ 0.4952381,   0.61052632],
 [ 0.6952381,   0.88421053],
 [ 0.91428571,  1.13684211],
 [ 2.77142857,  4.09473684],
 [ 8.56190476, 11.33684211]])


if model == "BphiX":
    x = np.array(df_obs["mass"]).reshape(-1,13)
    y = np.array(df_obs["ctau"]).reshape(-1,13)
    z = np.array(df_obs["BRbrULobs"]).reshape(-1,13)

    # f = interp2d(x, y, z, bounds_error=True)
    # x_new = np.arange(0.3,5,0.001)
    # y_new = np.arange(0.1,100,0.1)
    # X, Y = np.meshgrid(y_new,x_new)

    # print("shape of:", " x", x.shape, " y", y.shape, " z", z.shape, " x_new", x_new.shape, " y_new", y_new.shape, " X", X.shape, " Y", Y.shape)
    # print(x.flatten(),y.flatten())
    # Z = griddata((x.flatten(), y.flatten()), z.flatten(), (X, Y), method='cubic')

elif model == "HZdZd_2mu4mu":
    x = np.array(df_obs["mass"]).reshape(-1,21)
    y = np.array(df_obs["ctau"]).reshape(-1,21)
    z = np.array(df_obs["BRbrULobs"]).reshape(-1,21)

else:
    x = np.array(df_obs["mass"]).reshape(-1,21)
    y = np.array(df_obs["eps"]).reshape(-1,21)
    z = np.array(df_obs["BRbrULobs"]).reshape(-1,21)

fig = plt.figure()

from matplotlib.colors import LogNorm
norm = LogNorm()

# cb = plt.contourf(x,y,z,cmap="Set3",norm=norm)

if model == "BphiX":
    levels = [1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5]
    levvals = ["1em11","1em10","1em9","1em8","1em7","1em6","1em5"]    
    # levels = [1e-11,1e-10]
else:
    levels = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
    levvals = ["1em6","1em5","1em4","1em3","1em2","1em1","1e0"]

cb = plt.contourf(x,y,z,cmap="Set3",norm=norm,levels=levels)

# cb = plt.contourf(X,Y,Z,cmap="Set3",norm=norm)
# cb = plt.pcolormesh(x,y,z,cmap="Set3",norm=norm)
# cb = plt.scatter(x,y,c=z,cmap="Set3",norm=norm)
# cb = plt.contour(x,y,z,levels=[0.00000000001,0.0000000001,0.000000001,0.00000001,0.0000001,0.000001,0.00001],cmap="Set3",norm=norm)
# cb = plt.contour(x,y,z,levels=6,cmap="tab20",norm=norm)


plt.yscale("log")
plt.xscale("log")

if model == "BphiX":
    plt.xlim([0.3,5])
else:
    plt.xlim([0.5,50])

ylim = plt.ylim()
for low,high in mask_ranges:
    plt.fill_betweenx(ylim, low, high, color="0.8", zorder=2)
    
plt.xlabel(r"m$_{\phi}$ [GeV]")
if "eps" in model:
    plt.ylabel(r"$\epsilon$")
else:
    plt.ylabel(r"$c\tau$ [mm]")

plt.title(r"${}_\mathbf{CMS}$ ${}_\mathit{Supplementary}$                         101 fb$^{-1}$(13 TeV)")

if model == "BphiX":
    fig.colorbar(cb, label=r"BR(B $\rightarrow$ $\phi$ X).BR($\phi$ $\rightarrow$ $\mu\mu$)")
else:
    fig.colorbar(cb, label=r"BR(H $\rightarrow$ Z$_{D}$Z$_{D}$)")

fig.patch.set_alpha(1)

fig.savefig("mass_ctau_{}_exclusion.png".format(model), dpi=300)
fig.savefig("mass_ctau_{}_exclusion.pdf".format(model))

# exit()

graphs = dict()

print("The number of levels", len(cb.allsegs))
for lev in range(len(cb.allsegs)):
    print("")
    print("level",lev,"levval",levvals[lev])
    print("The number of curves in this level", len(cb.allsegs[lev]))
    print("")

    for cur in range(len(cb.allsegs[lev])):
        # print("The number of points in this curve", len(cb.allsegs[lev][cur]))
        xval_ = cb.allsegs[lev][cur][:,0]
        yval_ = cb.allsegs[lev][cur][:,1]
        # print(xval_,yval_)

        xval = []
        yval = []
        for val in range(len(xval_)):
            if any((sublist[0] <= xval_[val] <= sublist[1]) for sublist in mask_ranges):
                continue
            xval.append(xval_[val])
            yval.append(yval_[val])
        # print(xval,yval)

        g = ROOT.TGraph(len(xval))
        g.SetName("gr_{}_{}".format(levvals[lev],cur+1))
        for i, (xv,yv) in enumerate(zip(xval,yval)):
            g.SetPoint(i,xv,yv)
        graphs["gr_{}_{}".format(levvals[lev],cur+1)] = g

if model == "BphiX":
    outfile = TFile('BphiX_UL.root', 'recreate')
elif model == "HZdZd_2mu4mu":
    outfile = TFile('HZdZd_UL.root', 'recreate')
else:
    outfile = TFile('HZdZd_eps_UL.root', 'recreate')
outfile.cd()
 
for k,v in graphs.items():
    v.Write(k)
outfile.Close()

exit()
#######################store values in root file######################

if model == "BphiX":
    outfile = TFile('BphiX_UL.root', 'recreate')
else:
    outfile = TFile('HZdZd_UL.root', 'recreate')

outfile.cd()

if model == "BphiX":
    ttree = TTree('ttree', 'ttree with BphiX UL (mass in GeV, ctau in mm)')
else:
    ttree = TTree('ttree', 'ttree with HZdZd UL (mass in GeV, ctau in mm)')

mass = array("f", [0.0])
ctau = array("f", [0.0])
ul = array("d", [0.0])

ttree.Branch("mass", mass, "mass/F")
ttree.Branch("ctau", ctau, "ctau/F")
ttree.Branch("ul", ul, "ul/D")

masses = df_obs["mass"].tolist()
ctaus = df_obs["ctau"].tolist()
uls = df_obs["BRbrULobs"].tolist()

num = 0
for ev in range(len(masses)):
    num+= 1

    mass[0] = round(masses[ev],3)
    ctau[0] = round(ctaus[ev],1)
    ul[0] = round(uls[ev],15)
    # ul[0] = uls[ev]
    
    # print(mass[0],ctau[0],ul[0])
    print(round(masses[ev],3),round(ctaus[ev],1),round(uls[ev],15))
        
    ttree.Fill()

outfile.Write()
outfile.Close()
