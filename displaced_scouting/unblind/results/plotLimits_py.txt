import os,sys
import numpy as np
import ROOT
from collections import OrderedDict

mask_ranges = [
    [0.43,0.49],
    [0.52,0.58],
    [0.73,0.84],
    [0.96,1.08],
    [2.91,3.27],
    [3.47,3.89],
    [8.99,9.87],
    [9.61,10.77],
    ]
mask_ranges = [[ mr[0]/(1+0.05), mr[1]/(1-0.05) ] for mr in mask_ranges]

nRows = 15 # number of values per mass point
nCols = 6 # number of columns per limit
if len(sys.argv)<4:
    print("Please, specify input csv file, model (B,H2,H4) and lifetime (1,10,100)")
    exit(1)

fin   = sys.argv[1]
model = sys.argv[2]
ctau  = int(sys.argv[3])

if model=="B":
    model = "pp #rightarrow B #rightarrow #phi X #rightarrow 2#mu X (c#tau_{0}^{#phi}=%d mm)"%ctau
elif model=="H2":
    model = "gg #rightarrow h #rightarrow Z_{D}Z_{D} #rightarrow 2#mu 2X (c#tau_{0}^{Z_{D}}=%d mm)"%ctau
elif model=="H4":
    model = "gg #rightarrow h #rightarrow Z_{D}Z_{D} #rightarrow 4#mu (c#tau_{0}^{Z_{D}}=%d mm)"%ctau
else:
    print("Unknown model (B,H2,H4)")
    exit(1)

def smooth_limits(limits):
    mex = []
    window = 3 # times 1%*mass
    nsigma = 1.5 # if |br - median(window-except-central)|/std(window-except-central) > nsigma, then filter out
    whichlimits  = ['BRbrUL','BRbrULup','BRbrULupup','BRbrULdown','BRbrULdowndown','BRbrULobs']
    for w in whichlimits:
        thisl={}
        for k in limits.keys():
            thisl[k]=limits[k][w]

        for k in thisl.keys():
            subl = [vv for kk,vv in thisl.items() if ((k > kk*(1-window*0.011)) and (k < kk*(1+window*0.011)) and not (kk==k))]
            if abs(np.median(subl))>0.0:
                subl = [vv for vv in subl if (vv/np.median(subl)<10.0 and vv/np.median(subl)>1.0/10.0)]
            if abs(np.std(subl))>0.0 and abs(thisl[k] - np.median(subl))/np.std(subl) > nsigma:
                mex.append(k)

    return mex

def plotSignificance(mv,lv,twosided=True,outname=""):

    if twosided:
        minS = -5.0
        maxS = +5.0
        outname = outname+"significance_twosided"
    else:
        minS = 0.0
        maxS = 5.0
        outname = outname+"significance_onesided"
        
    haxiss=ROOT.TH2D("haxiss","",1,np.amin(mv)*0.99,np.amax(mv)+1,1,minS,maxS)
    if "H" in sys.argv[2]:
        haxiss.GetXaxis().SetTitle("m_{Z_{D}} [GeV]")
    else:
        haxiss.GetXaxis().SetTitle("m_{#phi} [GeV]")
    haxiss.GetYaxis().SetTitle("Significance [#sigma]")
    if not twosided:
        haxiss.GetYaxis().SetNdivisions(505)
    haxiss.GetXaxis().SetMoreLogLabels()
    haxiss.GetXaxis().SetNoExponent()
    haxiss.GetXaxis().SetLabelSize(0.035)
    haxiss.GetXaxis().SetLabelOffset(0.003)
    haxiss.GetYaxis().SetLabelSize(0.04)
    haxiss.GetXaxis().SetTitleSize(0.047)
    haxiss.GetYaxis().SetTitleSize(0.047)
    haxiss.GetXaxis().SetTitleOffset(0.85)
    
    box=[]
    for r,m in enumerate(mask_ranges):
        box.append(ROOT.TBox(mask_ranges[r][0],minS+0.001*abs(minS),mask_ranges[r][1],maxS-0.001*maxS))
        box[r].SetFillColor(ROOT.kGray)
        box[r].SetLineColor(ROOT.kGray)
        box[r].SetLineWidth(0)
        box[r].SetFillStyle(1001)

    go = ROOT.TGraph(len(mv),mv,lv)

    ROOT.gStyle.SetOptStat(0)
    can=ROOT.TCanvas("can","",1200,600)
    can.cd()
    ROOT.gPad.SetLogx()
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()

    haxiss.Draw()
    go.Draw("PLsame")

    for b in box:
        if b.GetX1()>=(np.amin(mv)*0.99):
            b.Draw("same")
    
    ROOT.gPad.RedrawAxis()

    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(0.04)
    latex.SetNDC(True)
    latexCMS = ROOT.TLatex()
    latexCMS.SetTextFont(61)
    latexCMS.SetTextSize(0.055)
    latexCMS.SetNDC(True)
    latexCMSExtra = ROOT.TLatex()
    latexCMSExtra.SetTextFont(52)
    latexCMSExtra.SetTextSize(0.03)
    latexCMSExtra.SetNDC(True) 

    yearenergy="101 fb^{-1} (13 TeV)"
    latex.DrawLatex(0.9, 0.91, yearenergy)

    latexCMS.DrawLatex(0.1,0.91,"CMS")
    cmsExtra="Preliminary"    
    latexCMSExtra.DrawLatex(0.16,0.91,"%s"%(cmsExtra))				

    latexModel = ROOT.TLatex()
    latexModel.SetTextFont(42)
    latexModel.SetTextSize(0.05)
    latexModel.SetNDC(True)
    latexModel.DrawLatex(0.13, 0.825, "%s"%(model))

    can.Update()
    can.SaveAs("%s_%s_%smm.png" % (outname,sys.argv[2],sys.argv[3]))
    can.SaveAs("%s_%s_%smm.pdf" % (outname,sys.argv[2],sys.argv[3]))

###

def plotUpperLimit(mv,lov,lev,lp1v,lp2v,lm1v,lm2v,outname="upperlimit_br"):
    
    lumi=101

    maxL = 2.0*np.amax(lp2v)
    minL = 0.5*np.amin(lm2v)
    
    haxisl=ROOT.TH2D("haxiss","",1,np.amin(mv)*0.99,np.amax(mv)+1,1,minL,maxL)
    ytitle=""
    if "H" in sys.argv[2]:
        haxisl.GetXaxis().SetTitle("m_{Z_{D}} [GeV]")
        if "H2" in sys.argv[2]:
            ytitle="#it{B}(h #rightarrow Z_{D}Z_{D}) #it{B}(Z_{D} #rightarrow #mu#mu)"
        else:
            ytitle="#it{B}(h #rightarrow Z_{D}Z_{D})"
    else:
        haxisl.GetXaxis().SetTitle("m_{#phi} [GeV]")
        ytitle="#it{B}(B #rightarrow #phi X) #it{B}(#phi #rightarrow #mu#mu)"
    haxisl.GetYaxis().SetTitle(ytitle)
    haxisl.GetXaxis().SetMoreLogLabels()
    haxisl.GetXaxis().SetNoExponent()
    haxisl.GetXaxis().SetLabelSize(0.035)
    haxisl.GetXaxis().SetLabelOffset(0.003)
    haxisl.GetYaxis().SetLabelSize(0.04)
    haxisl.GetXaxis().SetTitleSize(0.047)
    haxisl.GetYaxis().SetTitleSize(0.047)
    haxisl.GetXaxis().SetTitleOffset(0.85)
    
    box=[]
    for r,m in enumerate(mask_ranges):
        box.append(ROOT.TBox(mask_ranges[r][0],minL+0.001*abs(minL),mask_ranges[r][1],maxL-0.001*maxL))
        box[r].SetFillColor(ROOT.kGray)
        box[r].SetLineColor(ROOT.kGray)
        box[r].SetLineWidth(0)
        box[r].SetFillStyle(1001)

    ROOT.gStyle.SetOptStat(0)
    can = ROOT.TCanvas("can", "", 1200, 600)
    can.cd()
    ROOT.gPad.SetLogx()
    ROOT.gPad.SetLogy() 
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()

    haxisl.Draw()
        
    gr_s2b = ROOT.TGraphAsymmErrors(len(mv),mv,lev,mv-mv,mv-mv,lev-lm2v,lp2v-lev)
    gr_s2b.SetFillColor(ROOT.kOrange)
    gr_s2b.SetLineColor(0)
    gr_s2b.Draw("3")
      
    gr_s1b = ROOT.TGraphAsymmErrors(len(mv),mv,lev,mv-mv,mv-mv,lev-lm1v,lp1v-lev) 
    gr_s1b.SetFillColor(ROOT.kGreen+1)
    gr_s1b.SetLineColor(0)
    gr_s1b.Draw("3")
    
    gexp = ROOT.TGraph(len(mv),mv,lev)
    gexp.SetLineStyle(7)
    gexp.SetLineWidth(3)
    gexp.SetLineColor(ROOT.kBlack)
    gexp.SetLineStyle(3)
    gexp.SetLineWidth(3)
    gexp.SetLineColor(ROOT.kRed)
    gexp.Draw("L")
    
    gobs = ROOT.TGraph(len(mv),mv,lov)
    gobs.SetMarkerStyle(ROOT.kFullCircle)
    gobs.SetMarkerSize(1.5)
    gobs.SetMarkerColor(ROOT.kBlack)
    gobs.SetLineWidth(3)
    gobs.SetLineWidth(1)
    gobs.SetLineColor(ROOT.kBlack)
    gobs.Draw("L")
   
    for b in box:
        if b.GetX1()>=(np.amin(mv)*0.99):
            b.Draw("same")

    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(0.04)
    latex.SetNDC(True)
    latexCMS = ROOT.TLatex()
    latexCMS.SetTextFont(61)
    latexCMS.SetTextSize(0.055)
    latexCMS.SetNDC(True)
    latexCMSExtra = ROOT.TLatex()
    latexCMSExtra.SetTextFont(52)
    latexCMSExtra.SetTextSize(0.03)
    latexCMSExtra.SetNDC(True) 

    yearenergy="101 fb^{-1} (13 TeV)"
    latex.DrawLatex(0.9, 0.91, yearenergy)

    latexCMS.DrawLatex(0.1,0.91,"CMS")
    cmsExtra="Preliminary"    
    latexCMSExtra.DrawLatex(0.16,0.91,"%s"%(cmsExtra))				

    latexModel = ROOT.TLatex()
    latexModel.SetTextFont(42)
    latexModel.SetTextSize(0.05)
    latexModel.SetNDC(True)
    latexModel.DrawLatex(0.13, 0.825, "%s"%(model))
    
    l1 = ROOT.TLegend(0.625, 0.575, 0.875, 0.875)
    if sys.argv[2]=="B" and int(sys.argv[3])==1:
        l1 = ROOT.TLegend(0.525, 0.125, 0.775, 0.425)
    l1.SetTextFont(42)
    l1.SetTextSize(0.036)
    l1.SetLineColor(ROOT.kWhite)
    l1.SetShadowColor(ROOT.kWhite)
    l1.SetFillColor(ROOT.kWhite)

    l1.AddEntry(gobs , "Observed limit (95% CL)", "l")
    l1.AddEntry(gexp , "Median expected limit", "l")
    l1.AddEntry(gr_s1b , "68% expected", "f")
    l1.AddEntry(gr_s2b , "95% expected", "f")
    l1.Draw()

    ROOT.gPad.RedrawAxis()
        
    can.Update()

    can.SaveAs("%s_%s_%smm.png" % (outname,sys.argv[2],sys.argv[3]))
    can.SaveAs("%s_%s_%smm.pdf" % (outname,sys.argv[2],sys.argv[3]))

###

fi = open(fin)
li = fi.readlines()

nP=0 # number of mass points
nR=0 # number of rows per mass point
limits = dict()

C=0
if ctau==1:
    C=1
elif ctau==10:
    C=2
elif ctau==100:
    C=3
elif ctau==1000:
    C=4

for nl,l in enumerate(li):
    if "what" in l or l.startswith("#"):
        continue
    nR = nR+1

    if nR==1:
        mass=float(l.strip("\n").split(",")[0])
        limits[mass] = {}
        
    what=str(l.strip("\n").split(",")[nCols-1])
    try:
        lim =float(l.strip("\n").split(",")[C])
    except:
        print("value missing in this line for some reason: ", l)
    limits[mass][what] = lim

    if nR>=nRows:
        nR=0
        nP=nP+1

limits  = OrderedDict(sorted(limits.items()))

### Plot
ml=[]
stl=[]
sol=[]
lol=[]
lel=[]
lp1l=[]
lp2l=[]
lm1l=[]
lm2l=[]

dropLargeSig=True
largeSig=2.5
dropSpikes=True

whichLimit  ='BRbrUL'
whichLimit1up = whichLimit+"up"
whichLimit2up = whichLimit+"upup"
whichLimit1dn = whichLimit+"down"
whichLimit2dn = whichLimit+"downdown"
whichLimitObs = whichLimit+"obs"

mex = smooth_limits(limits)

for k in limits.keys():
    if dropLargeSig and abs(limits[k]['Significance'])>largeSig:
        continue
    if dropSpikes and k in mex:
        continue
    ml.append(k)
    if abs(limits[k]['Significance'])>1e-3:
        stl.append(limits[k]['Significance'])
        sol.append(max(0.0,limits[k]['Significance']))
    else:
        stl.append(0.0)
        sol.append(0.0)
    lel .append(limits[k][whichLimit   ])
    lp1l.append(limits[k][whichLimit1up])
    lp2l.append(limits[k][whichLimit2up])
    lm1l.append(limits[k][whichLimit1dn])
    lm2l.append(limits[k][whichLimit2dn])
    lol .append(limits[k][whichLimitObs])

mv   = np.array(ml  )
lev  = np.array(lel )
lp1v = np.array(lp1l)
lp2v = np.array(lp2l)
lm1v = np.array(lm1l)
lm2v = np.array(lm2l)
lov  = np.array(lol )
stv  = np.array(stl )           
sov  = np.array(sol )           
    
plotSignificance(mv,stv,1)
plotSignificance(mv,sov,0)
plotUpperLimit(mv,lov,lev,lp1v,lp2v,lm1v,lm2v)


