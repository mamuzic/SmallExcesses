##################################################################
# Program for cathegorisation of small excesses using pMSSM scan #
# To run Bino, Wino or Higgsino case:                            # 
# python small_excesses_in_pmssm -m <bino,wino,higgsino>         #
# Contact: judita.mamuzic@cern.ch and brian.petersen@cern.ch     #
##################################################################
#!/usr/bin/python
import sys, getopt
import copy
import math,os
import os.path
import json
import ROOT
from ROOT import TFile,TGraph, TCanvas, TH1, TH1I, TH1F, TH1D, TH2F, TH2D, TH3D, TH3F, gStyle, TLegend, TLine, TF1, TMath, TColor, THStack, TAttMarker, gROOT
from math import sqrt
import subprocess
import time
from scipy import arange
import numpy as np
gROOT.SetBatch()

plotsdir = "plots/"

fnames = {}
fnames["bino"]     = "input/Reduced-Bino-NoLongLived.root"
fnames["wino"]     = "input/Reduced-Wino-NoLongLived.root"
fnames["higgsino"] = "input/Reduced-Higgsino-NoLongLived.root"

# Add later range and units
parameters = ["mdR", "muR", "msR", "mcR", "mbR", "mtR",       # right-handed squark mass scales
              "mqL1", "mqL2", "mqL3",                         # left-handed squark mass scales
              "meR", "mmuR", "mtauR", "meL", "mmuL", "mtauL", # slepton mass scales
              "M_1", "M_2", "M_3", "mu",                      # gaugino and Higgsino mass scales
              "Ab", "At", "Atau",                             # trilinear couplings
              "mA", "tanb",                                   # Higgs sector parameters
              ]

particle_masses = ["m_h", "m_W", "m_b",                                  # Higgs, W-boson and b-quark mass
                   "m_d_L", "m_d_R", "m_u_L", "m_u_R", "m_s_L", "m_s_R",
                   "m_c_L", "m_c_R", "m_b_1", "m_b_2", "m_t_1", "m_t_2", # squark masses
                   "m_gl",                                               # gluino mass
                   "m_e_L", "m_e_R", "m_mu_L", "m_mu_R", "m_tau_1",
                   "m_tau_2", "m_nu_eL", "m_nu_muL", "m_nu_tauL",        # slepton masses
                   "m_chi_10", "m_chi_20", "m_chi_30", "m_chi_40",       # neutralino masses
                   "m_chi_1p", "m_chi_2p",                               # chargino masses
                   "m_A, m_H", "m_Hp",                                   # Heavy Higgs masses 
                   ]

# There can be only one (production): "QCD", "EW", "Mix", "test" 
# Production, Analysis, SRs, p0
def getProduction(prod):
    doEW, doQCD, doMix, doTest = False, False, False, False
    if prod == "QCD":    doQCD, doEW, doMix, doTest = True,  False, False, False
    elif prod == "EW":   doQCD, doEW, doMix, doTest = False, True,  False, False
    elif prod == "Mix":  doQCD, doEW, doMix, doTest = False, False, True,  False
    elif prod == "test": doQCD, doEW, doMix, doTest = False, False, False, True
    else: print "ERROR: Unknown production", prod
    
    production = {} 
    if doQCD:
        production["QCD"] = {} # Analysis, SR, p0
        production["QCD"]["ZeroLeptonStrong"] = {"SR2jl": 0.50 , "SR2jm": 0.49 , "SR2jt": 0.29 , "SR3jt": 0.24  ,
                                                 "SR4jl": 0.50 , "SR4jm": 0.50 , "SR4jt": 0.50 , "SR4jvl": 0.35 , "SR5jt": 0.50 ,
                                                 "SR6jl": 0.27 , "SR6jm": 0.25 , "SR6jt": 0.50 , "SR6jvt": 0.36 }
        production["QCD"]["MultijetAnalysis"] = {"SR_10j50": 0.13   ,"SR_10j50_MJ340": 0.8 ,"SR_10j50_MJ420": 0.7 ,"SR_7j80_0b":    0.5  ,"SR_7j80_1b": 0.6    ,"SR_7j80_2b": 0.8 ,
                                                 "SR_8j50_0b": 0.24 ,"SR_8j50_1b": 0.5     ,"SR_8j50_2b": 0.7     ,"SR_8j50_MJ340": 0.60 ,"SR_8j50_MJ420": 0.7 ,"SR_8j80_0b": 0.19 ,"SR_8j80_1b": 0.6 ,"SR_8j80_2b": 0.5 ,
                                                 "SR_9j50_0b": 0.21 ,"SR_9j50_1b": 0.28    ,"SR_9j50_2b": 0.6     ,"SR_9j50_MJ340": 0.7 ,"SR_9j50_MJ420": 0.6 }
        production["QCD"]["ThreeBjets"]       = {"SR_0l_4j_A": 0.40 ,"SR_0l_4j_B": 0.15 ,"SR_0l_4j_C": 0.5 ,
                                                 "SR_0l_7j_A": 0.50 ,"SR_0l_7j_B": 0.50 ,"SR_0l_7j_C": 0.46 ,
                                                 "SR_1l_6j_A": 0.50 ,"SR_1l_6j_B": 0.50 ,"SR_1l_6j_C": 0.50 }
        production["QCD"]["OneLeptonStrong"]  = {"Hard": 0.27 ,"Soft": 0.48 } # Taking the smallest p0
        production["QCD"]["MonoJetFinal"]     = {"SR1": 0.50 , "SR2": 0.50 , "SR3": 0.50 , "SR4": 0.50 , "SR5": 0.50 , "SR6": 0.50 , "SR7": 0.50 , "SR8": 0.41 , "SR9": 0.05} # No p0, I calculated with HF
        production["QCD"]["SameSignThreeLepton"] = {"SameSignThreeLepton": 0.03 } # Taking the smallest p0
        production["QCD"]["TauStrong"]        = {"TauStrong": 0.04 } # Taking the smallest p0
        production["QCD"]["StopZeroLepton"]   = {"SRA1": 0.50 ,"SRA2": 0.50,"SRA3": 0.35 ,"SRA4": 0.20 ,"SRB1": 0.50 ,"SRC1": 0.50 ,"SRC2": 0.50 ,"SRC3": 0.50} # No p0, I calculated with HF
        production["QCD"]["StopOneLepton"]    = {"SRbC0": 0.27 ,"SRbC1_shape": 0.29, "SRbC4": 0.09 ,"SRbC5": 0.36,
                                                 "SRbCvW_shape": 0.12, "SRbWN_shape": 0.5, "SRtN1p_shape": 0.35, "SRtN2x": 0.5, "SRtN3p":  0.5, "SRtNbC": 0.13 } # Sophio SRbC0 is bCc_diag, bC1_shape is bCd_bulk, bC4 and bC5 are bCd_high1 and bCd_high2, SRbCvW_shape is bCb_med1, SRbWN_shape is 3body, SRtN1p_shape is tN_diag, SRtN2x is tN_med, SRtN3p is tN_high, SRtNbC is tNbC_mix    
        production["QCD"]["TwoBjets"]         = {"SRA_mCT150": 0.31 ,"SRA_mCT200": 0.16 ,"SRA_mCT250": 0.50 ,"SRA_mCT300": 0.34 ,"SRA_mCT350": 0.39 ,"SRB": 0.47} # No p0, I calculated with HF
        production["QCD"]["StopTwoLepton"]    = {"SR_2": 0.5,"SR_3": 0.31 ,"SR_4": 0.50 ,"SR_5": 0.50} # No p0, I calculated with HF
        production["QCD"]["StopMonoJet"]      = {"M1": 0.50 ,"M2": 0.50 ,"M3": 0.49 } 
        production["QCD"]["StopZ"]            = {"StopZ": 0.30 } # Taking the smallest p0
        production["QCD"]["StopTB"]           = {"SRexA": 0.11,"SRinA": 0.12,"SRinB": 0.24,"SRinC": 0.06} # No p0, I calculated with HF
        production["QCD"]["EwkThreeLepton"]   = {"EwkThreeLepton": 0.02 } 
        production["QCD"]["EwkFourLepton"]    = {"EwkFourLepton": 0.10 }
    elif doEW:
        production["EW"] = {} # Analysis, SR, p0
        production["EW"]["ZeroLeptonStrong"]  = {"SR2jl": 0.50 , "SR2jm": 0.49 , "SR2jt": 0.29 , "SR3jt": 0.24  ,
                                                 "SR4jl": 0.50 , "SR4jm": 0.50 , "SR4jt": 0.50 , "SR4jvl": 0.35 , "SR5jt": 0.50 ,
                                                 "SR6jl": 0.27 , "SR6jm": 0.25 , "SR6jt": 0.50 , "SR6jvt": 0.36 }
        production["EW"]["MonoJetFinal"]      = {"SR1": 0.50 , "SR2": 0.50 , "SR3": 0.50 , "SR4": 0.50 , "SR5": 0.50 , "SR6": 0.50 , "SR7": 0.50 , "SR8": 0.41 , "SR9": 0.05} # No p0, I calculated with HF
        production["EW"]["EwkTwoLepton"]      = {"SR_WWa": 0.50 ,"SR_WWb": 0.50 ,"SR_WWc": 0.31 ,"SR_Zjets": 0.50 ,"SR_mT2a": 0.50 ,"SR_mT2b": 0.27 ,"SR_mT2c": 0.21 }
        production["EW"]["EwkThreeLepton"]    = {"EwkThreeLepton": 0.02 }
        production["EW"]["EwkFourLepton"]     = {"EwkFourLepton": 0.10 }
        production["EW"]["StopTwoLepton"]     = {"SR_2": 0.5,"SR_3": 0.31 ,"SR_4": 0.50 ,"SR_5": 0.50} # No p0, I calculated with HF
    elif doMix:
        production["Mix"] = {} # Analysis, SR, p0
        production["Mix"]["ZeroLeptonStrong"] = {"SR2jl": 0.50 , "SR2jm": 0.49 , "SR2jt": 0.29 , "SR3jt": 0.24  ,
                                                 "SR4jl": 0.50 , "SR4jm": 0.50 , "SR4jt": 0.50 , "SR4jvl": 0.35 , "SR5jt": 0.50 ,
                                                 "SR6jl": 0.27 , "SR6jm": 0.25 , "SR6jt": 0.50 , "SR6jvt": 0.36 }
        production["Mix"]["MonoJetFinal"]     = {"SR1": 0.50 , "SR2": 0.50 , "SR3": 0.50 , "SR4": 0.50 , "SR5": 0.50 , "SR6": 0.50 , "SR7": 0.50 , "SR8": 0.41 , "SR9": 0.05} # No p0, I calculated with HF 
    elif doTest:
        production = {} # Analysis, SR, p0
        production["QCD"] = {}
        production["QCD"]["EwkThreeLepton"]   = {"EwkThreeLepton": 0.02 }
    return production

        
def listAllBranches(fname):
    if not os.path.isdir(plotsdir): os.system("mkdir "+plotsdir)
    f = openFile(fname)
    t = f.Get("susy")
    branches = t.GetListOfBranches()
    bnames = []
    for br in branches:
        bname = br.GetName()
        bnames.append(bname)
        bname.split("_")        
        print bname
    
def findSmallExcesses(fname, lsp, prod):
    if not os.path.isdir(plotsdir): os.system("mkdir "+plotsdir)
    outfile = TFile(plotsdir+prod+"_"+lsp+".root", "RECREATE")
    # Use production
    production = getProduction(prod)
    f = openFile(fname)
    susy = f.Get("susy")
    susy.SetBranchStatus("*Cat*",0)
    susy.SetBranchStatus("*LL*",0)
    susy.SetBranchStatus("*RecoEvents*",0)
    susy.SetBranchStatus("*Higgs*",0)
    
    # Root gSystem->Load("libSusyFitter.so"), StatTools::GetProbFromSigma()
    #p0threshold = 1.0 # All, should be identical for Z = 0 sigma for the small excesses condition
    #p0threshold = 0.5                      # Z = 0 sigma, Model independent 
    #p0threshold = 1.58655253931457074e-01 # Z = 1.0 sigma, Model independent 
    #p0threshold = 6.68072012688580852e-02 # Z = 1.5 sigma, Model independent 
    p0threshold = 3.59303191129258792e-02 # Z = 1.8 sigma, Model independent 
    #p0threshold = 2.27501319481792086e-02 # Z = 2.0 sigma, Model independent 

    smallexcess = {} #smallexcess[prod+"_"+lsp+"_"+ana+"_"+sr]=TTree ==> Should be smallexcesses[seid] = {"lspl": model, "prod": prod, "lsp": lsp, "ana": ana, "sr": sr, "entry": entry}
    for prod in production:
        for ana in production[prod]:
            for sr in production[prod][ana]:
                p0 = production[prod][ana][sr]
                if lsp != "bino" and prod == "EW" and sr == "EwkThreeLepton": continue
                if lsp != "bino" and prod == "QCD" and sr == "StopZ": continue
                ##JM##print "...  Processing LSP:", lsp, "production:", prod, "analysis:", ana, "SR:", sr, "p0 =", production[prod][ana][sr]
                # Naming of CLs_exp and CLs_obs branches
                if prod == "QCD" and ("ZeroLeptonStrong" in ana or "MonoJetF" in ana):
                    e_tr_CLs_exp = "("+prod+"_CLs_exp_"+ana+"_"+sr+"+(1+x"+prod+"_Truth_CLs_exp_"+ana+"_"+sr+")*("+prod+"_CLs_exp_"+ana+"_"+sr+"==-1)+2*(x"+prod+"_Truth_CLs_exp_"+ana+"_"+sr+"==-1))"
                    e_tr_CLs_obs = "("+prod+"_CLs_obs_"+ana+"_"+sr+"+(1+x"+prod+"_Truth_CLs_obs_"+ana+"_"+sr+")*("+prod+"_CLs_obs_"+ana+"_"+sr+"==-1)+2*(x"+prod+"_Truth_CLs_obs_"+ana+"_"+sr+"==-1))"
                elif (prod == "EW" or prod == "Mix") and ("ZeroLeptonStrong" in ana or "MonoJetF" in ana):
                    e_tr_CLs_exp = "(x"+prod+"_Truth_CLs_exp_"+ana+"_"+sr+")"
                    e_tr_CLs_obs = "(x"+prod+"_Truth_CLs_obs_"+ana+"_"+sr+")"
                else:
                    e_tr_CLs_exp = "("+prod+"_CLs_exp_"+ana+"_"+sr+")"
                    e_tr_CLs_obs = "("+prod+"_CLs_obs_"+ana+"_"+sr+")"
                    if len(production[prod][ana]) == 1:
                        e_tr_CLs_exp = "("+prod+"_CLs_exp_"+ana+")"
                        e_tr_CLs_obs = "("+prod+"_CLs_obs_"+ana+")"
                        
                # Prepare a list of lspls and Analysis+SR that pass the condition
                condition_se = "("+e_tr_CLs_exp+"<=0.05) && ("+e_tr_CLs_obs+">0.05) && ("+e_tr_CLs_exp+"!=-1) && ("+e_tr_CLs_obs+"!=-1)" # Small excesses, sensitivity and not excluded
                condition_ex = "("+e_tr_CLs_exp+"<=0.05) && ("+e_tr_CLs_obs+"<=0.05) && ("+e_tr_CLs_exp+"!=-1) && ("+e_tr_CLs_obs+"!=-1)" # Sensitivity and excluded
                #print condition
                # Small excesses for p0 < threshold
                if p0 <= p0threshold: 
                    setree = susy.CopyTree(condition_se)
                    # Loop on models in a tree
                    for entry in setree:
                        model = int(entry.model)
                        seid = str(model)+"_"+prod+"_"+lsp+"_"+ana+"_"+sr
                        smallexcess[seid] = {"model": model,
                                             "prod": prod,
                                             "lsp": lsp,
                                             "ana": ana,
                                             "sr": sr,
                                             "entry": copy.deepcopy(entry)}
                    print "--->",len(smallexcess), "entries have small excess in production:", prod, "LSP:", lsp, "analysis:", ana, "SR:", sr
    f.Close()
    #for sk in smallexcess: print smallexcess[sk].GetEntries()
    return smallexcess

def getInfoFromKey(key):
    info = key.split("_")
    prod = info[0]
    mode = info[1]
    ana = info[2]
    strippart = prod+"_"+mode+"_"+ana+"_"
    sr = key.replace(strippart,"")
    #print "JM",prod,mode,ana,sr
    return prod,mode,ana,sr

def inspectSmallExcesses(smallexcesses):
    if not os.path.isdir(plotsdir): os.system("mkdir "+plotsdir)
    print "... Inspecting models passing criteria."
    se = {}
    count = 0
    lmodels = []
    for key in smallexcesses:
        model = smallexcesses[key]["model"]
        #prod  = smallexcesses[key]["prod"]
        #lsp   = smallexcesses[key]["lsp"]
        #ana   = smallexcesses[key]["ana"]
        #sr    = smallexcesses[key]["sr"]
        #entry = smallexcesses[key]["entry"]
        lmodels.append(model)
    lmodels_no = set(lmodels)
    print "===> There are", len(lmodels), "small excess models found, in total", len(lmodels_no), "non overlapping." 
    if True:
        "... Models with small excess:"
        # Use histogram to find doubly found models
        hist = TH1I("models","Models not excluded by analysis with small excess",220000,0,220000)
        for i in smallexcesses:
            ##JM##print smallexcesses[i]["model"], smallexcesses[i]
            hist.Fill(smallexcesses[i]["model"])            
        cm = TCanvas()
        cm.cd()
        hist.Draw()
        cm.Print(plotsdir+"models.pdf")
        lmodels_once, lmodels_twice, lmodels_trice, lmodels_four, lmodels_five = [],[],[],[],[]
        # Inspect models found by more analyses, at least 2, at least 3
        for bin in range(1,hist.GetNbinsX()):
            model = int(hist.GetBinLowEdge(bin))
            if hist.GetBinContent(bin) == 1:
                #print "===> ", model, "found once."
                lmodels_once.append(model)
            if hist.GetBinContent(bin) == 2:
                #print "===> ", model, "found twice."
                #for key in se:
                #    if se[key]["model"] == model:
                #        print se[key]
                lmodels_twice.append(model)
            if hist.GetBinContent(bin) == 3:
                #print "===> ", model, "found three times."
                #for key in se:
                #    if se[key]["model"] == model:
                #        print se[key]
                lmodels_trice.append(model)
            if hist.GetBinContent(bin) == 4:
                #print "===> ", model, "found four times."
                #for key in se:
                #    if se[key]["model"] == model:
                #        print se[key]
                lmodels_four.append(model)
            if hist.GetBinContent(bin) >= 5:
                #print "===> ", model, "found five or more times."
                #for key in se:
                #    if se[key]["model"] == model:
                #        print se[key]
                lmodels_five.append(model)
        print "There were", len(lmodels_five), "models found five times,", len(lmodels_four), "models found four times,", len(lmodels_trice), "models found three times,", len(lmodels_twice), "found twice, and", len(lmodels_once), "found once."
    #return lmodels_once, lmodels_twice, lmodels_trice, lmodels_four, lmodels_five
    return

def inspectExclusion(smallexcesses):
    print "... Inspecting if models being excluded."
    se = {}
    count = 0
    lmodels = []
    for key in smallexcesses:
        model = smallexcesses[key]["model"]
        #prod  = smallexcesses[key]["prod"]
        #lsp   = smallexcesses[key]["lsp"]
        #ana   = smallexcesses[key]["ana"]
        #sr    = smallexcesses[key]["sr"]
        #tree  = smallexcesses[key]["entry"]
        lmodels.append(model)
    lmodels_no = set(lmodels)

    print "===> There are", len(lmodels), "small excess models found, in total", len(lmodels_no), "non overlapping."#, list(lmodels_no)

    # Now loop on original input, in the condition loop on each model
    print "... Exclusion for small excess models."
    if not os.path.isdir(plotsdir): os.system("mkdir "+plotsdir)
    excluded = {}
    lsps = ["bino","wino","higgsino"]
    prods = ["QCD","EW","Mix"]
    # Now loop on original input susy trees, but condition only models from the list
    count = 0
    for lsp in lsps:
        fname = fnames[lsp]
        f = openFile(fname)
        susy = f.Get("susy")        
        susy.SetBranchStatus("*Cat*",0)
        susy.SetBranchStatus("*LL*",0)
        susy.SetBranchStatus("*RecoEvents*",0)
        susy.SetBranchStatus("*Higgs*",0)
        for prod in prods:
            # Use production
            production = getProduction(prod)
            for ana in production[prod]:
                # Merge conditions per sr
                condition_ex = ""
                # Determine how many chunks for SR and models should be used
                # ncond = 4 * nSR + nmodels
                # nm = 50*M / (4*SR+M) sub-list of models
                # nsr = 50*SR / (4*SR+M) sub-list of SRs
                nSR = len(production[prod][ana])
                nM  = len(list(lmodels_no))
                nsr = 40 * nSR / (4*nSR + nM)
                nm  = 40 * nM / (4*nSR + nM)
                print ana, production[prod][ana].keys(), list(lmodels_no)
                subs_list_sr = np.array_split(production[prod][ana].keys(), int(math.ceil(float(len(production[prod][ana].keys()))/nsr)))
                subs_list_m  = np.array_split(list(lmodels_no), int(math.ceil(float(len(list(lmodels_no)))/nm)))
                print subs_list_sr, subs_list_m
                if nM == 0: continue 
                for sub_list_sr in subs_list_sr:
                    for sr in list(sub_list_sr):
                        p0 = production[prod][ana][sr]
                        if lsp != "bino" and prod == "EW" and sr == "EwkThreeLepton": continue
                        if lsp != "bino" and prod == "QCD" and sr == "StopZ": continue
                        ##JMprint "...  Processing LSP:", lsp, "production:", prod, "analysis:", ana, "SR:", sr, "p0 =", production[prod][ana][sr]
                        # Naming of CLs_exp and CLs_obs branches
                        if prod == "QCD" and ("ZeroLeptonStrong" in ana or "MonoJetF" in ana):
                            e_tr_CLs_exp = "("+prod+"_CLs_exp_"+ana+"_"+sr+"+(1+x"+prod+"_Truth_CLs_exp_"+ana+"_"+sr+")*("+prod+"_CLs_exp_"+ana+"_"+sr+"==-1)+2*(x"+prod+"_Truth_CLs_exp_"+ana+"_"+sr+"==-1))"
                            e_tr_CLs_obs = "("+prod+"_CLs_obs_"+ana+"_"+sr+"+(1+x"+prod+"_Truth_CLs_obs_"+ana+"_"+sr+")*("+prod+"_CLs_obs_"+ana+"_"+sr+"==-1)+2*(x"+prod+"_Truth_CLs_obs_"+ana+"_"+sr+"==-1))"
                        elif (prod == "EW" or prod == "Mix") and ("ZeroLeptonStrong" in ana or "MonoJetF" in ana):
                            e_tr_CLs_exp = "(x"+prod+"_Truth_CLs_exp_"+ana+"_"+sr+")"
                            e_tr_CLs_obs = "(x"+prod+"_Truth_CLs_obs_"+ana+"_"+sr+")"
                        else:
                            e_tr_CLs_exp = "("+prod+"_CLs_exp_"+ana+"_"+sr+")"
                            e_tr_CLs_obs = "("+prod+"_CLs_obs_"+ana+"_"+sr+")"
                            if len(production[prod][ana]) == 1:
                                e_tr_CLs_exp = "("+prod+"_CLs_exp_"+ana+")"
                                e_tr_CLs_obs = "("+prod+"_CLs_obs_"+ana+")"

                        # Prepare a list of models and Analysis+SR that pass the condition of not being excluded
                        condition_ex += "(("+e_tr_CLs_exp+"<=0.05) && ("+e_tr_CLs_obs+"<=0.05) && ("+e_tr_CLs_exp+"!=-1) && ("+e_tr_CLs_obs+"!=-1)) || " # Sensitivity and excluded
                    
                        # Loop on models in chunks of 20
                        # Split list of models into chunks of nm
                        #subs_lmodels_no = np.array_split(list(lmodels_no), int(math.ceil(float(len(list(lmodels_no)))/nm)))
                        for sub_list_m in subs_list_m:
                            cond_mo = ""
                            condition = ""
                            for mod in list(sub_list_m):
                                cond_mo += "model=="+str(int(mod))+" || "
                            condition = "("+condition_ex[:-3]+") && ("+cond_mo[:-3]+")"
                            print condition
                            extree = susy.CopyTree(condition)
                            excluded[count] = copy.deepcopy(extree)
                            count = count+1

    # Loop on not excluded models
    for key in excluded:
        tree  = excluded[key]
        for entry in tree:
            model = entry.model
            lexmodels.append(model)
    lexmodels_no = set(lexmodels)
    print "===> There are", len(lexmodels), "excluded, in total", len(lexmodels_no), "non overlapping:", list(lexmodels_no)
    
    if True:
        "... Excluded models with small excess:"
        # Use histogram to find doubly found models
        hist = TH1I("excluded","excluded",220000,0,220000)
        for model in list(lexmodels_no):
            hist.Fill(model)            
        cm = TCanvas()
        cm.cd()
        cm.SetLogy()
        hist.SetMinimum(0.1)
        hist.Draw()
        cm.Print(plotsdir+"excluded.pdf")

        # Inspect models found by more analyses, at least 2, at least 3
        exmodels = []
        for bin in range(1,hist.GetNbinsX()):
            model = int(hist.GetBinLowEdge(bin))
            if hist.GetBinContent(bin) > 0:
                exmodels.append(model)
                print model, "excluded", int(hist.GetBinContent(bin)), "times."
        print "Models", exmodels, "were excluded at least once."
    hope = (list(set(lmodels_no)-set(lexmodels_no)))
    print "Models", hope, "not excluded by any analysis."

    return 
        

def findModelExclusion(fname, mode, prod, model): # model id, prod: QCD, EW, Mix, mode: bino, wino, higgsino
    print "... Exclusion summary for model", model
    if not os.path.isdir(plotsdir): os.system("mkdir "+plotsdir)
    # Use production
    production = getProduction(prod)
    f = openFile(fname)
    susy = f.Get("susy")        
    susy.SetBranchStatus("*Cat*",0)
    susy.SetBranchStatus("*LL*",0)
    susy.SetBranchStatus("*RecoEvents*",0)
    susy.SetBranchStatus("*Higgs*",0)
    
    # Check if models excluded by other analyses, reverse findSmallExcesses
    # Take sub tree for model
    modeltree = susy.CopyTree("model=="+str(model))
    if modeltree.GetEntries() < 1: print "... WARNING: No info for model", model, "in the file", fname
    elif modeltree.GetEntries() > 1: print "... WARNING: Duplicated info for model", model, "in the file", fname
    # Loop on all input files, check if excluded by any other analysis
    smallexcess = {}
    models_smallexcess = {}
    models_excluded = {}
    for prod in production:
        for ana in production[prod]:
            for sr in production[prod][ana]:
                p0 = production[prod][ana][sr]
                if mode != "bino" and prod == "EW" and sr == "EwkThreeLepton": continue
                if mode != "bino" and prod == "QCD" and sr == "StopZ": continue
                #print "...  Processing LSP:", mode, "production:", prod, "analysis:", ana, "SR:", sr, "p0 =", production[prod][ana][sr]
                # Naming of CLs_exp and CLs_obs branches
                if prod == "QCD" and ("ZeroLeptonStrong" in ana or "MonoJetF" in ana):
                    e_tr_CLs_exp = "("+prod+"_CLs_exp_"+ana+"_"+sr+"+(1+x"+prod+"_Truth_CLs_exp_"+ana+"_"+sr+")*("+prod+"_CLs_exp_"+ana+"_"+sr+"==-1)+2*(x"+prod+"_Truth_CLs_exp_"+ana+"_"+sr+"==-1))"
                    e_tr_CLs_obs = "("+prod+"_CLs_obs_"+ana+"_"+sr+"+(1+x"+prod+"_Truth_CLs_obs_"+ana+"_"+sr+")*("+prod+"_CLs_obs_"+ana+"_"+sr+"==-1)+2*(x"+prod+"_Truth_CLs_obs_"+ana+"_"+sr+"==-1))"
                elif (prod == "EW" or prod == "Mix") and ("ZeroLeptonStrong" in ana or "MonoJetF" in ana):
                    e_tr_CLs_exp = "(x"+prod+"_Truth_CLs_exp_"+ana+"_"+sr+")"
                    e_tr_CLs_obs = "(x"+prod+"_Truth_CLs_obs_"+ana+"_"+sr+")"
                else:
                    e_tr_CLs_exp = "("+prod+"_CLs_exp_"+ana+"_"+sr+")"
                    e_tr_CLs_obs = "("+prod+"_CLs_obs_"+ana+"_"+sr+")"
                    if len(production[prod][ana]) == 1:
                        e_tr_CLs_exp = "("+prod+"_CLs_exp_"+ana+")"
                        e_tr_CLs_obs = "("+prod+"_CLs_obs_"+ana+")"
                #print "--->",e_tr_CLs_exp,e_tr_CLs_obs
                
                condition_ex = "("+e_tr_CLs_exp+"<=0.05) && ("+e_tr_CLs_obs+"<=0.05) && ("+e_tr_CLs_exp+"!=-1) && ("+e_tr_CLs_obs+"!=-1)" # Sensitivity and excluded
                condition_se = "("+e_tr_CLs_exp+"<=0.05) && ("+e_tr_CLs_obs+">0.05) && ("+e_tr_CLs_exp+"!=-1) && ("+e_tr_CLs_obs+"!=-1)" # Small excesses, sensitivity and not excluded
                
                # Loop on tree
                #for i in range(1,entries+1):
                #    modeltree.Scan("model:"+e_tr_CLs_exp+":"+e_tr_CLs_obs)

                # Exclusion
                modeltree.Draw("model:"+e_tr_CLs_exp+":"+e_tr_CLs_obs,condition_ex,"goff")
                for i in range(modeltree.GetSelectedRows()):
                    testmodel = modeltree.GetV1()[i]
                    CLs_exp, CLs_obs = modeltree.GetV2()[i], modeltree.GetV3()[i]
                    print "Excluded:", modeltree.GetV1()[i], modeltree.GetV2()[i], modeltree.GetV3()[i], "LSP:", mode, "production:", prod, "analysis:", ana, "SR:", sr
                    models_excluded[int(modeltree.GetV2()[i])] = {"prod": prod, "mode": mode, "ana": ana, "sr": sr, "CLs_exp": modeltree.GetV2()[i],  "CLs_exp": modeltree.GetV3()[i]}
                    
                # Excess
                modeltree.Draw("model:"+e_tr_CLs_exp+":"+e_tr_CLs_obs,condition_se,"goff")
                for i in range(modeltree.GetSelectedRows()):
                    testmodel = modeltree.GetV1()[i]
                    CLs_exp, CLs_obs = modeltree.GetV2()[i], modeltree.GetV3()[i]
                    print "Excess:", modeltree.GetV1()[i], modeltree.GetV2()[i], modeltree.GetV3()[i], "LSP:", mode, "production:", prod, "analysis:", ana, "SR:", sr
                    models_smallexcess[int(modeltree.GetV2()[i])] = {"prod": prod, "mode": mode, "ana": ana, "sr": sr, "CLs_exp": modeltree.GetV2()[i],  "CLs_exp": modeltree.GetV3()[i]}
                    
    # Plot their paramters (migrate 1D, 2D and 3D plots here)
    dumpJSON(models_excluded,"excluded")
    dumpJSON(models_smallexcess,"smallexcess")
                        
def openFile(fname):
    print "... Oppening file", fname
    return TFile.Open(fname)

def dumpJSON(result,mode):
    plotsdir = "plots/"
    if not os.path.isdir(plotsdir): os.system("mkdir "+plotsdir)
    with open(plotsdir+"_"+mode+"_pmssm.json", "w") as outfile:
        json.dump(result, outfile)

def main(argv):
    gStyle.SetOptStat(0)
    mode = "bino" # bino, wino, higgsino
    prod = "test" # QCD, EW, Mix
    inputdir = "input"
    doAll = False
    try:
       opts, args = getopt.getopt(argv,"hi:m:p:a",["--input","--mode","--prod"])
    except getopt.GetoptError:
       print 'python small_excesses_in_pmssm.py -i <inputdir> -m <mode:bino,wino,higgsino> -p <QCD,EW,Mix,test> -a' 
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'python small_excesses_in_pmssm.py -i <inputdir> -m <mode:bino,wino,higgsino> -p <QCD,EW,Mix,test> -a' 
          sys.exit()
       elif opt in ("-i"):
          inputdir = arg
       elif opt in ("-m"):
          mode = arg
       elif opt in ("-p"):
          prod = arg
       elif opt in ("-a"):
          doAll = True

    if False: listAllBranches(fnames[mode])
    if doAll:
        
        modes = ["bino","wino","higgsino"]
        prods = ["QCD","EW","Mix"]
        smallexcesses = {}
        excludeds     = {}
        for m in modes:    
            for p in prods:    
                se = findSmallExcesses(fnames[m], m, p)
                smallexcesses.update(se)
                # Explore exclusion for each mode and prod (duplicating a bit the scan, but no crash)
        inspectSmallExcesses(smallexcesses)
        inspectExclusion(smallexcesses) 
    else:
        smallexcesses = findSmallExcesses(fnames[mode], mode, prod)
        #inspectSmallExcesses(smallexcesses)
        #inspectExclusion(smallexcesses)

        
if __name__ == "__main__":
    main(sys.argv[1:])



# Explore in root with
# skim = susy->CopyTree("model==156814")


#skim->Scan("model:mdR:muR:msR:mcR:mbR:mtR:mqL1:mqL2:mqL3:meR:mmuR:mtauR:meL:mmuL:mtauL:M_1:M_2:M_3")
#skim->Scan("mu:A_b:A_t:A_tau:mA:tanb:m_h:m_W:m_b:m_d_L:m_d_R:m_u_L:m_u_R:m_s_L:m_s_R:m_c_L:m_c_R")
#skim->Scan("m_b_1:m_b_2:m_t_1:m_t_2:m_gl:m_e_L:m_e_R:m_mu_L:m_mu_R:m_tau_1:m_tau_2:m_nu_eL:m_nu_muL:m_nu_tauL")
#skim->Scan("m_chi_10:m_chi_20:m_chi_30:m_chi_40:m_chi_1p:m_chi_2p:m_A:m_H:m_Hp")
#skim->Scan("Cross_section_EW:Cross_section_QCD:Cross_section_Mix")
