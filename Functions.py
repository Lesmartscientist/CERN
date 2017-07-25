#!/usr/bin/python

from ROOT import *
import math
import operator
import itertools
import sys 
from Dictionaries import *

#Writing objects

class Photon:

    def __init__(self, mass, eta, phi, pt):
	    self.mass = mass
	    self.eta = eta
	    self.phi = phi
	    self.pt  = pt
                
class Jet:

    def __init__(self,mass,eta,phi,pt):
         self.mass = mass
         self.eta = eta
         self.phi = phi
         self.pt = pt

def recojet(mass , eta, phi, pt):
    vectJ = TLorentzVector(0.,0.,0.,0.)
    vectJ.SetPtEtaPhiM(pt, eta, phi, mass)

    return vectJ;

def getCompPhoton(vect1, vect2):

    cVect = TLorentzVector(0.,0.,0.,0.)
    cVect = vect1+vect2

    return cVect;

def getCompJet(jvect1 ,jvect2):
    
    cjVect = TLorentzVector(0.,0.,0.,0.)
    cjVect = jvect1+jvect2

    return cjVect;
def getCompJetReweight(jvect1, jvect2):
    cjVect = TLorentzVector(0.,0.,0.,0.)
    cjVect = jvect1+jvect2
    cjVect *= 125.0/cjVect.M()
    return cjVect;

def recoDiH(hvect1, hvect2):

    divect = TLorentzVector(0.,0.,0.,0.)
    divect = hvect1+hvect2
    
    return divect;

def deltaStuff(vect1, vect2):
    dic = {}
    dic["dphi"] = abs(vect1.DeltaPhi(vect2))
    dic["dR"] = vect1.DeltaR(vect2)
    dic["dEta"] = abs(vect1.Eta()-vect2.Eta())
    dic["dpt"] = abs(vect1.Pt()-vect2.Pt())

    return dic;       

def Jets125(jlist):
    Jet125 = [] 
    for JetPair in itertools.combinations(jlist,r=2):
        cjets = getCompJet(*JetPair)
        Jet125.append(cjets)
    d = {}
    for RecoStuff in Jet125:
        d[abs(RecoStuff.M()-125.0)] = RecoStuff
    s = sorted(d.items())
    HiggsCan = s[0]
    HiggsCan2 = HiggsCan[1]

    return HiggsCan2;

def fill(objectlist,pthistos,etahistos,phihistos,mhistos):
    for key in objectlist:
        pthistos[key].Fill(objectlist[key].Pt())
        etahistos[key].Fill(objectlist[key].Eta())
        phihistos[key].Fill(objectlist[key].Phi())  
        mhistos[key].Fill(objectlist[key].M())

def getdiff(diffdictionary,dphihistos,detahistos,dpthistos,dRhistos):
    for key in diffdictionary:        
        k = deltaStuff(diffdictionary[key][0],diffdictionary[key][1])
        dphihistos[key].Fill(k.get("dphi"))
        detahistos[key].Fill(k.get("dEta"))
        dpthistos[key].Fill(k.get("dpt"))
        dRhistos[key].Fill(k.get("dR"))
        k = 0

def matrix(yybbvector):

    transmatrix = TMatrix(4,4)
    Qsquared = abs(yybbvector.Dot(yybbvector))
    ptsquared = yybbvector.Pt()*yybbvector.Pt()
    Xtrans = math.sqrt(Qsquared + ptsquared)
    transmatrix[0][0] = -1*(yybbvector.Pt()/math.sqrt(Qsquared))
    transmatrix[0][1] = 0
    transmatrix[0][2] = -1*(yybbvector.E()/math.sqrt(Qsquared))
    transmatrix[0][3] = yybbvector.Px()/math.sqrt(Qsquared)
    transmatrix[1][0] = Xtrans/math.sqrt(Qsquared)
    transmatrix[1][1] = 0
    transmatrix[1][2] = (yybbvector.Pt()*yybbvector.Px())/(math.sqrt(Qsquared)*Xtrans)
    transmatrix[1][3] = -1*(yybbvector.Pt()*yybbvector.Px())/(math.sqrt(Qsquared) *Xtrans)
    transmatrix[2][0] = 0
    transmatrix[2][1] = 1
    transmatrix[2][2] = 0
    transmatrix[2][3] = 0
    transmatrix[3][0] = 0
    transmatrix[3][1] = 0
    transmatrix[3][2] = yybbvector.Px()/Xtrans
    transmatrix[3][3] = -1*(yybbvector.E()/Xtrans)
#    print "Q^2",Qsquared
#    print "P_t",ptsquared
#    print "X_t",Xtrans
#    print "element 00",transmatrix[0][0]
#    print "element 02",transmatrix[0][2]
#    print "element 03",transmatrix[0][3]
#    print "element 10",transmatrix[1][0]
#    print "element 12",transmatrix[1][2]
#    print "element 13",transmatrix[1][3]
#    print "element 32",transmatrix[3][2]
#    print "element 33",transmatrix[3][3]

    return transmatrix;

def matrixmult(yybbmatrix, lorentzvect):
	colinsoper = TMatrix(4,1)
	colinsopertlorentz = TLorentzVector(0.,0.,0.,0.)
	labframe = TMatrix(4,1)
	labframe[0][0] = lorentzvect.Px()
	labframe[1][0] = lorentzvect.Py()
	labframe[2][0] = lorentzvect.Pz()
	labframe[3][0] = lorentzvect.Energy()
	colinsoper.Mult(yybbmatrix,labframe)
	colinsopertlorentz.SetPx(colinsoper[0][0])
	colinsopertlorentz.SetPy(colinsoper[1][0])
	colinsopertlorentz.SetPz(colinsoper[2][0])
	colinsopertlorentz.SetE(colinsoper[3][0])

	return colinsopertlorentz;

def colinsopercombo(colinsopervect1,colinsopervect2):
	recovector = TLorentzVector(0.,0.,0.,0.)
	recovector = colinsopervect1+colinsopervect2
	
	return recovector;

def fillcolinsoper(colinsopervect1,colinsopervect2,plot, tree):
	deltaphi = colinsopervect1.DeltaPhi(colinsopervect2)
	weight = tree.hgam_weight*tree.yybb_high_weight*tree.lumiXsecWeight
	plot.Fill(deltaphi,weight)

        
def analyze(infileName, yybb_bcat_arg,listofobjectname,differentiallist,pthistos,etahistos,phihistos,dpthistos,dphihistos,detahistos,dRhistos,mhistos,cans):
    gROOT.SetBatch(True)
#    infile = TFile("/atlas/data/caleb/lxplus-home/cmosakow/NIU_Data/tiny/"+infileName+".root","READ")
    infile = TFile("/atlas/data/caleb/lxplus-home/cmosakow/NIU_Data/v4/muon_corrections/"+infileName+".root")
    infile.ls()
    tree = infile.Get("mini")
    HPTH = TH1F("HPTH","Highest Pt Photon",100,25,350)
    SLPTH = TH1F("SLPTH", "Subleading Pt Photon",100,15,190)
    JetLH = TH1F("JLH","Jet Leading Pt",100,20,250)
    JetSLH = TH1F("JSLH","Jet Sub Leading Pt",100,0,100)
    JSSLH = TH1F("JSSLH", "3rd Jet Leading Pt", 100,15, 100)
    JSSSLH = TH1F("JSSSLH", "4th Jet Leading Pt",100,15,100)
    Jet_nHist = TH1F("Jet_n","Jet Multiplicity",200,0,20)
    Jet_ptHist = TH1F("Jet_pt","Jet Pt",100,15,200)
    Jet_etaHist = TH1F("Jet_eta", "Jet Eta",100,-5,5)
    Jet_phiHist = TH1F("Jet_phi","Jet Phi",100,-3.5,3.5)
    Jet_mHist = TH1F("Jet_m","Jet Mass",100,0,50)
    Gam_nHist = TH1F("Gam_n","Photon Multiplicity",100,0,4)
    Gam_ptHist = TH1F("Gam_pt","Photon Pt",100,0,250)
    Gam_etaHist = TH1F("Gam_eta","Photon Eta",100,-3,3)
    Gam_phiHist = TH1F("Gam_phi","Photon Phi",100,-3.5,3.5)
    Gam_mHist = TH1F("Gam_m","Photon Mass",100,-3,3)
    mgammagammaHist = TH1F("m(gamma,gamma)", "m(gamma,gamma)",100,100, 180)
    mjjHist = TH1F("mjj", "mjj", 100, 0,400)
    DHH = TH1F("diHm", "diHm",100,0,1000)
    mjj125H = TH1F("mjj125", "", 100, 0, 200)
    DHCMH = TH1F("DHCM", "",100,0,1000)
    dphiyyH = TH1F("dphiyyH", "",30,0,5)
    dEtayyH = TH1F("dEtayyH", "",30,0,5)
    dRyyH = TH1F("dRyyH", "",30,0,5)
    dphijjH = TH1F("dphijjH", "",30,0,5)
    dEtajjH = TH1F("dEtajjH", "",30,0,5)
    dRjjH = TH1F("dRjjH","",30,0,5)
    dphiyyjjH = TH1F("dphiyyjjH", "",30,0,5)
    dEtayyjjH = TH1F("dEtayyjjH", "",30,0,5)
    dRyyjjH = TH1F("dRyyjjH", "",30,0,5)
    dphiyjHPTBH = TH1F("dphiyjHPTBH","",30,0,5)
    dEtayjHPTBH = TH1F("dEtayjHPTBH","",30,0,5)
    dRyjHPTBH = TH1F("dRyjHPTBH","",30,0,5)
    dRyjHPTPH = TH1F("dRyjHPTPH","",30,0,5)
    dphiyjHPTPH = TH1F("dphiyjHPTPH","",30,0,5)
    dEtayjHPTPH = TH1F("dEtayjHPTPH","",30,0,5)
    dRyjHPTJH = TH1F("dRyjHPTJH","",30,0,5)
    dEtayjHPTJH = TH1F("dEtayjHPTJH","",30,0,5)
    dphiyjHPTJH = TH1F("dphiyjHPTJH","",30,0,5)
    dphiyjSPTBH = TH1F("dphiyjSPTBH","",30,0,5)
    dEtayjSPTBH = TH1F("dEtayjSPTBH","",30,0,5)
    dRyjSPTBH = TH1F("dRyjSPTBH","",30,0,5)
    dphiyjjHPTPH = TH1F("dphiyjjHPTPH","",30,0,5)
    dEtayjjHPTPH = TH1F("dEtayjjHPTPH","",30,0,5)
    dRyjjHPTPH = TH1F("dRyjjHPTPH","",30,0,5)
    dphiyjjSPTPH = TH1F("dphiyjjSPTPH","",30,0,5)
    dEtayjjSPTPH = TH1F("dEtayjjSPTPH","",30,0,5)
    dRyjjSPTPH = TH1F("dRyjjSPTPH","",30,0,5)
    dRyyjHPTJH = TH1F("dRyyjHPTJH","",30,0,5)
    dphiyyjHPTJH = TH1F("dphiyyjHPTJH","",30,0,5)
    dEtayyjHPTJH = TH1F("dEtayyjHPTJH","",30,0,5)
    dphiyyjSPTJH = TH1F("dphiyyjSPTJH","",30,0,5)
    dEtayyjSPTJH = TH1F("dEtayyjSPTJH","",30,0,5)
    dRyyjSPTJH = TH1F("dRyyjSPTJH","",30,0,5)
    dEtayyj_yyjjH = TH2F("dEtayyj_yyjj","",30,0,5,30,0,5)
    dEtayyj_yyjSH = TH2F("dEtayyj_yyjSH","",30,0,5,30,0,5)
    dEtayyjj_yyjSH = TH2F("dEtayyjj_yyjSH","",30,0,5,30,0,5)
    dEtayyj_yyjjH.GetXaxis().SetTitle("dEta(yy,j) High Pt Jet")
    dEtayyj_yyjjH.GetYaxis().SetTitle("dEta(yy,jj)")
    dEtayyj_yyjSH.GetXaxis().SetTitle("dEta(yy,j) Leading Pt Jet")
    dEtayyj_yyjSH.GetYaxis().SetTitle("dEta(yy,j) Sub Pt Jet")
    dEtayyjj_yyjSH.GetXaxis().SetTitle("dEta(yy,jj)")
    dEtayyjj_yyjSH.GetYaxis().SetTitle("dEta(yy,j) Sub Pt Jet")
    deltaphi_colinsoper = TH1F("deltaphi_colinsoper","",30,0,4)
    deltaphi_colinsoper.GetXaxis().SetTitle("Phi")
    deltaphi_colinsoper.GetYaxis().SetTitle("Events")
    momentum_Higgs = TH1F("deltamomentum","",40,0,1000)
    test_CSyybb_p = TH1F("CollinsSoperyybb_p","",40,0,1000)
    test_CSyybb_m = TH1F("CollinsSoperyybb_m","",40,0,1000)
    for hist in listofobjectname:
        pthistos[hist] = TH1F(hist+"_pt",hist+"_pt",30,0,300)
        pthistos[hist].GetXaxis().SetTitle("pt [GeV]")
        pthistos[hist].GetYaxis().SetTitle("Events")
        etahistos[hist] = TH1F(hist+"_eta",hist+"_eta",20,-5,5)
        etahistos[hist].GetXaxis().SetTitle("Eta")
        etahistos[hist].GetYaxis().SetTitle("Events")
        phihistos[hist] = TH1F(hist+"_phi",hist+"_phi",20,-4,4)
        phihistos[hist].GetXaxis().SetTitle("Phi")
        phihistos[hist].GetYaxis().SetTitle("Events")
        mhistos[hist] = TH1F(hist+"_m",hist+"_m",100,0,500)
    for hist in differentiallist:
        dpthistos[hist] = TH1F(hist+"_d(pt)",hist+"_d(pt)",30,0,100)
        dpthistos[hist].GetXaxis().SetTitle("d(pt) [GeV]")
        dpthistos[hist].GetYaxis().SetTitle("Events")
        detahistos[hist] = TH1F(hist+"_d(eta)",hist+"_d(eta)",30,0,5)
        detahistos[hist].GetXaxis().SetTitle("d(Eta)")
        detahistos[hist].GetYaxis().SetTitle("Events")
        dphihistos[hist] = TH1F(hist+"_d(phi)",hist+"_d(phi)",30,0,5)
        dphihistos[hist].GetXaxis().SetTitle("d(Phi)")
        dphihistos[hist].GetYaxis().SetTitle("Events")
        dRhistos[hist] = TH1F(hist+"_dR",hist+"_dR",30,0,5)
        dRhistos[hist].GetXaxis().SetTitle("dR")
        dRhistos[hist].GetYaxis().SetTitle("Events")
    for hist in listofobjectname:
        cans[hist] = TCanvas(hist,hist,800,600)
    for hist in differentiallist:
        gROOT.SetBatch(True)
        cans[hist] = TCanvas(hist,hist,800,600)
#Get NumOfEntries
    NumOfEntries = tree.GetEntries()
#Loop
    NumOfEventsPassed = 0
    wrong = 0
    for i in xrange(NumOfEntries):
        tree.GetEntry(i)
        if (i%10000) == 0:
            print "Mismatched: "+str(wrong)
        if (tree.hgam_isPassed == False or tree.yybb_high_btagCat != yybb_bcat_arg or tree.yybb_cutList_high<5 or tree.photon_n != 2 or tree.jet_n <2 or tree.photon_isTight[0]==0 or tree.photon_isTight[1]==0 or tree.photon_iso_Loose[0]==0 or tree.photon_iso_Loose[1]==0 or (tree.photon_pt[0]/tree.m_yy)<0.35 or (tree.photon_pt[1]/tree.m_yy)<0.25):
            continue
        else:
            Jetlist = []
            NumOfEventsPassed +=1
            Jet_nHist.Fill(tree.jet_n)
            for j in xrange(tree.jet_n):
                jet_r = recojet(tree.jet_m_allMu[j], tree.jet_eta_allMu[j], tree.jet_phi_allMu[j], tree.jet_pt_allMu[j]) 
                Jetlist.append(jet_r) 
                Jet_ptHist.Fill(tree.jet_pt_allMu[j])
                Jet_phiHist.Fill(tree.jet_phi_allMu[j])
                Jet_mHist.Fill(tree.jet_m_allMu[j])
                Jet_etaHist.Fill(tree.jet_eta_allMu[j])
#            Jetlist.sort(key= operator.attrgetter('Pt'),reverse=True)
            if Jetlist[0].Pt()<Jetlist[1].Pt():
                jet1 = Jetlist.pop(1)
                jet2 = Jetlist.pop(0)
                Jetlist.append(jet1)
                Jetlist.append(jet2)
                wrong +=1
#                print "Event Number:",tree.eventNumber
#            print "Jet Pt's: ", Jetlist[0].Pt(),Jetlist[1].Pt()
            HiggsCan = Jets125(Jetlist)
            DHCMweight = tree.hgam_weight*tree.yybb_low_weight*tree.lumiXsecWeight
            HiggsCanM = HiggsCan.M()
            mjj125H.Fill(HiggsCanM)
            dStuff1 = deltaStuff(Jetlist[0],Jetlist[1])
            dRjjH.Fill(dStuff1.get('dR'))
            dEtajjH.Fill(dStuff1.get('dEta'))
            dphijjH.Fill(dStuff1.get('dphi'))
            if len(Jetlist) == 4:
                JSSSLH.Fill(Jetlist[3].Pt())
            elif len(Jetlist) == 3:
                JSSLH.Fill(Jetlist[2].Pt())
            else:
                JetLH.Fill(Jetlist[0].Pt())
                JetSLH.Fill(Jetlist[1].Pt())
            Higgs2 = getCompJet(Jetlist[0], Jetlist[1])
            Higgs2_m = Higgs2.M()
            mjjHist.Fill(Higgs2_m)
            photonlist = []
            Gam_nHist.Fill(tree.photon_n)
            for j in xrange(tree.photon_n):
                photon_r = recojet(tree.photon_m[j], tree.photon_eta[j], tree.photon_phi[j], tree.photon_pt[j])
                photonlist.append(photon_r) 
                Gam_ptHist.Fill(tree.photon_pt[j])
                Gam_etaHist.Fill(tree.photon_eta[j])
                Gam_phiHist.Fill(tree.photon_phi[j])
                Gam_mHist.Fill(tree.photon_m[j])
            photonlist.sort(key= operator.attrgetter('Pt'),reverse=True)
            dStuff2 = deltaStuff(photonlist[0], photonlist[1])
            dRyyH.Fill(dStuff2.get('dR'))
            dEtayyH.Fill(dStuff2.get('dEta'))
            dphiyyH.Fill(dStuff2.get('dphi'))
            Higgs = getCompPhoton(photonlist[0], photonlist[1])
            Higgs_m = Higgs.M()
            mgammagammaHist.Fill(Higgs_m)
            HPTH.Fill(photonlist[0].Pt())
            SLPTH.Fill(photonlist[1].Pt())
            dStuff3 = deltaStuff(Higgs,Higgs2)
            dphiyyjjH.Fill(dStuff3.get('dphi'))
            dEtayyjjH.Fill(dStuff3.get('dEta'))
            dRyyjjH.Fill(dStuff3.get('dR'))	
            dStuff4 = deltaStuff(photonlist[0],Jetlist[0])
            dRyjHPTBH.Fill(dStuff4.get('dR'))
            dEtayjHPTBH.Fill(dStuff4.get('dEta'))
            dphiyjHPTBH.Fill(dStuff4.get('dphi'))
            dStuff5 = deltaStuff(photonlist[0],Jetlist[1])
            dRyjHPTPH.Fill(dStuff5.get('dR'))
            dEtayjHPTPH.Fill(dStuff5.get('dEta'))
            dphiyjHPTPH.Fill(dStuff5.get('dphi'))
            dStuff6 = deltaStuff(photonlist[1],Jetlist[0])
            dRyjHPTJH.Fill(dStuff6.get('dR'))
            dEtayjHPTJH.Fill(dStuff6.get('dEta'))
            dphiyjHPTJH.Fill(dStuff6.get('dphi'))
            dStuff7 = deltaStuff(photonlist[1],Jetlist[1])
            dRyjSPTBH.Fill(dStuff7.get('dR'))
            dEtayjSPTBH.Fill(dStuff7.get('dEta'))
            dphiyjSPTBH.Fill(dStuff7.get('dphi'))
            dStuff8 = deltaStuff(photonlist[0],Higgs2)
            dRyjjHPTPH.Fill(dStuff8.get('dR'))
            dEtayjjHPTPH.Fill(dStuff8.get('dEta'))
            dphiyjjHPTPH.Fill(dStuff8.get('dphi'))
            dStuff9 = deltaStuff(photonlist[1],Higgs2)
            dRyjjSPTPH.Fill(dStuff9.get('dR'))
            dEtayjjSPTPH.Fill(dStuff9.get('dEta'))
            dphiyjjSPTPH.Fill(dStuff9.get('dphi'))
            dStuff10 = deltaStuff(Higgs,Jetlist[0])
            dRyyjHPTJH.Fill(dStuff10.get('dR'))
            dEtayyjHPTJH.Fill(dStuff10.get('dEta'))
            dphiyyjHPTJH.Fill(dStuff10.get('dphi'))
            dStuff11 = deltaStuff(Higgs,Jetlist[1])
            dRyyjSPTJH.Fill(dStuff11.get('dR'))
            dEtayyjSPTJH.Fill(dStuff11.get('dEta'))
            dphiyyjSPTJH.Fill(dStuff11.get('dphi'))
            diHiggs = recoDiH(Higgs, Higgs2)
            diHm = diHiggs.M()
            DHH.Fill(diHm)
            DHC = recoDiH(Higgs, HiggsCan)
            DHCM = DHC.M()
            DHCMH.Fill(DHCM,DHCMweight)
            dEtayyj_yyjjH.Fill(dStuff10.get('dEta'),dStuff3.get('dEta'))
            dEtayyj_yyjSH.Fill(dStuff10.get('dEta'),dStuff11.get('dEta'))
            dEtayyjj_yyjSH.Fill(dStuff3.get('dEta'),dStuff11.get('dEta'))
            y1j1 = getCompJet(photonlist[0],Jetlist[0])
            y1j2 = getCompJet(photonlist[0],Jetlist[1])
            y2j1 = getCompJet(photonlist[1],Jetlist[0])
            y2j2 = getCompJet(photonlist[1],Jetlist[1])
            yyj1 = getCompJet(Higgs,Jetlist[0])
            yyj2 = getCompJet(Higgs, Jetlist[1])
            y1jj = getCompJet(photonlist[0],Higgs2)
            y2jj = getCompJet(photonlist[1], Higgs2)
            dictofobjects = {}
            dictofobjects["y1j1"] = y1j1
            dictofobjects["y1j2"] = y1j2
            dictofobjects["y2j1"] = y2j1
            dictofobjects["y2j2"] = y2j2
            dictofobjects["yyj1"] = yyj1
            dictofobjects["yyj2"] = yyj2
            dictofobjects["y1jj"] = y1jj
            dictofobjects["y2jj"] = y2jj
            dictofobjects["yyjj"] = diHiggs
            dictofobjects["yy"] = Higgs
            dictofobjects["jj"] = Higgs2
            differential = {}
            differential["yyj1_j2"] = [yyj1,Jetlist[1]]
            differential["yyj2_j1"] = [yyj2,Jetlist[0]]
            differential["y1jj_y2"] = [y1jj,photonlist[1]]
            differential["y2jj_y1"] = [y2jj,photonlist[0]]
            differential["yy_jj"] = [Higgs,Higgs2]
            differential["y1j1_y2j2"] = [y1j1,y2j2]
            differential["y2j1_y1j2"] = [y2j1,y1j2]
            differential["y1_j1"] = [photonlist[0],Jetlist[0]]
            differential["y1_j2"] = [photonlist[0],Jetlist[1]]
            differential["y2_j1"] = [photonlist[1],Jetlist[0]]
            differential["y2_j2"] = [photonlist[1],Jetlist[1]]
            differential["y1_y2"] = [photonlist[0],photonlist[1]]
            differential["j1_j2"] = [Jetlist[0],Jetlist[1]]
            differential["yy_j1"] = [Higgs,Jetlist[0]]
            differential["yy_j2"] = [Higgs,Jetlist[1]]
            differential["jj_y1"] = [Higgs2,photonlist[0]]
            differential["jj_y2"] = [Higgs2,photonlist[1]]
            differential["y1j1_j2"] = [y1j1,Jetlist[1]]
            differential["y1j1_y2"] = [y1j1,photonlist[1]]
            differential["y1j2_y2"] = [y1j2,photonlist[1]]
            differential["y1j2_j1"] = [y1j2,Jetlist[0]]
            differential["y2j1_j2"] = [y2j1, Jetlist[1]]
            differential["y2j1_y1"] = [y2j1, photonlist[0]]
            differential["y2j2_y1"] = [y2j2, photonlist[0]]
            differential["y2j2_j1"] = [y2j2, Jetlist[0]]
            getdiff(differential,dphihistos,detahistos,dpthistos,dRhistos)
            fill(dictofobjects,pthistos,etahistos,phihistos,mhistos)
            transmatrix = matrix(diHiggs)
            CSgam1 = matrixmult(transmatrix, photonlist[0])
            CSgam2 = matrixmult(transmatrix, photonlist[1])
            CSjet1 = matrixmult(transmatrix, Jetlist[0])
            CSjet2 = matrixmult(transmatrix, Jetlist[1])
            CSyyHiggs = colinsopercombo(CSgam1, CSgam2)
            CSbbHiggs = colinsopercombo(CSjet1, CSjet2)
#            CSdiHiggs = matrixmult(transmatrix,diHiggs)
#            test_CSyybb_p.Fill(CSdiHiggs.P())
#            test_CSyybb_m.Fill(CSdiHiggs.M())
            fillcolinsoper(CSyyHiggs, CSbbHiggs, deltaphi_colinsoper, tree)
#            momentum_Higgs.Fill(abs(CSbbHiggs.P()+ CSyyHiggs.P()))
    print "Number of events pass: "+str(NumOfEventsPassed)
    print "Number of mismatched: "+str(wrong)
    #Creating new file
    outfile = TFile("plots_"+infileName+".root","RECREATE")
    print str(outfile)+" was created!"
    outfile.cd()
    #Write/Close
    test_CSyybb_p.Write()
    momentum_Higgs.Write()
    test_CSyybb_m.Write()
    deltaphi_colinsoper.Write()
    for hist in pthistos:
        pthistos[hist].Write()
    for hist in etahistos:
        etahistos[hist].Write()
    for hist in phihistos:
        phihistos[hist].Write()
    for hist in dpthistos:
        dpthistos[hist].Write()
    for hist in detahistos:
        detahistos[hist].Write()
    for hist in dphihistos:
        dphihistos[hist].Write()
    for hist in dRhistos:
        dRhistos[hist].Write()
    for hist in mhistos:
        mhistos[hist].Write()
    dEtayyj_yyjjH.Write()
    dEtayyj_yyjSH.Write()
    dEtayyjj_yyjSH.Write()
    dRyjHPTPH.Write()
    dEtayjHPTPH.Write()
    dphiyjHPTPH.Write()
    dRyjHPTJH.Write()
    dEtayjHPTJH.Write()
    dphiyjHPTJH.Write()
    dRyjSPTBH.Write()
    dEtayjSPTBH.Write()
    dphiyjSPTBH.Write()
    dRyjjHPTPH.Write()
    dEtayjjHPTPH.Write()
    dphiyjjHPTPH.Write()
    dRyjjSPTPH.Write()
    dEtayjjSPTPH.Write()
    dphiyjjSPTPH.Write()
    dRyyjHPTJH.Write()
    dEtayyjHPTJH.Write()
    dphiyyjHPTJH.Write()
    dRyyjSPTJH.Write()
    dEtayyjSPTJH.Write()
    dphiyyjSPTJH.Write()
    dRyjHPTBH.Write()
    dEtayjHPTBH.Write()
    dphiyjHPTBH.Write()
    dEtayyjjH.Write()
    dphiyyjjH.Write()
    dRyyjjH.Write()
    dEtajjH.Write()
    dRjjH.Write()
    dphijjH.Write()
    dEtayyH.Write()
    dRyyH.Write()
    dphiyyH.Write()
    DHCMH.Write()
    mjj125H.Write()
    DHH.Write()
    mjjHist.Write()
    mgammagammaHist.Write()
    JSSLH.Write()
    JSSSLH.Write()
    HPTH.Write()
    SLPTH.Write()
    JetLH.Write()
    JetSLH.Write()
    Jet_nHist.Write()
    Jet_ptHist.Write()
    Jet_mHist.Write()
    Jet_phiHist.Write()
    Jet_etaHist.Write()
    Gam_mHist.Write()
    Gam_phiHist.Write()
    Gam_etaHist.Write()
    Gam_ptHist.Write()
    Gam_nHist.Write()
    print str(outfile)+" was created!"
    outfile.Close()
    
    return outfile;

def analyzelist(flist,yybb_bcat_arg,listofobjectname,differentiallist):
    for Sample in flist:
        pthistos = {}
        etahistos = {}
        phihistos = {}
        dpthistos = {}
        detahistos = {}
        dphihistos = {}
        dRhistos = {}
        mhistos = {}
        cans = {}
        analyzedsample = analyze(Sample,yybb_bcat_arg,listofobjectname,differentiallist,pthistos,etahistos,phihistos,dpthistos,dphihistos,detahistos,dRhistos,mhistos,cans)

def plotall(plotlist,Histolist):
    files = {}
    stacks = {}
    legs = {}
    cans = {}
    norm = 1
    for histo in Histolist:
        gROOT.SetBatch(True)
        stacks[histo] = THStack("Stack_%s" %(histo),histo)
        legs[histo] = TLegend(0.2,0.3,0.9,0.9)  
        legs[histo].SetFillColor(0)  
    for Analyzed in plotlist:
        files[Analyzed] = TFile.Open("plots_"+Analyzed+".root", "READ")
        plots = {}
        for hist in Histolist:
            plots[hist] = files[Analyzed].Get(hist)
            plots[hist].Scale(norm/plots[hist].Integral())
            plots[hist].SetStats(0)
            plots[hist].SetLineColor(linecolors[Analyzed])
            legs[hist].AddEntry(plots[hist], leglabels[Analyzed], "L")
            stacks[hist].Add(plots[hist])
    for histos in Histolist:
        cans[histos] = TCanvas("Canvas_%s" %(histos),"Canvas_%s" %(histos),800,600)
        cans[histos].Divide(2,1,.01,0.1)
        cans[histos].cd(1)
        stacks[histos].Draw("nostack")
        stacks[histos].GetHistogram().GetXaxis().SetTitle("Mass [GeV]")
#        xaxistitle = stacks[histos].GetTitle().split("_")
#        if len(xaxistitle) == 2:
#            stacks[histos].GetHistogram().GetXaxis().SetTitle("%s %s" %(xaxistitle[0],xaxistitle[1]))
#        else:
#            stacks[histos].GetHistogram().GetXaxis().SetTitle("%s,%s  %s" %(xaxistitle[0],xaxistitle[1],xaxistitle[2]))
        stacks[histos].GetHistogram().GetYaxis().SetTitle("Events")
        stacks[histos].SetTitle("")
        cans[histos].cd(2)
        legs[histos].Draw("same")
#        cans[histos].SaveAs("/afs/cern.ch/user/c/cmosakow/Phase1/PythonScripts/plots_"+Analyzed+"/vars/"+hist+Analyzed+".pdf","pdf")
        cans[histos].SaveAs(histos+Analyzed+".pdf","pdf")
    for f in files.values():
        f.Close()


def plotallCS(plotlist,Histolist):
    files = {}
    stacks = {}
    ntuples = {}
    cssamps = {}
    legs = {}
    cans = {}
    cutflows = {}
    eventweight = {}
    mxaod = {}
    dxaod = {}
    xaod = {}
    scalingweight = {}
    eventweight = {}
    hgamweight = {}
    normweight = {}
    minitree = {}
    weight = {}
    massweights = {}
    crosssectbrratio = {}
    lumicross = {}
    entries = {}
    for histo in Histolist:
        gROOT.SetBatch(True)
        stacks[histo] = THStack("Stack_%s" %(histo),histo)
        legs[histo] = TLegend(0.2,0.3,0.9,0.9)  
        legs[histo].SetFillColor(0)  
    for Analyzed in plotlist:
        print "File:",Analyzed
        files[Analyzed] = TFile.Open("plots_"+Analyzed+".root", "READ")
        plots = {}
        for hist in Histolist:
            plots[hist] = files[Analyzed].Get(hist)
            print plots[hist].Integral()
            plots[hist].SetStats(0)
            plots[hist].SetLineColor(linecolors[Analyzed])
            legs[hist].AddEntry(plots[hist], leglabels[Analyzed], "L")
            stacks[hist].Add(plots[hist])
    for histos in Histolist:
        cans[histos] = TCanvas("Canvas_%s" %(histos),"Canvas_%s" %(histos),800,600)
        cans[histos].Divide(2,1,.01,0.1)
        cans[histos].cd(1)
        stacks[histos].Draw("nostack, hist")
        stacks[histos].GetHistogram().GetXaxis().SetTitle("Mass [GeV]")
        stacks[histos].GetHistogram().GetYaxis().SetTitle("Events")
        stacks[histos].SetTitle("")
        cans[histos].cd(2)
        legs[histos].Draw("same")
#        cans[histos].SaveAs("/afs/cern.ch/user/c/cmosakow/Phase1/PythonScripts/plots_"+Analyzed+"/vars/"+hist+Analyzed+".pdf","pdf")
        cans[histos].SaveAs(histos+Analyzed+".pdf","pdf")
    for f in files.values():
        f.Close()

def Pretty2D(fileslist,Histolist):
    gROOT.SetBatch(True)
    dfiles = {}
    cans = {}
    for files in fileslist:
        dfiles[files] = TFile.Open("plots_"+files+".root","READ")
        plots = {}
        for hist in Histolist:
            plots[hist] = dfiles[files].Get(hist)
            plots[hist].SetStats(0)
    for histos in Histolist:
        cans[histos] = TCanvas("Canvas_%s" %(histos), "Canvas_%s" %(histos),800,600)
        plots[histos].Draw("COLZ")
        cans[histos].SaveAs(histos+files+".pdf","pdf")
    for f in dfiles.values():
        f.Close()

def figofmerit(arg1,arg2,arg3):
    gROOT.SetBatch(True)
    bgplots = {}
    splots = {}
    fom = {}
    fom2 = {}
    fommax = {}
    fom2max = {}
    stacks = {}
    cans = {}
    legs = {}
    axes = {}
    for histos in arg3:
        stacks[histos] = THStack("stack %s" %(histos),"")
        legs[histos] = TLegend(.2,.3,.9,.9)
        legs[histos].SetFillColor(0)
    sampfile = TFile.Open("plots_"+arg1+".root","READ")
    bgfile = TFile.Open("plots_"+arg2+".root","READ")
    signtuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+arg1+".root")
    tree = signtuple.Get("cutflow_weighted")
    sigEvents = tree.GetBinContent(1)
    bkgntuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+arg2+".root")
    tree2 = bkgntuple.Get("cutflow_weighted")
    bkgEvents = tree2.GetBinContent(1)
    for hist in arg3:
        splots[hist] = sampfile.Get(hist)
#        print splots[hist]
        bgplots[hist] = bgfile.Get(hist)
#        print bgplots[hist]
        bgscale = (((3.9722*math.pow(10,1))*4.9729*math.pow(10,-1))/bkgEvents)*30*math.pow(10,3) 
        sigscale = (5.0/sigEvents)*30*math.pow(10,3)
        splots[hist].Scale(sigscale)
        bgplots[hist].Scale(bgscale)
        splots[hist].SetStats(0)
        bgplots[hist].SetStats(0)
        splots[hist].SetLineColor(linecolors[arg1])
        bgplots[hist].SetLineColor(linecolors[arg2])
        legs[hist].AddEntry(splots[hist], leglabels[arg1], "L")
        legs[hist].AddEntry(bgplots[hist], leglabels[arg2], "L")
        stacks[hist].Add(splots[hist])
        stacks[hist].Add(bgplots[hist])
        s = splots[hist].GetXaxis().GetXmax()
        fom[hist] = TH1F(hist,hist,30,0,s)
        fom2[hist] = TH1F(hist,hist,30,0,s)
        fom2[hist].SetLineColor(kBlack)
        fom[hist].SetTitle("")
        fom2[hist].SetTitle("")
        fom[hist].SetStats(0)
        fom2[hist].SetStats(0)
        for i in xrange(1,300):
            s = 0
            b = 0
            if bgplots[hist].Integral(0,i) != 0 and splots[hist].Integral(0,i) != 0:
                s=splots[hist].Integral(0,i) 
                b = bgplots[hist].Integral(0,i)
                eff= s/math.sqrt(s+b)
                eff2 = s/math.sqrt(b)
                fom[hist].SetBinContent(i,eff)
                fom2[hist].SetBinContent(i,eff2)
                s = 0
                b = 0
            else:
                continue
        legs[hist].AddEntry(fom[hist], "s/(s+b)^1/2","L")
        legs[hist].AddEntry(fom2[hist], "s/(b)^1/2","L")

    for hist in arg3:
        gROOT.SetBatch(True)
        cans[hist] = TCanvas("Canvas %s" %(hist),"Canvas %s" %(hist),1000,800)
        cans[hist].Divide(2,2,.01,.01)
        cans[hist].cd(1)
        stacks[hist].Draw("nostack")
        maxfombin = fom[hist].GetMaximumBin()
        center_of_fombin = stacks[hist].GetHistogram().GetXaxis().GetBinCenter(maxfombin)
        maxfom2bin = fom2[hist].GetMaximumBin()
        center_of_fom2bin = stacks[hist].GetXaxis().GetBinCenter(maxfom2bin)
        fommax[hist] = TLine(center_of_fombin,gPad.GetUymin(),center_of_fombin,gPad.GetUymax())
        fom2max[hist] = TLine(center_of_fom2bin,gPad.GetUymin(),center_of_fom2bin,gPad.GetUymax())
        fommax[hist].SetLineColor(kBlue)
        fom2max[hist].SetLineColor(kBlack)
        fommax[hist].Draw("same")
        fom2max[hist].Draw("same")
        stacks[hist].GetHistogram().GetXaxis().SetTitle(labels[hist])
        stacks[hist].GetHistogram().GetYaxis().SetTitle("Events")
        cans[hist].cd(2)
        legs[hist].Draw("same")
        cans[hist].cd(3)
        fom2[hist].SetMinimum(0)
        fom[hist].SetMinimum(0)
        fom2[hist].Draw("same")
        fom[hist].Draw("same")
        fom2[hist].GetXaxis().SetTitle(labels[hist])
        fom2[hist].GetYaxis().SetTitle("Figure of Merit")
        cans[hist].SaveAs("/afs/cern.ch/user/c/cmosakow/Phase1/PythonScripts/plots_"+arg1+"/"+hist+arg1+".pdf","pdf")

    sampfile.Close()
    bgfile.Close()

        
def listofcombos(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12):
    results = list(itertools.product(list1,list2))
    print "\n"
    print "Results:",results
    r2 = list(itertools.product(list3,list8))
    r3 = list(itertools.product(list3,list9))
    r4 = list(itertools.product(list3,list10))
    r5 = list(itertools.product(list4,list7))
    r6 = list(itertools.product(list4,list9))
    r7 = list(itertools.product(list4,list10))
    r8 = list(itertools.product(list5,list7))
    r9 = list(itertools.product(list5,list8))
    r10 = list(itertools.product(list5,list10))
    r11= list(itertools.product(list6,list7))
    r12= list(itertools.product(list7,list8))
    r13= list(itertools.product(list7,list9))
    sumr = results +r2 +r3+r4+r5+r6+r7+r8+r9+r10+r11+r12+r13
    r14 = list(itertools.product(sumr,list11))
    r15 = list(itertools.product(sumr,list12))
    totalresults = sumr + r14+r15
    print "\n"
    print "RESULTS:",totalresults
    return totalresults;

def listgen(listofobjects,otherlist):
    masterlist = []    
    for vector in listofobjects:
        masterlist.append(vector+"_pt")
        masterlist.append(vector+ "_eta")
        masterlist.append(vector+ "_phi")
        masterlist.append(vector+ "_m")
    for vector in otherlist:
        masterlist.append(vector+"_d(pt)")
        masterlist.append(vector+"_dR")
        masterlist.append(vector+"_d(eta)")
        masterlist.append(vector+ "_d(phi)")
    return masterlist;

def KStest(sig300reweight,sig350reweight,sig275reweight,bkg,histlist,bkgscale,sig350,sig275,sig300):
    gROOT.SetBatch(True)
    bgplots = {}
    splots = {}
    stacks = {}
    cans = {}
    legs = {}
    splots350 = {}
    splots275 = {}
    c1 = TCanvas("c1","",800,600)
    c2 = TCanvas("c2","",800,600)
    c3 = TCanvas("c3","",800,600)
    h1 = TH1F("ADT","ADT",100,0,7000)
    h2 = TH1F("ADT350","ADT350",100,0,7000)
    for histos in histlist:
        legs[histos] = TLegend(.2,.3,.9,.9)
        legs[histos].SetFillColor(0)
    sampfile = TFile.Open("plots_"+sig300reweight+".root","READ")
    bgfile = TFile.Open("plots_"+bkg+".root","READ")
    sampfile350 = TFile.Open("plots_"+sig350reweight+".root")
    sampfile275 = TFile.Open("plots_"+sig275reweight+".root")
    signtuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+sig300+".root")
    ntuple350 = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+sig350+".root")
    ntuple275 = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+sig275+".root")
    tree = signtuple.Get("cutflow_weighted")
    xAOD = tree.GetBinContent(1)
    DxAOD = tree.GetBinContent(2)
    MxAOD = tree.GetBinContent(3)
    scalefact = (DxAOD/xAOD)*(1/MxAOD)
    tree350 = ntuple350.Get("cutflow_weighted")
    xAOD350 = tree350.GetBinContent(1)
    DxAOD350= tree350.GetBinContent(2)
    MxAOD350 = tree350.GetBinContent(3)
    scalefact350 = (DxAOD350/xAOD350)*(1/MxAOD350)
    tree275 = ntuple275.Get("cutflow_weighted")
    xAOD275 = tree275.GetBinContent(1)
    DxAOD275 = tree275.GetBinContent(2)
    MxAOD275 = tree275.GetBinContent(3)
    scalefact275 = (DxAOD275/xAOD275)*(1/MxAOD275)
    bkgntuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+bkgscale+".root")
    tree2 = bkgntuple.Get("cutflow_weighted")
    bgxAOD = tree2.GetBinContent(1)
    bgDxAOD = tree2.GetBinContent(2)
    bgMxAOD = tree2.GetBinContent(3)
    bgscalefact = (bgDxAOD/bgxAOD)*(1/bgMxAOD)
    maxminlistofhistos = {}
    minmaxsamples = 0
    goodlistofhistos = []
    nsamples = 0
    finallist = []
    for hist in histlist:
        splots[hist] = sampfile.Get(hist)
        bgplots[hist] = bgfile.Get(hist)
        splots350[hist] = sampfile350.Get(hist)
        splots275[hist] = sampfile275.Get(hist)
        splots350[hist].SetTitle("")
        splots275[hist].SetTitle("")
        splots[hist].SetTitle("")
        bgplots[hist].SetTitle("")
        splots350[hist].Scale((5.0*scalefact350)*30*math.pow(10,3))
        splots275[hist].Scale((5.0*scalefact275)*30*math.pow(10,3))
        splots[hist].SetStats(0)
        splots[hist].Scale((5.0*scalefact)*30*math.pow(10,3))
        bgplots[hist].Scale(((3.9722*math.pow(10,1))*4.9729*math.pow(10,-1)*bgscalefact)*(30*math.pow(10,3)))
        bgplots[hist].SetStats(0)
        splots[hist].SetLineColor(linecolors[sig300reweight])
        bgplots[hist].SetLineColor(linecolors[bkgscale])
        legs[hist].AddEntry(splots[hist], leglabels[sig300reweight], "L")
        legs[hist].AddEntry(bgplots[hist], leglabels[bkgscale], "L")
        adt = splots[hist].AndersonDarlingTest(bgplots[hist],"T")
        minmaxadt = splots350[hist].AndersonDarlingTest(splots275[hist],"T")
        if minmaxadt <150 and adt >70:
            maxminlistofhistos[hist] = [minmaxadt,adt]
            minmaxsamples +=1
        h1.Fill(adt)
        h2.Fill(minmaxadt)
    s = sorted(maxminlistofhistos.items(), key=lambda i: i[1][1])
    print s
    print len(s)
#    c1.cd()
#    h1.Draw()
#    h1.SetTitle("ADT Between 300 GeV Signal Sample and Data")
#    c1.SaveAs("ADTdata.pdf","pdf")
#    c2.cd()
#    h2.Draw()
#    h2.SetTitle("ADT Between 350 GeV and 275 GeV Signal Samples")
#    c2.SaveAs("ADT350.pdf","pdf")
#    for hist in histlist:
#        cans[hist] = TCanvas("Canvas %s" %(hist),"Canvas %s" %(hist),1000,800)
#        cans[hist].Divide(2,1,.01,.01)
#        cans[hist].cd(1)
#        splots[hist].Draw()
#        bgplots[hist].Draw("same")
#        cans[hist].cd(2)
#        legs[hist].Draw("same")
#        cans[hist].SaveAs(hist+".pdf","pdf")
    return maxminlistofhistos;


def figofmeritcuts(sig,sigMC,data,listofhistos):
    gROOT.SetBatch(True)
    bgplots = {}
    splots = {}
    fom = {}
    fom2 = {}
    fommax = {}
    fom2max = {}
    stacks = {}
    cans = {}
    legs = {}
    axes = {}
    for histos in listofhistos:
        stacks[histos] = THStack("stack %s" %(histos),"")
        legs[histos] = TLegend(.2,.3,.9,.9)
        legs[histos].SetFillColor(0)
    sampfile = TFile.Open("plots_"+sigMC+".root","READ")
    bgfile = TFile.Open("plots_"+data,"READ")
    signtuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+sig+".root")
    tree = signtuple.Get("cutflow_weighted")
    xAOD = tree.GetBinContent(1)
    DxAOD = tree.GetBinContent(2)
    MxAOD = tree.GetBinContent(3)
    scalefact = (DxAOD/xAOD)*(1/MxAOD)  
    tree2 = bgfile.Get("jj_pt")
    bkgEvents = tree2.GetEntries()
    for hist in listofhistos:
        splots[hist] = sampfile.Get(hist)
        bgplots[hist] = bgfile.Get(hist)
        bgscale = 28343
        sigscale = (5.0*scalefact)*12.171*math.pow(10,3)
#        splots[hist].Scale(sigscale)
#        bgplots[hist].Scale(bgscale/bkgEvents)
        splots[hist].SetStats(0)
        bgplots[hist].SetStats(0)
        splots[hist].SetLineColor(linecolors[sigMC])
        bgplots[hist].SetLineColor(linecolors[data])
        legs[hist].AddEntry(splots[hist], leglabels[sigMC], "L")
        legs[hist].AddEntry(bgplots[hist], leglabels[data], "L")
        stacks[hist].Add(splots[hist])
        stacks[hist].Add(bgplots[hist])
        s = splots[hist].GetXaxis().GetXmax()
        fom[hist] = TH1F(hist,hist,30,0,s)
        fom2[hist] = TH1F(hist,hist,30,0,s)
        fom[hist].SetLineColor(kBlue)
        fom2[hist].SetLineColor(kBlack)
        fom[hist].SetTitle("")
        fom2[hist].SetTitle("")
        fom[hist].SetStats(0)
        fom2[hist].SetStats(0)
        for i in xrange(1,32):
            s = 0
            b = 0
            if bgplots[hist].Integral(0,i) != 0 and splots[hist].Integral(0,i) != 0:
                s=splots[hist].Integral(0,i) 
                b = bgplots[hist].Integral(0,i)
                eff= s/math.sqrt(s+b)
                eff2 = s/math.sqrt(b)
                fom[hist].SetBinContent(i,eff)
                fom2[hist].SetBinContent(i,eff2)
                
            else:
                continue
            if i == 31:
                print eff
                print eff2
                print "s:",s
                print "b:",b
        fom2[hist].SetLineColor(kBlack)
        fom[hist].SetTitle("")
        fom2[hist].SetTitle("")
        fom[hist].SetStats(0)
        fom2[hist].SetStats(0)
        legs[hist].AddEntry(fom[hist], "s/(s+b)^1/2","L")
        legs[hist].AddEntry(fom2[hist], "s/(b)^1/2","L")

    for hist in listofhistos:
        gROOT.SetBatch(True)
        cans[hist] = TCanvas("Canvas %s" %(hist),"Canvas %s" %(hist),1000,800)
        cans[hist].Divide(2,2,.01,.01)
        cans[hist].cd(1)
        stacks[hist].Draw("nostack")
        maxfombin = fom[hist].GetMaximumBin()

        center_of_fombin = stacks[hist].GetHistogram().GetXaxis().GetBinCenter(maxfombin)
        maxfom2bin = fom2[hist].GetMaximumBin()
        center_of_fom2bin = stacks[hist].GetXaxis().GetBinCenter(maxfom2bin)
        fommax[hist] = TLine(center_of_fombin,gPad.GetUymin(),center_of_fombin,gPad.GetUymax())
        fom2max[hist] = TLine(center_of_fom2bin,gPad.GetUymin(),center_of_fom2bin,gPad.GetUymax())
        fommax[hist].SetLineColor(kBlue)
        fom2max[hist].SetLineColor(kBlack)
        fommax[hist].Draw("same")
        fom2max[hist].Draw("same")
        stacks[hist].GetHistogram().GetXaxis().SetTitle(labels[hist])
        stacks[hist].GetHistogram().GetYaxis().SetTitle("Events")
        cans[hist].cd(2)
        legs[hist].Draw("same")
        cans[hist].cd(3)
        fom2[hist].SetMinimum(0)
        fom[hist].SetMinimum(0)
        fom2[hist].Draw("same")
        fom[hist].Draw("same")
        fom2[hist].GetXaxis().SetTitle(labels[hist])
        fom2[hist].GetYaxis().SetTitle("Figure of Merit")
        cans[hist].SaveAs("/afs/cern.ch/user/c/cmosakow/www/"+hist+data+sigMC+".pdf","pdf")
        print "Hist:", hist
        print "\n"
        print "s/root(s+b)",center_of_fombin
        print "s/root(b)",center_of_fom2bin
        maxfombin = 0
        maxfom2bin = 0
        center_of_fombin = 0
        center_of_fom2bin = 0
    sampfile.Close()
    bgfile.Close()


def compareMC_data(bgMCreweight, bgDatareweight, histolist,bgMC):
    gROOT.SetBatch(True)
    cans = {}
    bgMCplots = {}
    bgDataplots = {}
    legs = {}
    bgMCfile = TFile.Open("plots_"+bgMCreweight+".root","READ")
    bgDatafile = TFile.Open("plots_"+bgDatareweight+".root","READ")
    bkgntuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+bgMC+".root")
    tree2 = bkgntuple.Get("cutflow_weighted")
    bgxAOD = tree2.GetBinContent(1)
    bgDxAOD = tree2.GetBinContent(2)
    bgMxAOD = tree2.GetBinContent(3)
    bgscalefact = (bgDxAOD/bgxAOD)*(1/bgMxAOD)
    for hist in histolist:
        legs[hist] = TLegend(.2,.3,.9,.9)
        legs[hist].SetFillColor(0)
    for hist in histolist: 
        bgDataplots[hist] = bgDatafile.Get(hist)
        bgMCplots[hist] = bgMCfile.Get(hist)
#        bgscale = (((3.9722*math.pow(10,1))*4.9729*math.pow(10,-1))*bgscalefact)*12.171*math.pow(10,3) 
#        bgMCplots[hist].Scale(bgscale)
        bgMCplots[hist].Scale(1/bgMCplots[hist].Integral())
        bgDataplots[hist].Scale(1/bgDataplots[hist].Integral())
        bgDataplots[hist].SetStats(0)
        bgMCplots[hist].SetStats(0)
        bgDataplots[hist].SetLineColor(linecolors[bgDatareweight])
        bgMCplots[hist].SetLineColor(linecolors[bgMCreweight])
        legs[hist].AddEntry(bgDataplots[hist], leglabels[bgDatareweight], "L")
        legs[hist].AddEntry(bgMCplots[hist], leglabels[bgMCreweight], "L")
        bgDataplots[hist].SetStats(0)
        bgMCplots[hist].SetStats(0)
    for histos in histolist:
        cans[histos] = TCanvas("Canvas_%s" %(histos),histos,800,600)
        cans[histos].Divide(2,1,.01,0.1)
        cans[histos].cd(1)
        bgMCplots[histos].Draw()
        bgDataplots[histos].Draw("same")
#"E1"
#"LSAME"
#        xaxistitle = cans[histos].GetTitle().split("_")
#        if len(xaxistitle) == 2:
#            bgDataplots[histos].GetXaxis().SetTitle("%s %s" %(xaxistitle[0],xaxistitle[1]))
#        else:
#            bgDataplots[histos].GetXaxis().SetTitle("%s,%s  %s" %(xaxistitle[0],xaxistitle[1],xaxistitle[2]))
        bgDataplots[histos].GetYaxis().SetTitle("Events")
        bgDataplots[histos].GetXaxis().SetTitle("pT [GeV]")
        cans[histos].cd(2)
        legs[histos].Draw("same")
        cans[histos].SetTitle("")
        cans[histos].SaveAs(histos+bgMC+bgDatareweight+".pdf","pdf")

def MG_plot(mg1,mg2,mg3,mg4,mg5,data,histolist):
    gROOT.SetBatch(True)
    cans = {}
    Dataplots = {}
    legs = {}
    stacks = {}
    MG1 = {}
    MG2 = {}
    MG3 = {}
    MG4 = {}
    MG5 = {}
    for hist in histolist:
        legs[hist] = TLegend(.2,.3,.9,.9)
        legs[hist].SetFillColor(0)
        stacks[hist] = THStack(hist,hist)
    MG1file = TFile.Open("plots_"+mg1+".root","READ")
    MG2file = TFile.Open("plots_"+mg2+".root","READ")
    MG3file = TFile.Open("plots_"+mg3+".root","READ")
    MG4file = TFile.Open("plots_"+mg4+".root","READ")
    MG5file = TFile.Open("plots_"+mg5+".root","READ") 
    Datafile = TFile.Open("plots_"+data+".root","READ")
    MG1ntuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+mg1+".root")
    MG2ntuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+mg2+".root")
    MG3ntuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+mg3+".root")
    MG4ntuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+mg4+".root")
    MG5ntuple = TFile.Open("/afs/cern.ch/user/e/ebrost/public/tiny/"+mg5+".root")
    tree = MG1ntuple.Get("cutflow_weighted")
    MG1xAOD = tree.GetBinContent(1)
    MG1DxAOD = tree.GetBinContent(2)
    MG1MxAOD = tree.GetBinContent(3)
    MG1scale = (MG1DxAOD/MG1xAOD)*(1/MG1MxAOD)
    tree2 = MG1ntuple.Get("cutflow_weighted")
    MG2xAOD = tree2.GetBinContent(1)
    MG2DxAOD = tree2.GetBinContent(2)
    MG2MxAOD = tree2.GetBinContent(3)
    MG2scale = (MG2DxAOD/MG2xAOD)*(1/MG2MxAOD)
    tree3 = MG1ntuple.Get("cutflow_weighted")
    MG3xAOD = tree3.GetBinContent(1)
    MG3DxAOD = tree3.GetBinContent(2)
    MG3MxAOD = tree3.GetBinContent(3)
    MG3scale = (MG3DxAOD/MG3xAOD)*(1/MG3MxAOD)
    tree4 = MG1ntuple.Get("cutflow_weighted")
    MG4xAOD = tree4.GetBinContent(1)
    MG4DxAOD = tree4.GetBinContent(2)
    MG4MxAOD = tree4.GetBinContent(3)
    MG4scale = (MG4DxAOD/MG4xAOD)*(1/MG4MxAOD)
    tree5 = MG1ntuple.Get("cutflow_weighted")
    MG5xAOD = tree5.GetBinContent(1)
    MG5DxAOD = tree5.GetBinContent(2)
    MG5MxAOD = tree5.GetBinContent(3)
    MG5scale = (MG1DxAOD/MG1xAOD)*(1/MG1MxAOD)
    for hist in histolist:
        MG1[hist] = MG1file.Get(hist)
        MG2[hist] = MG2file.Get(hist)
        MG3[hist] = MG3file.Get(hist)
        MG4[hist] = MG4file.Get(hist)
        MG5[hist] = MG5file.Get(hist)
        Dataplots[hist] = Datafile.Get(hist)
        Dataplots[hist].SetStats(0)
        Dataplots[hist].SetLineColor(linecolors[data])
        MG1[hist].SetStats(0)
        MG2[hist].SetStats(0)
        MG3[hist].SetStats(0)
        MG4[hist].SetStats(0)
        MG5[hist].SetStats(0)
        MG1[hist].GetXaxis().SetTitle("")
        MG2[hist].GetXaxis().SetTitle("")
        MG3[hist].GetXaxis().SetTitle("")
        MG4[hist].GetXaxis().SetTitle("")
        MG5[hist].GetXaxis().SetTitle("")
        MG1[hist].SetLineColor(linecolors[mg1])
        MG2[hist].SetLineColor(linecolors[mg2])
        MG3[hist].SetLineColor(linecolors[mg3])
        MG4[hist].SetLineColor(linecolors[mg4])
        MG5[hist].SetLineColor(linecolors[mg5])
        MG1[hist].Scale(1.1327*math.pow(10,2)*2*(MG1scale)*12.171*math.pow(10,3))
        MG2[hist].Scale(2.0504*math.pow(10,2)*2*(MG2scale)*12.171*math.pow(10,3))
        MG3[hist].Scale(5.7727*math.pow(10,3)*2*(MG3scale)*12.171*math.pow(10,3))
        MG4[hist].Scale(3.0429*math.pow(10,-2)*2*(MG4scale)*12.171*math.pow(10,3))
        MG5[hist].Scale(1.2404*math.pow(10,-1)*2*(MG5scale)*12.171*math.pow(10,3))
        legs[hist].AddEntry(Dataplots[hist], leglabels[data], "L")
        legs[hist].AddEntry(MG1[hist],leglabels[mg1],"L")
        legs[hist].AddEntry(MG2[hist],leglabels[mg2],"L")
        legs[hist].AddEntry(MG3[hist],leglabels[mg3],"L")
        legs[hist].AddEntry(MG4[hist],leglabels[mg4],"L")
        legs[hist].AddEntry(MG5[hist],leglabels[mg5],"L")
    for hist in histolist:
        cans[hist] = TCanvas(hist,hist,800,600)
        cans[hist].Divide(2,1,.01,.01)
        cans[hist].cd(1)
        Dataplots[hist].Draw("E1")
        xaxistitle = cans[hist].GetTitle().split("_")
        if len(xaxistitle) == 2:
            Dataplots[hist].GetXaxis().SetTitle("%s %s" %(xaxistitle[0],xaxistitle[1]))
        else:
            Dataplots[hist].GetXaxis().SetTitle("%s,%s  %s" %(xaxistitle[0],xaxistitle[1],xaxistitle[2]))
        cans[hist].cd(2)
        legs[hist].Draw("same")
        cans[hist].SaveAs(hist+"MadGraph_Data.pdf","pdf")
        
