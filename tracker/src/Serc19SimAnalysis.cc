
//////////////////////////////////////////////////////////////////////////////
//
//
//
//
////////////////////////////////////////////////////////////////////////////////

#include "Serc19SimAnalysis.hh"
#include "Serc19PrimaryGeneratorAction.hh"

Serc19SimAnalysis *Serc19SimAnalysis::AnPointer;

Serc19SimAnalysis::Serc19SimAnalysis() {

  AnPointer = this;
  //  text_inputFile = sprefix.data();  // file name without .root
  irun = 100;
  ievt		= 0;
  ievent		= 0;
  
  pRootFile	=0;
  inputRootFile=0;
  
  pEventTree	=0;
  inputEventTree=0;
  
  rtree=0;


}

void Serc19SimAnalysis::OpenRootfiles(G4String outfile) {

  ievent=0;
  //  G4String titleop= outfile;
  //  titleop.append(".root");
  //cout <<"title top=========================="<<titleop<<endl;
  pRootFile = new TFile(outfile, "RECREATE"); //VALGRIND
  
  if (!pRootFile) {
    G4cout << "Error opening root file !" << G4endl;
    exit(-1);
  } else {
    cout<< "Output stored in root file:"<< outfile <<endl;
  }
  
  pEventTree = new TTree("T1", "SERC2019_CMS");
  
  pEventTree->Branch("irun",&irun,"irun/i");
  pEventTree->Branch("ievt",&ievt,"ievt/i");
  
  pEventTree->Branch("ngent",&ngent,"ngent/i");
  pEventTree->Branch("pidin",pidin,"pidin[ngent]/I");
  pEventTree->Branch("ievt_wt",&ievt_wt,"ievt_wt/F"); //Generated events with different weight
  pEventTree->Branch("momin",momin,"momin[ngent]/F");
  pEventTree->Branch("thein",thein,"thein[ngent]/F");
  pEventTree->Branch("phiin",phiin,"phiin[ngent]/F");
  pEventTree->Branch("posxin",posxin,"posxin[ngent]/F");
  pEventTree->Branch("posyin",posyin,"posyin[ngent]/F");
  pEventTree->Branch("poszin",poszin,"poszin[ngent]/F");
  
  //SImulation output of tracker
  pEventTree->Branch("nsimhtTk", &nsimhtTk, "nsimhtTk/i");
  pEventTree->Branch("detidTk", detidTk, "detidTk[nsimhtTk]/i");
  pEventTree->Branch("simpdgidTk", simpdgidTk, "simpdgidTk[nsimhtTk]/I");
  pEventTree->Branch("simtimeTk", simtimeTk, "simtimeTk[nsimhtTk]/F");
  pEventTree->Branch("simenrTk", simenrTk, "simenrTk[nsimhtTk]/F");
  pEventTree->Branch("simvxTk", simvxTk, "simvxTk[nsimhtTk]/F");
  pEventTree->Branch("simvyTk", simvyTk, "simvyTk[nsimhtTk]/F");
  pEventTree->Branch("simvzTk", simvzTk, "simvzTk[nsimhtTk]/F");
  pEventTree->Branch("simpxTk", simpxTk, "simpxTk[nsimhtTk]/F");
  pEventTree->Branch("simpyTk", simpyTk, "simpyTk[nsimhtTk]/F");
  pEventTree->Branch("simpzTk", simpzTk, "simpzTk[nsimhtTk]/F");

  //SImulation output of ECAL
  pEventTree->Branch("nsimhtEC", &nsimhtEC, "nsimhtEC/i");
  pEventTree->Branch("detidEC", detidEC, "detidEC[nsimhtEC]/i");
  pEventTree->Branch("energyEC", energyEC, "energyEC[nsimhtEC]/F");
  pEventTree->Branch("thetaEC", thetaEC, "thetaEC[nsimhtEC]/F");
  pEventTree->Branch("phiEC", phiEC, "phiEC[nsimhtEC]/F");

  //SImulation output of HCAL
  pEventTree->Branch("nsimhtHL", &nsimhtHL, "nsimhtHL/i");
  pEventTree->Branch("detidHL", detidHL, "detidHL[nsimhtHL]/l");
  pEventTree->Branch("timeHL", timeHL, "timeHL[nsimhtHL]/i");

  double pibytwo = acos(-1.)/2.;
  char name[100];
  char title[100];
  pPosRT = new TH1F("pPosRT", "pPosRT", 120, 0.0, 900.0);
  pPosET = new TH1F("pPosET", "pPosET", 120, -pibytwo, pibytwo);
  pPosPT = new TH1F("pPosPT", "pPosPT", 120, -2.0, 2.0);
  
  pPosRE = new TH1F("pPosRE", "pPosRE", 120, 900.0, 1600.0);
  pPosEE = new TH1F("pPosEE", "pPosEE", 120, -pibytwo, pibytwo);
  pPosPE = new TH1F("pPosPE", "pPosPE", 120, -2.0, 2.0);
  
  pPosRH = new TH1F("pPosRH", "pPosRH", 120, 1800.0, 3000.0);
  pPosEH = new TH1F("pPosEH", "pPosEH", 120, -pibytwo, pibytwo);
  pPosPH = new TH1F("pPosPH", "pPosPH", 120, -2.0, 2.0);

  for (int ij=0; ij<(int)nsilayer; ij++) {
    sprintf(name, "etaphi_trkenergy_L%i", ij);
    sprintf(title, "Trk Energy in 2D for L%i", ij);
    h_etaphi_trkenergy[ij] = new TH2F(name, title, 120, -1.2, 1.2, 120, -pibytwo, pibytwo);
    
    sprintf(name, "trkenergy_L%i", ij);
    sprintf(title, "Trk Energy in L%i (log10(KeV))", ij);
    h_trkenergy[ij] = new TH1F(name, title, 120, -1.0, 4.0); //log scale

    sprintf(name, "trkabsenergy_L%i", ij);
    sprintf(title, "Trk Abs Energy in L%i (log10(KeV))", ij);
    h_trkabsenergy[ij] = new TH1F(name, title, 120, -1.0, 4.0); //log scale
  }

  h_parfenergy = new TH1F("parfenergy", "Energy in modulator (log10(MeV))", 120, -2.0, 3.0);
  h_ecalenergy = new TH1F("ecalenergy", "Energy in ECAL (log10(GeV))", 120, -2.0, 3.0);
  h_2decalenergy = new TH2F("h2decalenergy", "Energy in ECAL (GeV)", 120, -1.2, 1.2, 120, -pibytwo, pibytwo);
  h_g10energy = new TH1F("g10energy", "Energy in G10 (log10(MeV))", 120, -3.0, 3.0);

  for (int ij=0; ij<(int)nhcalLayer; ij++) {
    sprintf(name, "hcalenergy_L%i", ij);
    sprintf(title, "Hcal Energy in L%i (log10(keV))", ij);
    h_hcalenergy[ij] = new TH1F(name, title, 120, 0.0, 6.0);

    sprintf(name, "hcalabsenergy_L%i", ij);
    sprintf(title, "Hcal Abs Energy in L%i (log10(MeV))", ij);
    h_hcalabsenergy[ij] = new TH1F(name, title, 120, 0.0, 6.0);

    sprintf(name, "hcal2denergy_L%i", ij);
    sprintf(title, "Hcal Energy in L%i (keV)", ij);
    h2d_hcalenergy[ij] = new TH2F(name, title, 120, -2.0, 2.0, 120, -pibytwo, pibytwo);

    sprintf(name, "hcalabs2denergy_L%i", ij);
    sprintf(title, "Hcal Abs Energy in L%i (MeV)", ij);
    h2d_hcalabsenergy[ij] = new TH2F(name, title, 120, -2.0, 2.0, 120, -pibytwo, pibytwo);
  }
}

void Serc19SimAnalysis::CloseRootfiles() {

  if (pRootFile) {
    pRootFile->cd();

    pRootFile->Write(); //VALGRIND

    pRootFile->Close();
    delete pRootFile; pRootFile=0;
    
    cout << "Root file is writen and closed properly" << endl;

  } else {
    cout << "Root file not made !" << endl;
  }
}

Serc19SimAnalysis::~Serc19SimAnalysis() {
  cout <<"end of Serc19simanalysis detructor "<<endl;
}
