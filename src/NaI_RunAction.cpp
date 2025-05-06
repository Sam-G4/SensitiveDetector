#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "NaI_RunAction.h"
NaI_RunAction::NaI_RunAction() {}

NaI_RunAction::~NaI_RunAction() {}

void NaI_RunAction::BeginOfRunAction(const G4Run* aRun) {
	//G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  G4AnalysisManager *analMan = G4AnalysisManager::Instance();
  analMan->OpenFile("output.root");
  G4cout << "Opened output.root for writing analysis data." << G4endl;
  /*analMan->CreateH1("edep", "Energy Deposited", 50, 0, 1000);
  analMan->CreateH1("edepSmeared", "Energy Deposited with Detector resolution", 50, 0, 1000);
  analMan->CreateNtuple("ftree", "A simple tree for NaI Data");
  analMan->CreateNtupleDColumn("eDep");
  analMan->CreateNtupleDColumn("eDepSmeared");
  analMan->FinishNtuple();*/
  
  //New
  	analMan->CreateH2("dE-E", "dE v/s E", 1000, 0, 20, 1000, 0, 20);
  	//analMan->CreateH1("edepSmeared_E", "Energy Deposited with Detector resolution", 50, 0, 1000);
  	//analMan->CreateH1("edepSmeared_dE", "Energy Deposited with Detector resolution", 50, 0, 1000);
  	analMan->CreateNtuple("ftree", "A simple tree for Si-CD Data");
  analMan->CreateNtupleDColumn("eDep_E");
  analMan->CreateNtupleDColumn("eDep_dE");
  //analMan->CreateNtupleDColumn("eDepSmeared_E");
  //analMan->CreateNtupleDColumn("eDepSmeared_dE");
  analMan->FinishNtuple();
}

void NaI_RunAction::EndOfRunAction(const G4Run *) {

  G4AnalysisManager *analMan = G4AnalysisManager::Instance();
  analMan->Write();
  analMan->CloseFile();
  G4cout << "Analysis data written and file closed." << G4endl;
}
