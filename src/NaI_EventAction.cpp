#include "G4AnalysisManager.hh"
#include "G4EventManager.hh"
#include "G4RandomTools.hh"
#include "G4TrackingManager.hh"
#include "NaI_EventAction.h"

#include "G4SystemOfUnits.hh"
NaI_EventAction::NaI_EventAction() {}

NaI_EventAction::~NaI_EventAction() {}

//static int eventsWithDeposition = 0; //New-old
 int dE_dep=0;
 int E_dep = 0;

void NaI_EventAction::BeginOfEventAction(const G4Event *event) {
  // TODO : All the required logic you want to do at the start
  //        of each event
  eDep_E = 0;
  eDep_dE = 0;
 // std::cout << "===================== EVENT ID : " << event->GetEventID() << "=========================" << std::endl;
}

void NaI_EventAction::EndOfEventAction(const G4Event *event) {
  // TODO : All the required logic you want to do at the end
  //        of each event

 // std::cout << "Total Energy deposited in event : " << event->GetEventID() << " :in E = " << eDep_E << "& in dE = "<< eDep_dE << std::endl;

  /*G4AnalysisManager *analMan = G4AnalysisManager::Instance();
  if (eDep > 0) {
    //eventsWithDeposition++;               //New-old
    double smeardEDep = G4RandGauss::shoot(eDep, 0.04385);
    analMan->FillH1(0, eDep * 1000.);
    analMan->FillH1(1, smeardEDep * 1000.);

    analMan->FillNtupleDColumn(0, 0, eDep * 1000.);
    analMan->FillNtupleDColumn(0, 1, smeardEDep * 1000.);
    analMan->AddNtupleRow(0);
    
  }*/
  
  	/*G4AnalysisManager *analMan = G4AnalysisManager::Instance();
  	 {
  	   This fills even zero energy deposition events, so the plot is not neat. Need to zoom.
  	   if(eDep_dE>0) {
  	   dE_dep++;
  	   analMan->FillNtupleDColumn(1, eDep_dE);}
  	   if(eDep_E>0){ 
  	   E_dep++ ; 
  	   analMan->FillNtupleDColumn(0, eDep_E );}
  	   
  	   analMan->FillH2(0,eDep_dE, eDep_E);
  	   analMan->FillNtupleDColumn(1, eDep_dE);
  	   analMan->FillNtupleDColumn(0, eDep_E );
    	   analMan->AddNtupleRow();
    	   }*/
    	   G4AnalysisManager *analMan = G4AnalysisManager::Instance();{
		bool hasEnergyDeposited = false;

		if (eDep_dE > 0) {
		    dE_dep++;
		    analMan->FillNtupleDColumn(1, eDep_dE);
		    //double smeardEDep_dE = G4RandGauss::shoot(eDep_dE, 0.04385);
		    //analMan->FillNtupleDColumn(3,smeardEDep_dE * 1000.);
		    //analMan->FillH1(1, smeardEDep_dE * 1000.);
		    hasEnergyDeposited = true;}
		
		if (eDep_E > 0) {
		    E_dep++;
		    analMan->FillNtupleDColumn(0, eDep_E);
		    //double smeardEDep_E = G4RandGauss::shoot(eDep_E, 0.004385);
		   // analMan->FillNtupleDColumn(2,smeardEDep_E*1000);
		    //analMan->FillH1(0, smeardEDep_E * 1000.);
		    hasEnergyDeposited = true;}
		

		if (hasEnergyDeposited) {
		    analMan->FillH2(0, eDep_E , eDep_dE );
		    analMan->AddNtupleRow(0);}

    	 }
    G4cout << "End of Event " << event->GetEventID() << ":" << G4endl;
    G4cout << "Total energy deposit (E): " << eDep_E / MeV << " MeV" << G4endl;
    G4cout << "Total energy deposit (dE): " << eDep_dE / MeV << " MeV" << G4endl;
  	
  /*if (event->GetEventID() == 99999) {
    std::cout << "\n===== Summary: Events with energy deposited: " << eventsWithDeposition << " / 100000 =====\n";
    }*/

if (event->GetEventID() == 999999) {
    std::cout << "\n===== Summary: Events with energy deposited in dE: " << dE_dep << " / 1000000 =====\n";
    std::cout << "\n===== Summary: Events with energy deposited in E: " << E_dep << " / 1000000 =====\n";
    }

}

//void NaI_EventAction::AddEnergy(double val) { eDep += val; }
void NaI_EventAction:: AddEdep_E(G4double edep) { eDep_E += edep;} 
void NaI_EventAction::AddEdep_dE(G4double edep) { eDep_dE += edep;}

	
	
