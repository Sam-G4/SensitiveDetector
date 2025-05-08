#include "NaI_SensitiveDetector.h"
#include "NaI_SensitiveDetector_Hit.h"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
//#include "G4Alpha.hh"

//#include "G4UserEventAction.hh"
//#include "NaI_EventAction.h"
//#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"


// Initialize static collection ID
G4int NaI_SensitiveDetector::fCollectionID = -1;

NaI_SensitiveDetector::NaI_SensitiveDetector(const G4String& name) : G4VSensitiveDetector(name), fHitsCollection(nullptr)
{
	collectionName.insert("HitsCollection");	
}

NaI_SensitiveDetector::~NaI_SensitiveDetector() {}

//called at begining of each event
void NaI_SensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new G4THitsCollection<NaI_SensitiveDetector_Hit>(
        SensitiveDetectorName, collectionName[0]);

    if (fCollectionID < 0) {
        fCollectionID = GetCollectionID(0);
    }

    hce->AddHitsCollection(fCollectionID, fHitsCollection);
}

//called every time a step is made
G4bool NaI_SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep == 0) return false;

    auto hit = new NaI_SensitiveDetector_Hit();
    hit->SetEnergy(edep);

    G4int detID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
    hit->SetDetectorID(detID);

    G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
    hit->SetPosition(pos);
    
    //to get particle name
    G4Track* track = step->GetTrack();
    G4String particleName = track->GetDefinition()->GetParticleName();
    hit->SetParticleName(particleName);
  
    //hit->Print( ); //prints hits in one event(here usually 4) for all events
  
    fHitsCollection->insert(hit);
    return true;
}


//called at end of event
void NaI_SensitiveDetector::EndOfEvent(G4HCofThisEvent*) {

std::cout << " =======================================" << std::endl;
 G4int nHits= fHitsCollection->entries();
 std::cout <<"No. of Hits : "<< nHits <<std::endl;
 
 //to fill histogram without using EventAction
 	double dE = 0, E = 0;
 	for (unsigned int i = 0; i < nHits; i++) {
 	(*fHitsCollection)[i]->Print();   //prints the SD_Hit Print() output
 	(*fHitsCollection)[i]->Draw(); 
 	if ((*fHitsCollection)[i]->GetDetectorID() < 100) {
	dE += (*fHitsCollection)[i]->GetEnergy();
    } else {
	E += (*fHitsCollection)[i]->GetEnergy();
    }
  }
 	
std::cout << "dE: " <<dE<< " & E: "<<E <<std::endl;
G4AnalysisManager *analMan = G4AnalysisManager::Instance();analMan->FillH2(0,E, dE);
analMan->FillNtupleDColumn(0, E);
analMan->FillNtupleDColumn(1, dE );
analMan->AddNtupleRow();


}


/*G4bool NaI_SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    G4Track* track = step->GetTrack();
     G4double energyDeposit = step->GetTotalEnergyDeposit();
  if (energyDeposit == 0.0) return false;
   // G4double energy = track->GetKineticEnergy();
   //new G4cout << "Detected energy: " << energy / MeV << " MeV" << G4endl;
 
    
    //TODO :  Add Whatall information you want to get from steps
    
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
  G4int detectorID = step->GetTrack()->GetVolume()->GetCopyNo();  // Strip number
	
   /*	G4UserEventAction *evAction = const_cast<G4UserEventAction *>(G4RunManager::GetRunManager()->GetUserEventAction());
  NaI_EventAction *myEvAction = static_cast<NaI_EventAction *>(evAction);
 	G4String volumeName = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

     G4double edep = step->GetTotalEnergyDeposit();
    if (edep == 0.) return false;
    
    if (volumeName == "logicSiCD_dE") {
  myEvAction-> AddEdep_dE(edep);
  }
  
  if (volumeName == "logicSiCD_E") {
  myEvAction-> AddEdep_E(edep);
  }
    
    
  new G4cout<< "edep by SD" << edep << G4endl; (removed star/ from here)
    
    // Create a new hit for this event and strip
  SensitiveDetectorHit* hit = new SensitiveDetectorHit();
  hit->SetDetectorID(detectorID);
  hit->SetPosition(position);
  hit->SetEnergy(energyDeposit);

    
    // Store the hit in the collection
  fHitsCollection->insert(hit); 
  return true;
}
    
   // End of event: Print out the hits for each strip
void SensitiveDetector::EndOfEvent(G4HCofThisEvent* hce) {
  PrintHits();
} 

// Print hits per strip
void SensitiveDetector::PrintHits() {
  // Store a count of hits for each strip (detectorID)
  std::map<G4int, G4int> hitCountMap;
  
  // Iterate through all hits
  for (G4int i = 0; i < fHitsCollection->GetSize(); i++) {
    SensitiveDetectorHit* hit = (*fHitsCollection)[i];
    G4int stripNumber = hit->GetDetectorID();
    
    // Increment the hit count for the strip
    hitCountMap[stripNumber]++;
  }

  // Output the number of hits for each strip
  for (auto& entry : hitCountMap) {
    G4cout << "Strip " << entry.first << " has " << entry.second << " hits." << G4endl;
  }
}

    
    
    //New
 /*  if (track->GetDefinition() == G4Alpha::Definition()) {
        G4double edep = step->GetTotalEnergyDeposit();
        G4int trackID = track->GetTrackID();
        // You can store trackID in a set to avoid double-counting
        // or increment a counter
        new G4cout << "Alpha entered silicon. TrackID: " << trackID << ", Edep: " << edep << G4endl;
    }*/
   // return true;
//}*/
