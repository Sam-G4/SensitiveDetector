#ifndef NaI_SENSITIVEDETECTOR_HH
#define NaI_SENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "NaI_SensitiveDetector_Hit.h"
#include "G4THitsCollection.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
//#include "G4VHit.hh"

//#include <set>

class NaI_SensitiveDetector : public G4VSensitiveDetector {
public:
    NaI_SensitiveDetector(const G4String& name);
    virtual ~NaI_SensitiveDetector();
    
     virtual void Initialize(G4HCofThisEvent* hce) override;
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;  //New override, was not there b4
    //void Reset();  //New
    virtual void EndOfEvent(G4HCofThisEvent* hce) override;
    // size_t GetAlphaCount() const { return alphaTrackIDs.size(); } //New
//private: //New
    //std::set<G4int> alphaTrackIDs;  //New
    //void PrintHits();
    
    // Retrieve hit collection(optional)
    NaI_SensitiveDetector_Hit* GetHit(G4int i);

 private:
    // A collection of hits for this detector
    G4THitsCollection<NaI_SensitiveDetector_Hit>* fHitsCollection;

    // Hits collection ID
    static G4int fCollectionID;   
    
    
};

#endif
