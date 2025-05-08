#ifndef NaI_SENSITIVEDETECTOR_HIT_HH
#define NaI_SENSITIVEDETECTOR_HIT_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"

class NaI_SensitiveDetector_Hit : public G4VHit {
private:
	 G4int fDetectorID; // The detector's copy number or ID
    	 G4ThreeVector fPosition; // Position of the hit
    	 G4double fEnergy; // Energy deposited in the detector
    	 G4String fParticle; 

public:
    NaI_SensitiveDetector_Hit();
    virtual ~NaI_SensitiveDetector_Hit();
    
     // Setters and getters
    void SetDetectorID(G4int id);
    G4int GetDetectorID() const;

    void SetPosition(const G4ThreeVector& position);
    G4ThreeVector GetPosition() const;

    void SetEnergy(G4double energy);
    G4double GetEnergy() const;
    
    void SetParticleName(const G4String& name);
    G4String GetParticleName() const;

    // Override virtual methods for hit visualization and printing
    void Draw() override;
    void Print() override;
    
    //int fChannelNo;
    //double fEnergy;
};

#endif
