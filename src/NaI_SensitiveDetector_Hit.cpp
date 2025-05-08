#include "NaI_SensitiveDetector_Hit.h"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"


//Constructor
NaI_SensitiveDetector_Hit::NaI_SensitiveDetector_Hit(): fDetectorID(-1), fPosition(G4ThreeVector()), fEnergy(0.0){}

//Destructor
NaI_SensitiveDetector_Hit::~NaI_SensitiveDetector_Hit() {}


//Setters and Getters
void NaI_SensitiveDetector_Hit::SetDetectorID(G4int id) {
    fDetectorID = id;
}

G4int NaI_SensitiveDetector_Hit::GetDetectorID() const {
    return fDetectorID;
}

void NaI_SensitiveDetector_Hit::SetPosition(const G4ThreeVector& position) {
    fPosition = position;
}

G4ThreeVector NaI_SensitiveDetector_Hit::GetPosition() const {
    return fPosition;
}

void NaI_SensitiveDetector_Hit::SetEnergy(G4double energy) {
    fEnergy = energy;
}

G4double NaI_SensitiveDetector_Hit::GetEnergy() const {
    return fEnergy;
}


void NaI_SensitiveDetector_Hit::SetParticleName(const G4String& name){
	fParticle = name;
}
G4String NaI_SensitiveDetector_Hit::GetParticleName() const {
	return fParticle;
}

//Visualization
void NaI_SensitiveDetector_Hit::Draw() {
G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();
    if (!visManager) return;

    G4Circle circle(fPosition);
    circle.SetScreenSize(4.0);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.0, 0.0, 0.0); // red
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);

    visManager->Draw(circle);
    }

void NaI_SensitiveDetector_Hit::Print() {
  /*  G4cout << "Hit: "
           << "Detector ID = " << fDetectorID
           << ", Energy = " << G4BestUnit(fEnergy, "Energy")
           << ", Position = " << fPosition
            << ", Particle = " << fParticle
           << G4endl;*/
}

