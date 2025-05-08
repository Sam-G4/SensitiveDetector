#ifndef NaI_DETECTORCONSTRUCTION_HH
#define NaI_DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4Colour.hh"
#include "globals.hh"

class NaI_DetectorConstruction : public G4VUserDetectorConstruction {
public:
    NaI_DetectorConstruction();
    virtual ~NaI_DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
   // G4LogicalVolume* GetFrontDetector(unsigned int numOfStrips,double thickness,G4Color color);
   void BuildSiCDDetector(G4Material* material, G4LogicalVolume* logicWorld,
                           G4double zPos_front,G4double zPos_back, G4double thickness,
                           G4int copyOffset, G4Colour colour,G4Colour colour2,
                           G4String baseName, G4VSensitiveDetector* detector);

};

#endif
