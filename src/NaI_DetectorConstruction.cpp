#include "NaI_DetectorConstruction.h"
#include "NaI_SensitiveDetector.h"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

NaI_DetectorConstruction::NaI_DetectorConstruction() {}

NaI_DetectorConstruction::~NaI_DetectorConstruction() {}


void  NaI_DetectorConstruction::BuildSiCDDetector(G4Material* material, G4LogicalVolume* logicWorld,
                       G4double zPos_front,G4double zPos_back, G4double thickness,
                       G4int copyOffset, G4Colour colour,G4Color colour2,
                       G4String baseName, G4VSensitiveDetector* detector){
                       
                    
    for (int q = 0; q < 4; ++q) {
        G4double phiStart = q * 0.5 * CLHEP::pi;  // 0, π/2, π, 3π/2
        
          G4RotationMatrix* rotZ = new G4RotationMatrix();
            rotZ->rotateZ(phiStart);

        for (int i = 0; i < 16; ++i) {
        
		std::string stripName="PhysicalSiCD-dE_"+std::to_string(q)+"_";
    		stripName += std::to_string(i);
    		
            double rad_in = 2.4 + 0.15 * i;          // cm
            double rad_out = 2.4 + 0.15 * (i + 1);   // cm
            G4double height = thickness / 2;
            

            G4Tubs* solid = new G4Tubs(baseName + "_Front", rad_in * cm, rad_out * cm, height, 0., 0.5 * M_PI);
            G4LogicalVolume* logic = new G4LogicalVolume(solid, material, "Logic_" + baseName + "_Front_");
            
            logic->SetSensitiveDetector(detector);

          

            new G4PVPlacement(rotZ,
                              G4ThreeVector(0, 0, zPos_front),
                              logic,
                              stripName + "_Front",
                              logicWorld,
                              false,
                              copyOffset + q * 16 + i,
                              true);

            G4VisAttributes* vis = new G4VisAttributes(colour);
            vis->SetVisibility(true);
            vis->SetForceSolid(true);
            logic->SetVisAttributes(vis);
        }
    }

    // Back (4 wedges per quadrant)
    for (int p = 0; p < 16; ++p) {
        G4double phiStart = (p * 22.5)*deg;
        
         std::string BstripName="PhysicalSiCD-E_"+std::to_string(p);
    	//stripName += std::to_string(i);
        
        
            G4RotationMatrix* rotZ = new G4RotationMatrix();
            rotZ->rotateZ(phiStart);

        //for (int k = 0; k < 4; ++k) {
            double rad_in = 2.4;  // cm
            double rad_out = 4.8; // cm
            G4double height = thickness / 2;

            //double theta_min = 22.5 * k; // degrees
            double theta_min = 0;
            double theta_span = 22.5;

            G4Tubs* solid = new G4Tubs(baseName + "_Back",
                                       rad_in * cm, rad_out * cm, height,
                                       theta_min * deg, theta_span * deg);

            G4LogicalVolume* logic = new G4LogicalVolume(solid, material, "Logic_" + baseName + "_Back");
             
            logic->SetSensitiveDetector(detector);


            new G4PVPlacement(rotZ,
                              G4ThreeVector(0, 0, zPos_back),  // slightly behind front
                              logic,
                              BstripName+ "_Back",
                              logicWorld,
                              false,
                              copyOffset + p ,  // continue copy number from front
                              true);

            G4VisAttributes* vis = new G4VisAttributes(colour2);
            vis->SetVisibility(true);
            vis->SetForceSolid(true);
            logic->SetVisAttributes(vis);
        //}
    }
}
/*

    G4LogicalVolume* NaI_DetectorConstruction::GetFrontDetector(unsigned int numOfStrips,double thickness,G4Color color){
    G4NistManager *nist = G4NistManager::Instance();
    G4Material *envMaterial = nist->FindOrBuildMaterial("G4_Galactic");
      G4Material *Si = nist->FindOrBuildMaterial("G4_Si");
      //double halfZ = 0.00325*cm/2.+0.0000001;
      double halfZ = thickness/2.+0.0000001;
    G4Tubs *envelopeTube = new G4Tubs("EnvelopTube",0.,5*cm,halfZ,0,2*M_PI);
    G4LogicalVolume *envLogical = new G4LogicalVolume(envelopeTube,envMaterial,"logicalEnvelopeTube");
    
    
    
    for (int q = 0; q < 4; ++q) {
    	G4double phiStart = q * 0.5 * CLHEP::pi;  // 0, π/2, π, 3π/2
	
	G4RotationMatrix* rotZ = new G4RotationMatrix();
        rotZ->rotateZ(phiStart);
    // Loop over 16 radial strips in each quadrant
    
    for (int i = 0; i < numOfStrips; ++i) {
    std::string stripName="PhysicalSiCD-dE_"+std::to_string(q)+"_";
    	stripName += std::to_string(i);
    	std::cout << "RAMAN : " << stripName << std::endl;
        double rad_in = 2.4 + 0.15 * i;          // inner radius in cm
        double rad_out = 2.4 + 0.15 * (i + 1);    // outer radius in cm
        //double height_dE = 0.0325 * mm;
        double height_dE = thickness;

        G4Tubs* solidSiCD_dE = new G4Tubs("SiCD_dE", rad_in * cm, rad_out * cm, height_dE/2, 0., 0.5 * M_PI);
        G4LogicalVolume* logicSiCD_dE = new G4LogicalVolume(solidSiCD_dE, Si, "LogicalSiCD_dE");

       
        new G4PVPlacement(rotZ,
                          //G4ThreeVector(0* cm, 0 * cm, 12*cm),
                          G4ThreeVector(),
                          logicSiCD_dE,
                          stripName,
                          envLogical,
                          false,
                          q * 16 + i,  // unique copy number
                          true);
                          
        G4VisAttributes* redVis = new G4VisAttributes(color);  // RGB: Red
	redVis->SetVisibility(true);
	redVis->SetForceSolid(true);  // Show as solid (not wireframe)
	logicSiCD_dE->SetVisAttributes(redVis);
	
	NaI_SensitiveDetector* detector = new NaI_SensitiveDetector("SensitiveDetector_"+stripName);
	G4SDManager::GetSDMpointer()->AddNewDetector(detector);
        logicSiCD_dE->SetSensitiveDetector(detector);
	
	
    }
}

return envLogical;
    }*/

G4VPhysicalVolume *NaI_DetectorConstruction::Construct() {
  G4NistManager *nist = G4NistManager::Instance();
  G4Material *worldMat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material *Si = nist->FindOrBuildMaterial("G4_Si");

  // Modify the world volume dimension as required
  G4Box *solidWorld = new G4Box("World", 0.5 * m, 0.5 * m, 0.5 * m);
  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
  G4VPhysicalVolume *physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0);

 //---------------------------------------------------------------
 
 
 NaI_SensitiveDetector* sd = new NaI_SensitiveDetector("SiCD_SD");
 G4SDManager::GetSDMpointer()->AddNewDetector(sd);

// Build dE (thin front layer)
BuildSiCDDetector(Si, logicWorld,
                  12.0 * cm,
                  (12.+0.00325+0.000001)*cm,           // z-position
                  0.0325 * mm,          // thickness (two sides of 0.0325 mm)
                  0,                   // copy number offset
                  G4Colour::Red(),
                  G4Colour::Yellow(),     // color
                  "dE",sd);
                 
 // NaI_SensitiveDetector* sd = new NaI_SensitiveDetector("SiCD_SD");
 //G4SDManager::GetSDMpointer()->AddNewDetector(sd_E);

//2nd detector E
BuildSiCDDetector(Si, logicWorld,
                  12.676625 * cm,
                   12.726625*cm,    // z-position (just behind dE)
                  0.05* mm,            // thicker
                  100,                 // copy number offset to avoid conflicts
                  G4Colour::Blue(),
                  G4Colour::White(),    // color
                  "E",sd);


//Front for dE
//;G4LogicalVolume *frontDet = GetFrontDetector(16,0.00325*cm,G4Colour(0., 1.0, 0.0));

/*;new G4PVPlacement(0, G4ThreeVector(0* cm, 0 * cm, 12*cm),
                          frontDet,
                          "FrontEnvelope-dE",
                          logicWorld,
                          false,
                          0,  // unique copy number
                          true);
                          
  G4LogicalVolume *frontDet = GetFrontDetector(16,0.05*cm,G4Colour(0., .0, 1.0));

new G4PVPlacement(0,
                          G4ThreeVector(0* cm, 0 * cm, 17*cm),
                          frontDet,
                          "FrontEnvelope-E",
                          logicWorld,
                          false,
                          0,  // unique copy number
                          true);        */                

















  // Logic to Attach sensitive detector to a logical volume
 // NaI_SensitiveDetector* detector = new NaI_SensitiveDetector("SensitiveDetector");
//  G4SDManager::GetSDMpointer()->AddNewDetector(detector);
  // logicWorld->SetSensitiveDetector(detector);
  //logicSiCD_E_f->SetSensitiveDetector(detector);
  ///logicSiCD_E_bk->SetSensitiveDetector(detector);
  //logicSiCD_dE_f->SetSensitiveDetector(detector);
  //logicSiCD_dE_f->SetSensitiveDetector(detector);


  return physWorld;
}
