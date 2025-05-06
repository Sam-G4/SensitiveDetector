#ifndef NaI_EVENTACTION_HH
#define NaI_EVENTACTION_HH

#include "G4Event.hh"
#include "G4UserEventAction.hh"
class NaI_EventAction : public G4UserEventAction {
  double eDep_E;
  double eDep_dE;

public:
  NaI_EventAction();
  virtual ~NaI_EventAction();
  virtual void BeginOfEventAction(const G4Event *event);
  virtual void EndOfEventAction(const G4Event *event);
  //void AddEnergy(double val);  //$
  void AddEdep_dE(G4double edep); //{ eDep_dE += edep; } //new $
  void AddEdep_E(G4double edep); //{ eDep_E += edep; }   //new $
};

#endif
