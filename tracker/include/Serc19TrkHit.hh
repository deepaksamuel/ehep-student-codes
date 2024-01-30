//

#ifndef Serc19TrkHit_h
#define Serc19TrkHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class Serc19TrkHit : public G4VHit
{
  public:

      Serc19TrkHit();
      ~Serc19TrkHit();
      Serc19TrkHit(const Serc19TrkHit &right);
      const Serc19TrkHit& operator=(const Serc19TrkHit &right);
      G4int operator==(const Serc19TrkHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4int    pdgid;  //Particle ID
      G4ThreeVector pos;
      G4double localthe;
      G4double localphi;
      G4ThreeVector mom; //Momentum of track at earliest energy deposite, has meaning only for muon track, not usefull at all for hadronic shower
      G4double toff;
      G4double tofx;
      G4double tofy;
      unsigned int    HitId;
  public:
      inline void SetpdgId(G4int id)
      { pdgid = id; }
      inline void SetEdep(G4double de)
      { edep = de; }
      inline void AddEdep(G4double de) 
      { edep +=de; }
      inline G4double GetEdep()
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos()
      { return pos; }
      inline void SetMom(G4ThreeVector xyz)
      { mom = xyz; }
      inline G4ThreeVector GetMom()
      { return mom; }
      inline void SetTime(G4double tf)
      { toff = tf; }
      inline G4double GetTime()
      { return toff; }

      inline void SetLocalTheta(G4double xyz)
     { localthe = xyz; }
      inline G4double GetLocalTheta()
      { return localthe; }
  
      inline void SetLocalPhi(G4double xyz)
      { localphi = xyz; }
      inline G4double GetLocalPhi()
      { return localphi; }
  

  //      inline void SetTimeX(G4double tf)
  //      { tofx = tf; }
  //      inline G4double GetTimeX()
  //      { return tofx; }
  //      inline void SetTimeY(G4double tf)
  //      { tofy = tf; }
  //     inline G4double GetTimeY()
  //      { return tofy; }
      inline void SetHitId (unsigned int id)
      { HitId = id; }
      inline unsigned int GetHitId()
      { return HitId; }
      inline G4int GetpdgId()
      { return pdgid; }
};

typedef G4THitsCollection<Serc19TrkHit> Serc19TrkHitsCollection;

extern G4Allocator<Serc19TrkHit> Serc19TrkHitAllocator;

inline void* Serc19TrkHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Serc19TrkHitAllocator.MallocSingle();
  return aHit;
}

inline void Serc19TrkHit::operator delete(void *aHit)
{
  Serc19TrkHitAllocator.FreeSingle((Serc19TrkHit*) aHit);
}

#endif
