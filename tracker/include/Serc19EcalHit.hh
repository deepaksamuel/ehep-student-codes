//

#ifndef Serc19EcalHit_h
#define Serc19EcalHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class Serc19EcalHit : public G4VHit
{
  public:

      Serc19EcalHit();
      ~Serc19EcalHit();
      Serc19EcalHit(const Serc19EcalHit &right);
      const Serc19EcalHit& operator=(const Serc19EcalHit &right);
      G4int operator==(const Serc19EcalHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;
      G4double toff;
      unsigned int  HitId;
  public:
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
      inline void SetTime(G4double tf)
      { toff = tf; }
      inline G4double GetTime()
      { return toff; }
      inline void SetHitId (unsigned int id)
      { HitId = id; }
      inline unsigned int GetHitId()
      { return HitId; }
};

typedef G4THitsCollection<Serc19EcalHit> Serc19EcalHitsCollection;

extern G4Allocator<Serc19EcalHit> Serc19EcalHitAllocator;

inline void* Serc19EcalHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Serc19EcalHitAllocator.MallocSingle();
  return aHit;
}

inline void Serc19EcalHit::operator delete(void *aHit)
{
  Serc19EcalHitAllocator.FreeSingle((Serc19EcalHit*) aHit);
}

#endif
