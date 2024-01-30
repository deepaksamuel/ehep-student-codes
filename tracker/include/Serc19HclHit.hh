//

#ifndef Serc19HclHit_h
#define Serc19HclHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class Serc19HclHit : public G4VHit
{
  public:

      Serc19HclHit();
      ~Serc19HclHit();
      Serc19HclHit(const Serc19HclHit &right);
      const Serc19HclHit& operator=(const Serc19HclHit &right);
      G4int operator==(const Serc19HclHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4int toff;
  G4ThreeVector pos;
      unsigned int  HitId;
  public:
      inline void SetEdep(G4double de)
      { edep = de; }
      inline void AddEdep(G4double de) 
      { edep +=de; }
      inline G4double GetEdep()
      { return edep; }
      inline void SetTime(G4int tf)
      { toff = tf; }
      inline G4int GetTime()
      { return toff; }
      inline void SetHitId (unsigned int id)
      { HitId = id; }
      inline unsigned int GetHitId()
      { return HitId; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos()
      { return pos; }

};

typedef G4THitsCollection<Serc19HclHit> Serc19HclHitsCollection;

extern G4Allocator<Serc19HclHit> Serc19HclHitAllocator;

inline void* Serc19HclHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) Serc19HclHitAllocator.MallocSingle();
  return aHit;
}

inline void Serc19HclHit::operator delete(void *aHit)
{
  Serc19HclHitAllocator.FreeSingle((Serc19HclHit*) aHit);
}

#endif
