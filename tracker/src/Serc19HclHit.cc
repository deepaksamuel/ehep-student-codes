
#include "Serc19HclHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


G4Allocator<Serc19HclHit> Serc19HclHitAllocator;

Serc19HclHit::Serc19HclHit()
{;}

Serc19HclHit::~Serc19HclHit()
{;}

Serc19HclHit::Serc19HclHit(const Serc19HclHit &right)
  : G4VHit()
{
  edep = right.edep;
  pos = right.pos;
  toff = right.toff;
  HitId = right.HitId;
}

const Serc19HclHit& Serc19HclHit::operator=(const Serc19HclHit &right) {
  edep = right.edep;
  pos = right.pos;
  toff = right.toff;
  HitId = right.HitId;
  
  return *this;
}

G4int Serc19HclHit::operator==(const Serc19HclHit &right) const {
  return (this==&right) ? 1 : 0;
}

void Serc19HclHit::Draw() {
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  //  if (!pVVisManager->FilterHit(*this)) return;
  if(pVVisManager) {
    G4Circle circle(pos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.,0.,1.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void Serc19HclHit::Print()
{;}


