
#include "Serc19EcalHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


G4Allocator<Serc19EcalHit> Serc19EcalHitAllocator;

Serc19EcalHit::Serc19EcalHit()
{;}

Serc19EcalHit::~Serc19EcalHit()
{;}

Serc19EcalHit::Serc19EcalHit(const Serc19EcalHit &right)
  : G4VHit()
{
  edep = right.edep;
  pos = right.pos;
  toff = right.toff;
  HitId = right.HitId;
}

const Serc19EcalHit& Serc19EcalHit::operator=(const Serc19EcalHit &right) {
  edep = right.edep;
  pos = right.pos;
  toff = right.toff;
  HitId = right.HitId;
  
  return *this;
}

G4int Serc19EcalHit::operator==(const Serc19EcalHit &right) const {
  return (this==&right) ? 1 : 0;
}

void Serc19EcalHit::Draw() {
//   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
//   if(pVVisManager) {
//     G4Circle circle(pos);
//     circle.SetScreenSize(0.04);
//     circle.SetFillStyle(G4Circle::filled);
//     G4Colour colour(0.,1.,0.);
//     G4VisAttributes attribs(colour);
//     circle.SetVisAttributes(attribs);
//     pVVisManager->Draw(circle);
//   }
}

void Serc19EcalHit::Print()
{;}


