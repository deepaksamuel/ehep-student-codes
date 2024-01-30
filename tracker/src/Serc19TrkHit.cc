
#include "Serc19TrkHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


G4Allocator<Serc19TrkHit> Serc19TrkHitAllocator;

Serc19TrkHit::Serc19TrkHit() {
  pdgid=-25;
  edep = 0;
  toff = 1000000;
  
}

Serc19TrkHit::~Serc19TrkHit()
{;}

Serc19TrkHit::Serc19TrkHit(const Serc19TrkHit &right)
  : G4VHit()
{
  pdgid  = right.pdgid;
  edep = right.edep;
  pos = right.pos;
  mom = right.mom;
  toff = right.toff;
  HitId = right.HitId;
}

const Serc19TrkHit& Serc19TrkHit::operator=(const Serc19TrkHit &right) {
  pdgid = right.pdgid;
  edep = right.edep;
  pos = right.pos;
  mom = right.mom;
  toff = right.toff;
  HitId = right.HitId;

  return *this;
}

G4int Serc19TrkHit::operator==(const Serc19TrkHit &right) const {
  return (this==&right) ? 1 : 0;
}

void Serc19TrkHit::Draw() {
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    G4Circle circle(pos);
    circle.SetScreenSize(2.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void Serc19TrkHit::Print()
{;}


