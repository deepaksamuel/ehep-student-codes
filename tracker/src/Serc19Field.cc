//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: Serc19Field.cc 78055 2013-12-03 08:27:48Z gcosmo $
//
/// \file parallel/ParN04/src/Serc19Field.cc
/// \brief Implementation of the Serc19Field class
//

#include "Serc19Field.hh"
#include "G4SystemOfUnits.hh"

Serc19Field::Serc19Field()
{
  Bx = 2.0*tesla;
  rmax_sq = sqr(50.*cm);
  zmax = 100.*cm;
}

Serc19Field::~Serc19Field()
{;}

void Serc19Field::GetFieldValue(const double Point[3],double *Bfield) const
{
  Bfield[0] = Bx;
  Bfield[1] = By;
  Bfield[2] = Bz;

  //  if(std::abs(Point[2])<zmax && (sqr(Point[0])+sqr(Point[1]))<rmax_sq)
  //  { Bfield[2] = Bz; }
  //  else
  //  { Bfield[2] = 0.; }
}

void Serc19Field::SetFieldValues(double a, double b, double c) {
  Bx = a;
  By = b;
  Bz = c;
}
