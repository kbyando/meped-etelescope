//
////////////////////////////////////////////////////////////////
//                                                            //
//  Title:  Physics List for Space Electronic Applications    //
//  Date:   20 March 2005                                     //
//  Author: D.H. Wright (SLAC)                                //
//                                                            //
////////////////////////////////////////////////////////////////
//

#ifndef SEPhysicsList_h
#define SEPhysicsList_h 1

#include "G4VModularPhysicsList.hh"


class SEPhysicsList: public G4VModularPhysicsList
{
public:
  SEPhysicsList();
  virtual ~SEPhysicsList();
  
  virtual void SetCuts();
};

#endif



