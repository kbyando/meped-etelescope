// NOAA POES Monte Carlo Simulation v1.3, 21/10/2008e
// MEPED :: Electron Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MEPED_DetectorConstruction_h
#define MEPED_DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"

class G4Box;
class G4Tubs;
class G4Cons;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
// G4VisAttributes (perhaps to be implemented someday)
class MEPED_DetectorMessenger;

class MEPED_DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    MEPED_DetectorConstruction();
   ~MEPED_DetectorConstruction();

  public:
     
     void SetDetectorMaterial 	(G4String);     
     void SetShieldingMaterial 	(G4String);
     void SetHousingMaterial 	(G4String);
     void SetContactMaterial 	(G4String);    
	void SetD1CoatingMaterial	(G4String);
	void SetEpoxyMaterial		(G4String); 
     void SetBelvilleMaterial	(G4String);
     void SetInsulatorMaterial 	(G4String);
     void SetFoilMaterial	(G4String);

     void SetDetectorThickness  (G4double);   
      
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintModuleParameters(); 
                    
     G4Material* GetDetectorMaterial() {return DetectorMaterial;};     
     G4Material* GetShieldingMaterial()     {return ShieldingMaterial;};
     G4Material* GetHousingMaterial()     {return HousingMaterial;};

     G4Material* GetContactMaterial()     {return ContactMaterial;};
     G4Material* GetD1CoatingMaterial()	{return D1CoatingMaterial;};
     G4Material* GetEpoxyMaterial() 	{return EpoxyMaterial;};

     G4Material* GetBelvilleMaterial()     {return BelvilleMaterial;};
     G4Material* GetInsulatorMaterial()     {return InsulatorMaterial;};
     G4Material* GetFoilMaterial()	{return FoilMaterial;};

     G4double GetDetectorThickness() {return DetectorThickness;};

     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetDete1()  {return physi_D1ActiveLayer;};
                
  private:
     
     G4Material*        DetectorMaterial;
     G4Material*	ShieldingMaterial;
     G4Material*        HousingMaterial;
     G4Material*        ContactMaterial;
     G4Material*	D1CoatingMaterial; // Al
     G4Material*	EpoxyMaterial;
     G4Material*        BelvilleMaterial;
     G4Material*        InsulatorMaterial;
     G4Material*	FoilMaterial; // Ni

     G4Material*	defaultMaterial;

// Detector Characteristics
     G4double           DetectorRadius;     
     G4double		DetectorThickness;
     G4double  		DActiveLayer;
     G4double		D1Silicon_rad;
// External Specifications
     G4double		etAssembly_rad;
     G4double  		etAssembly_length;
// Miscellaneous Components
     G4double		outerWall_thickness;
     G4double		partition_thickness;
     G4double		shield_thickness;
     G4double		frontCover_thickness;
     G4double		rearCover_thickness;
// Detector Housing
     G4double 		DeteHousing_length;
     G4double		DeteHousing_dia;
     G4double		inner_DH_dia;
// Inner Housing
     G4double 		BelvilleThickness;
     G4double 		ContactThickness;
     G4double		ContactDeposit;
     G4double		FoilDeposit;
     G4double		AlDensity;
     G4double		NiDensity;

     G4double		BelvilleInnerRadius;
     G4double		BelvilleOuterRadius;
     G4double		ContactInnerRadius;
     G4double		ContactOuterRadius;

     G4double 		InsulatorThickness;

// Collimator Characteristics (revised spec)
	G4double	fullWidthAngle;
	G4double	fullWidthICAngle;

	G4double	DIA1; // Front Cover Aperture
	G4double	DIA2inner; // Outer Collimator
	G4double	DIA2outer;
	G4double	DIA3inner; //	 "	"
	G4double	DIA3outer;
	G4double	DIA4;	// Structural Aperture
	G4double	DIA5;	// Front Shield Aperture
	G4double	DIA6;	// 	" " "

	G4double	AlColl1_dia; // Aluminium Collimator 1
	G4double	AlColl2_dia; // 	" " 2
	G4double	AlColl3_dia; //		" " 3
	G4double	AlColl_odia; // Outer Diameter for above

	G4double	CollimatorSep;
	G4double	CollimatorDia;
	G4double	CollimatorHeight;	

     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     G4Tubs*             solid_etAssembly;    //pointer to the solid assembly 
     G4LogicalVolume*   logic_etAssembly;    //pointer to the logical assembly
     G4VPhysicalVolume* physi_etAssembly;    //pointer to the physical assembly

	G4Tubs*			solidShellStruct;
	G4LogicalVolume*	logicShellStruct;
	G4VPhysicalVolume*	physiShellStruct;
	
     G4Tubs*             solidCollimatorChamber;    //pointer..
     G4LogicalVolume*   logicCollimatorChamber;    //pointer to the...
     G4VPhysicalVolume* physiCollimatorChamber;    //pointer to the//
         
     G4Cons*             solidCollimator; //pointer..
     G4LogicalVolume*   logicCollimator; //pointer..r
     G4VPhysicalVolume* physiCollimator; //pointer..

     G4Cons*		solidOCOuterCone;
     G4LogicalVolume*	logicOCOuterCone;
     G4VPhysicalVolume*	physiOCOuterCone;

     G4Cons*            solidOCInnerCone;
     G4LogicalVolume*   logicOCInnerCone;
     G4VPhysicalVolume* physiOCInnerCone;

     G4Tubs*		solidBChamber;
     G4LogicalVolume*	logicBChamber;
     G4VPhysicalVolume*	physiBChamber;

	G4Tubs*			solidSpacer;
	G4LogicalVolume*	logicSpacer;
	G4VPhysicalVolume*	physiSpacer;

// Misc. Shielding
	G4Tubs*			solidFrontCover;
	G4LogicalVolume*	logicFrontCover;
	G4VPhysicalVolume*	physiFrontCover;
	
	G4Cons*			solidPartition;
	G4LogicalVolume*	logicPartition;
	G4VPhysicalVolume*	physiPartition;

	G4Tubs*			solidRearCover;
	G4LogicalVolume*	logicRearCover;
	G4VPhysicalVolume*	physiRearCover;

// Detector Sub-Housing
     G4Tubs*		solidDeteHousing;
     G4LogicalVolume*	logicDeteHousing;
     G4VPhysicalVolume*	physiDeteHousing;

     G4Cons*		solidInnerCollimator;
     G4LogicalVolume*	logicInnerCollimator;
     G4VPhysicalVolume*	physiInnerCollimator;

     G4Tubs*		solidInnerHousing;
     G4LogicalVolume*   logicInnerHousing;
     G4VPhysicalVolume* physiInnerHousing;

     G4Tubs*            solidDStack;
     G4LogicalVolume*   logicDStack;
     G4VPhysicalVolume* physiDStack;

     G4Tubs*            solidBelville1;
     G4LogicalVolume*   logicBelville1;
     G4VPhysicalVolume* physiBelville1;

     G4Tubs*            solidContact1;
     G4LogicalVolume*   logicContact1;
     G4VPhysicalVolume* physiContact1;

     G4Tubs*            solid_Detector1;
     G4LogicalVolume*   logic_Detector1;
     G4VPhysicalVolume* physi_Detector1;

	G4Tubs*			solid_D1RingMntA;
	G4LogicalVolume*	logic_D1RingMntA;
	G4VPhysicalVolume*	physi_D1RingMntA;

	G4Tubs*			solid_D1RingMntB;
	G4LogicalVolume*	logic_D1RingMntB;
	G4VPhysicalVolume*	physi_D1RingMntB;

	G4Tubs*			solid_D1SurfDeposit;
	G4LogicalVolume*	logic_D1SurfDeposit;
	G4VPhysicalVolume*	physi_D1SurfDeposit;

     G4Tubs*		solid_D1ActiveLayer;
     G4LogicalVolume*	logic_D1ActiveLayer;
     G4VPhysicalVolume*	physi_D1ActiveLayer;

     G4Tubs*		solidContact2;
     G4LogicalVolume*   logicContact2;
     G4VPhysicalVolume* physiContact2;

     G4Tubs*		solidBelville2;
     G4LogicalVolume*	logicBelville2;
     G4VPhysicalVolume*	physiBelville2;

     G4Tubs*		solidDInsulator; // B. Detector Insulator
     G4LogicalVolume*	logicDInsulator; //  (also "Detector Holder")
     G4VPhysicalVolume*	physiDInsulator;

     G4Tubs*		solidFrontInsul; // i. Front
     G4LogicalVolume*	logicFrontInsul;
     G4VPhysicalVolume*	physiFrontInsul;

     G4Tubs*		solidLongiInsul; // ii. Longitudinal
     G4LogicalVolume*	logicLongiInsul;
     G4VPhysicalVolume*	physiLongiInsul;

     G4Tubs*		solidBackInsul; // iii. Back
     G4LogicalVolume*	logicBackInsul;
     G4VPhysicalVolume*	physiBackInsul;

     G4Tubs*		solidHouseShielding;
     G4LogicalVolume*	logicHouseShielding;
     G4VPhysicalVolume* physiHouseShielding;

     G4Tubs*		solidNiFoil;
     G4LogicalVolume*	logicNiFoil;
     G4VPhysicalVolume*	physiNiFoil;

	G4Tubs*			solidBackShield;
	G4LogicalVolume*	logicBackShield;
	G4VPhysicalVolume*	physiBackShield;

     G4RotationMatrix*	AssemblyRotation;
     G4double 		etAssemblyAngle;

     MEPED_DetectorMessenger* detectorMessenger;  //pointer to the Messenger
      
  private:
    
     void DefineMaterials();
     void ComputeModuleParameters();
     G4VPhysicalVolume* ConstructModule();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void MEPED_DetectorConstruction::ComputeModuleParameters()
{
  // Compute derived parameters of the calorimeter
  // (Empty) 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

