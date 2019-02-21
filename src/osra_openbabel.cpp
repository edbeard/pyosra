/******************************************************************************
 OSRA: Optical Structure Recognition Application

 Created by Igor Filippov, 2007-2013 (igor.v.filippov@gmail.com)

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110-1301, USA
 *****************************************************************************/

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>
#include <openbabel/alias.h>
#include <openbabel/stereo/tetrahedral.h>

#include <sstream> // std:ostringstream
#include <iostream> // std::cerr

#include "osra_common.h" // trim()
#include "osra_structure.h"
#include "osra_openbabel.h"
#include "osra.h"
#include "osra_stl.h"
#include "mcdlutil.h"

using namespace OpenBabel;

#define HYDROGEN_ATOMIC_NUM     1
#define LITHIUM_ATOMIC_NUM      3
#define CARBON_ATOMIC_NUM       6
#define OXYGEN_ATOMIC_NUM       8
#define FLUORINE_ATOMIC_NUM     9
#define SILICONE_ATOMIC_NUM  14
#define CHLORINE_ATOMIC_NUM     17
#define ARGON_ATOMIC_NUM        18
#define BROMINE_ATOMIC_NUM      35
#define IODINE_ATOMIC_NUM       53
#define URANIUM_ATOMIC_NUM 92

// Look at this issue: https://sourceforge.net/tracker/?func=detail&aid=3425216&group_id=40728&atid=428740
#define AROMATIC_BOND_ORDER     5

int osra_openbabel_init()
{
  OBConversion conv;
  OBMol mol;

  conv.SetInFormat("SMI");
  conv.ReadString(&mol, "[*]");

  if (mol.NumAtoms() == 0)
    return ERROR_UNKNOWN_OPENBABEL_FORMAT;

  return 0;
}

// Function: create_atom()
//
// For the atom represented by its OCR'ed label create a new atom in the given molecule.
//
// Parameters:
//      mol - the current molecule
//      atom - atom to add to molecule (will be updated)
//      scale - scale factor / coordinates multiplier
//      superatom - superatom dictionary, that maps atom labels to corresponding SMILES
//      verbose - print debug information
//
// Returns:
//      true in case the given atom is superatom
//
bool create_atom(OBMol &mol, atom_t &atom, double scale,
                 const std::map<std::string, std::string> &superatom, bool verbose)
{
  if (atom.label.empty() || atom.label == " ")
    {
      atom.anum = CARBON_ATOMIC_NUM;
    }
  else
    {
      // Lookup in superatom dictionary:
      std::map<std::string, std::string>::const_iterator it = superatom.find(atom.label);

      if (it != superatom.end())
        {
          // "superatom" case (e.g. "COOH")
          const std::string &smiles_superatom = it->second;

          OBConversion conv;
          OBMol superatom_mol;

          conv.SetInFormat("SMI");
          conv.ReadString(&superatom_mol, smiles_superatom);

          if (verbose)
            std::cout << "Considering superatom " << atom.label << "->" << smiles_superatom <<
                " vector: " << atom.x * scale << "x" << -atom.y * scale << '.' << std::endl;

          // This is the index of first atom in superatom in molecule:
          atom.n = mol.NumAtoms() + 1;

          OBAtomIterator atom_iter;

          // Transfer all atoms from "superatom" molecule to current molecule.


          for (OBAtom *a = superatom_mol.BeginAtom(atom_iter); a; a = superatom_mol.NextAtom(atom_iter))
            {
              if (verbose)
                std::cout << "Adding atom #" << mol.NumAtoms() + 1 << ", anum: " << a->GetAtomicNum() << std::endl;

              OBAtom *new_atom = mol.NewAtom();

              new_atom->SetAtomicNum(a->GetAtomicNum());
              new_atom->SetFormalCharge(a->GetFormalCharge());
	      if (!atom.label.empty() && atom.label != " ")
		{
		  // Unknown atom?
		  OBPairData *label = new OBPairData;
		  label->SetAttribute("UserLabel");
		  label->SetValue(atom.label);
		  label->SetOrigin(userInput); // set by user, not by Open Babel
		  new_atom->SetData(label);
		}
	      if (atom.anum == 0)
		{
		  AliasData* ad = new AliasData();
		  ad->SetAlias(atom.label);
		  ad->SetOrigin(external);
		  new_atom->SetData(ad);
		}
            }

          // Correct first atom meta-info:
          OBAtom *first_superatom = mol.GetAtom(atom.n);

          first_superatom->SetVector(atom.x * scale, -atom.y * scale, 0);

          atom.anum = first_superatom->GetAtomicNum();


          int first_bond_index = mol.NumBonds();

          OBBondIterator bond_iter;

          // Transfer all bonds from "superatom" molecule to current molecule:
          for (OBBond *b = superatom_mol.BeginBond(bond_iter); b; b = superatom_mol.NextBond(bond_iter))
            {
              if (verbose)
                std::cout << "Adding bond #" << mol.NumBonds() << " " << b->GetBeginAtomIdx() + atom.n - 1 << "->" << b->GetEndAtomIdx() + atom.n - 1
                          << ", order: " << b->GetBondOrder() << ", flags: " << b->GetFlags() << '.' << std::endl;

              mol.AddBond(b->GetBeginAtomIdx() + atom.n - 1, b->GetEndAtomIdx() + atom.n - 1,
                          b->GetBondOrder(), b->GetFlags());
            }

          // If at least one bond was added, the "superatom" coordinates should be recalculated:
          return first_bond_index != mol.NumBonds();
        }

      // If not found, lookup the atom number in periodic table of elements:
      atom.anum = etab.GetAtomicNum(atom.label.c_str());
    }

  atom.n = mol.NumAtoms() + 1;

  if (verbose)
    std::cout << "Creating atom #" << atom.n << " \"" << atom.label << "\", anum: " << atom.anum << '.' << std::endl;

  OBAtom *a = mol.NewAtom();

  a->SetAtomicNum(atom.anum);
  a->SetVector(atom.x * scale, -atom.y * scale, 0);

  if (atom.charge != 0)
    a->SetFormalCharge(atom.charge);

  if (atom.anum != 0 && !atom.label.empty() && atom.label != " ")
    {
      // Unknown atom?
      OBPairData *label = new OBPairData;
      label->SetAttribute("UserLabel");
      label->SetValue(atom.label);
      label->SetOrigin(userInput); // set by user, not by Open Babel
      a->SetData(label);
    }
  if (atom.anum == 0)
    {
      AliasData* ad = new AliasData();
      ad->SetAlias(atom.label);
      ad->SetOrigin(external); 
      a->SetData(ad);
    }

  return false;
}

// Function: confidence_function()
//
// Calculates confidence estimate based on molecular counts provided by <create_molecule()>
//
// Parameters:
//      C_Count, N_Count, O_Count, F_Count, S_Count, Cl_Count, Br_Count - number of carbon, nitrogen, oxygen, fluorine, sulfur, chlorine, and bromine atoms
//      R_Count - number of recognized Markush atomic labels, such as R1, R2....
//      Xx_Count - number of unrecognized atomic labels from <osra.cpp::remove_small_terminal_bonds()>
//      num_rings - number of rings
//      num_aromatic - number of aromatic rings
//      num_fragments - number of fragments
//      Num_Rings - vector of counts for number of 3,4,5,6,7-member rings
//
// Returns:
//      confidence estimate
double confidence_function(int* x, int n)
{
  double c[] = {-0.11469143725730054, 0.15723547931889853, 0.19765680222250673, 0.249101590474403, 0.1897669087341134, 0.19588348907301223, 0.3354622208036507, 0.16779269801176255, -0.21232000222198893, 0.016958281784354032, -0.08672059360133752, -0.05105752296619957, -0.349912750824004, 0.18836317536530647, 0.22316782354758827, 0.27741998968081166, 0.25710999274481955, 0.27968899280120096, 0.12695166847876285, -0.10020778884718293, 0.05150631410596443, 0.22283571763712148, 0.23130179826714167, 0.1049054095759948, 0.05333970810460394, -0.12491056666737535};
  double r = 0;
  for (int i=0; i<n; i++)
    r += c[i]*x[i];

  return (r);
}

void SetTetrahedtalUnknown(OBMolAtomIter atom)
{
  OBStereoFacade facade(atom->GetParent());
  OBTetrahedralStereo *stereo;
  if (facade.HasTetrahedralStereo(atom->GetId()))
    stereo = facade.GetTetrahedralStereo(atom->GetId());
  else
    stereo = new OBTetrahedralStereo(atom->GetParent());

  OBTetrahedralStereo::Config config = stereo->GetConfig();
  config.center = atom->GetId();
  config.specified = true;
  config.winding = OBStereo::UnknownWinding;
  //config.specified = false;
  stereo->SetConfig(config);

  if (!facade.HasTetrahedralStereo(atom->GetId()))
    atom->GetParent()->SetData(stereo);
}

// Function: create_molecule()
//
// Converts vectors of atoms and bonds into a molecular object and calculates the molecule statistics.
// Note: this function changes the atoms!
//
// Parameters:
//      atom - vector of <atom_s> atoms
//      bond - vector of <bond_s> bonds
//      n_bond - total number of bonds
//      avg_bond_length - average bond length as measured from the image (to be included into output if provided)
//      molecule_statistics - the molecule statistics (returned to the caller)
//      generate_2D_coordinates - generate 2D coordinates for chemical groups
//      confidence - confidence score (returned to the caller if provided)
//      superatom - dictionary of superatom labels mapped to SMILES
//      verbose - print debug info
void create_molecule(OBMol &mol, std::vector<atom_t> &atom,
                     const std::vector<bond_t> &bond, int n_bond, double avg_bond_length,
                     molecule_statistics_t &molecule_statistics,
                     bool generate_2D_coordinates, double * const confidence,
                     const std::map<std::string, std::string> &superatom, int n_letters,
                     std::string * const confidence_parameters, bool verbose,
		     const std::vector <bracket_t>&  brackets)
{
  std::string str;
  double scale = CC_BOND_LENGTH / avg_bond_length;
  // The indexes of "superatoms" and the bonds, that connect these superatoms to the molecule:
  std::vector<int> super_atoms, super_bonds;
  int anum;

  mol.SetDimension(2);
  mol.BeginModify();
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && i < MAX_ATOMS - 1 && bond[i].a < MAX_ATOMS - 1 && bond[i].b < MAX_ATOMS - 1)
      {
        atom_t* bond_atoms[] = { &atom[bond[i].a], &atom[bond[i].b] };

        for (int j = 0; j < 2; j++)
          if (bond_atoms[j]->n == 0)
            {
              if (create_atom(mol, *bond_atoms[j], scale, superatom, verbose))
                {
                  super_atoms.push_back(bond_atoms[j]->n);
                  // The current bond (next to be added) connects the super atom (bond_atoms[j]) with the molecule:
                  super_bonds.push_back(mol.NumBonds());
                }
            }

        if (bond[i].hash && !bond[i].wedge)
          {
            if (verbose)
              std::cout << "Creating hash bond #" << mol.NumBonds() << " " << atom[bond[i].a].n << "<->"
                        << atom[bond[i].b].n << ", order: " << bond[i].type << ", flags: " << OB_HASH_BOND << '.' << std::endl;

            if (atom[bond[i].a].anum == OXYGEN_ATOMIC_NUM || atom[bond[i].a].anum == HYDROGEN_ATOMIC_NUM || atom[bond[i].a].anum == FLUORINE_ATOMIC_NUM
                || atom[bond[i].a].anum == IODINE_ATOMIC_NUM || atom[bond[i].a].anum == CHLORINE_ATOMIC_NUM || atom[bond[i].a].anum == BROMINE_ATOMIC_NUM
                || atom[bond[i].a].anum == ARGON_ATOMIC_NUM || atom[bond[i].a].terminal)
              mol.AddBond(atom[bond[i].b].n, atom[bond[i].a].n, bond[i].type, OB_HASH_BOND);
            else
              mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_HASH_BOND);
          }
	else if (!bond[i].hash && bond[i].wedge)
          {
            if (atom[bond[i].a].anum == OXYGEN_ATOMIC_NUM || atom[bond[i].a].anum == HYDROGEN_ATOMIC_NUM || atom[bond[i].a].anum == FLUORINE_ATOMIC_NUM
                || atom[bond[i].a].anum == IODINE_ATOMIC_NUM || atom[bond[i].a].anum == CHLORINE_ATOMIC_NUM || atom[bond[i].a].anum == BROMINE_ATOMIC_NUM
                || atom[bond[i].a].anum == ARGON_ATOMIC_NUM || atom[bond[i].a].terminal)
              mol.AddBond(atom[bond[i].b].n, atom[bond[i].a].n, bond[i].type, OB_WEDGE_BOND);
            else
              mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_WEDGE_BOND);
          }
	else if (bond[i].arom)
          {
            if (verbose)
              std::cout << "Creating aromatic bond #" << mol.NumBonds() << " " << atom[bond[i].a].n << "->"
                        << atom[bond[i].b].n << ", order: " << AROMATIC_BOND_ORDER << '.' << std::endl;

            mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, AROMATIC_BOND_ORDER);
          }
        else
          {
            int bond_flags = 0;

            if (bond[i].up)
              bond_flags = OB_TORUP_BOND;
            else if (bond[i].down)
              bond_flags = OB_TORDOWN_BOND;

            if (verbose)
              std::cout << "Creating bond #" << mol.NumBonds() << " " << atom[bond[i].a].n << "->"
                        << atom[bond[i].b].n << ", type: " << bond[i].type << ", flags: " << bond_flags << '.' << std::endl;

	    if (bond[i].wedge && bond[i].hash) // wavy bonds
	      {
		mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type,  OB_WEDGE_OR_HASH_BOND);
		OBMolAtomIter a,b;
		FOR_ATOMS_OF_MOL(ai, mol)
		  if (ai->GetIdx() == atom[bond[i].b].n) b = ai;
		  else if (ai->GetIdx() == atom[bond[i].a].n) a = ai;
		SetTetrahedtalUnknown(a);
		SetTetrahedtalUnknown(b);
	      }
	    else
	      mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, bond_flags);
          }
      }
  mol.EndModify();

  mol.FindRingAtomsAndBonds();

  // Clear the counters of created OBAtom objects:
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      {
        atom[bond[i].a].n = 0;
        atom[bond[i].b].n = 0;
      }

  // The logic below calculates the information both for molecule statistics and for confidence function:

  OBBondIterator bond_iter;
  int Num_HashBonds = 0;
  int Num_WedgeBonds = 0;
  int Num_DoubleBonds = 0;

  // This block modifies the molecule:
  for (OBBond *b = mol.BeginBond(bond_iter); b; b = mol.NextBond(bond_iter))
    {
      if (b->IsInRing())
        {
          // Clear any indication of "/" and "\" double bond stereochemistry:
          b->UnsetUp();
          b->UnsetDown();
        }
      else
        // Clear all aromaticity information for the bond (affects "num_aromatic" variable below):
        b->UnsetAromatic();
      if (b->IsHash())
	Num_HashBonds++;
      if (b->IsWedge())
	Num_WedgeBonds++;
      if (b->IsDouble() && (b->GetBeginAtom()->IsOxygen() ||  b->GetEndAtom()->IsOxygen()))
	Num_DoubleBonds++;
    }

  std::vector<int> Num_Rings(8, 0); // number of rings of the given size (e.g. "Num_Rings[2]" = number of rings of size 2)
  int num_rings = 0; // total number of rings
  int num_aromatic = 0; // total number of aromatic rings

  {
    OBMol mol_without_R(mol);
    mol_without_R.BeginModify();
    OBAtomIterator atom_iter_R;
    for (OBAtom *a = mol_without_R.BeginAtom(atom_iter_R); a; a = mol_without_R.NextAtom(atom_iter_R))
      if (a->GetAtomicNum() == 0)
	a->SetAtomicNum(CARBON_ATOMIC_NUM);
    mol_without_R.EndModify();
    mol_without_R.FindRingAtomsAndBonds();
    // Get the Smallest Set of Smallest Rings:
    std::vector<OBRing*> vr = mol_without_R.GetSSSR();

    for (std::vector<OBRing*>::iterator iter = vr.begin(); iter != vr.end(); iter++)
      {
	num_rings++;
	if ((*iter)->IsAromatic())
	  num_aromatic++;
	if ((*iter)->Size() < 8)
	  Num_Rings[(*iter)->Size()]++;
      }
  }

  // Get a list of contiguous fragments sorted by size from largest to smallest:
  std::vector<std::vector<int> > cfl;
  mol.ContigFragList(cfl);

  molecule_statistics.rotors = mol.NumRotors();
  molecule_statistics.fragments = cfl.size();
  molecule_statistics.rings456 = Num_Rings[5] + Num_Rings[6] + Num_Rings[4];
  molecule_statistics.rings56 = Num_Rings[5] + Num_Rings[6];
  molecule_statistics.num_atoms = mol.NumAtoms();
  molecule_statistics.num_bonds = mol.NumBonds();

  if (confidence)
    {
      int C_Count = 0;
      int N_Count = 0;
      int O_Count = 0;
      int F_Count = 0;
      int S_Count = 0;
      int Cl_Count = 0;
      int Br_Count = 0;
      int R_Count = 0;
      int Xx_Count = 0;
      int Si_Count = 0;
      int Me_Count = 0;
      int U_Count = 0;
      int Li_Count = 0;
      int PositiveCharge = 0;

      OBAtomIterator atom_iter;

      for (OBAtom *a = mol.BeginAtom(atom_iter); a; a = mol.NextAtom(atom_iter))
        {
          if (a->IsCarbon())
	    {
	      C_Count++;
	      AliasData *ad = (AliasData *) a->GetData("UserLabel");
              if (ad != NULL && ad->GetAlias() != " " && !ad->GetAlias().empty())
		{
		  C_Count--;
		  Me_Count++;
		}
	    }
          else if (a->IsNitrogen())
            N_Count++;
          else if (a->IsOxygen())
            O_Count++;
          else if (a->IsSulfur())
            S_Count++;
          else if (a->GetAtomicNum() == FLUORINE_ATOMIC_NUM)
            F_Count++;
          else if (a->GetAtomicNum() == CHLORINE_ATOMIC_NUM)
            Cl_Count++;
          else if (a->GetAtomicNum() == BROMINE_ATOMIC_NUM)
            Br_Count++;
	  else if (a->GetAtomicNum() == SILICONE_ATOMIC_NUM)
            Si_Count++;
	  else if (a->GetAtomicNum() == URANIUM_ATOMIC_NUM)
	    U_Count++;
	  else if (a->GetAtomicNum() == LITHIUM_ATOMIC_NUM)
	    Li_Count++;
          else if (a->GetAtomicNum() == 0)
            {
              AliasData *ad = (AliasData *) a->GetData("UserLabel");
              if (ad != NULL && ad->GetAlias() != "Xx" && ad->GetAlias() != "*")
                R_Count++;
              else
                Xx_Count++;
            }
	  else
	    {
              AliasData *ad = (AliasData *) a->GetData("UserLabel");
              if (ad != NULL && ad->GetAlias() != "Xx" && ad->GetAlias() != " "  && ad->GetAlias() != "*" && !ad->GetAlias().empty())
                R_Count++;
	    }
	  if (a->GetFormalCharge() > 0) PositiveCharge++;
        }

      molecule_statistics.num_organic_non_carbon_atoms = N_Count + F_Count + Cl_Count + Br_Count;
      int num_angles = 0;
      //vector<double> angles;
      if (super_atoms.empty())
	{
	  OBBondIterator bond_iter1;
	  for (OBBond *b1 = mol.BeginBond(bond_iter1); b1; b1 = mol.NextBond(bond_iter1))
	    {
	      unsigned int a = b1->GetBeginAtomIdx();
	      unsigned int b = b1->GetEndAtomIdx();
	      if ( b1->GetBeginAtom()->IsHydrogen() || b1->GetEndAtom()->IsHydrogen() || b1->GetLength() == 0)
		continue;
	      OBBondIterator bond_iter2;
	      for (OBBond *b2 = mol.BeginBond(bond_iter2); b2; b2 = mol.NextBond(bond_iter2))
		{
		  if (b1->GetIdx() < b2->GetIdx())
		    {
		      unsigned int c = b2->GetBeginAtomIdx();
		      unsigned int d = b2->GetEndAtomIdx();
		      if ( b2->GetBeginAtom()->IsHydrogen() || b2->GetEndAtom()->IsHydrogen()  || b2->GetLength() == 0)
			continue;
		      double angle;
		      bool found = false;
		      if (a == c && b != d)
			{
			  angle = mol.GetAngle(b1->GetEndAtom(), b1->GetBeginAtom(), b2->GetEndAtom());
			  found = true;
			}
		      else if (a == d && b != c)
			{
			  angle = mol.GetAngle(b1->GetEndAtom(), b1->GetBeginAtom(), b2->GetBeginAtom());
			  found = true;
			}
		      else if (b == c && a != d)
			{
			  angle = mol.GetAngle(b1->GetBeginAtom(), b1->GetEndAtom(), b2->GetEndAtom());
			  found = true;
			}
		      else if (b == d && a != c)
			{
			  angle = mol.GetAngle(b1->GetBeginAtom(), b1->GetEndAtom(), b2->GetBeginAtom());
			  found = true;
			}
		      if (found )
			{
			  //angles.push_back(angle);
			  if (angle < 20)
			    num_angles++;
			}
		    }
		}
	    }
	}
      molecule_statistics.num_small_angles = num_angles;
      /*if (!angles.empty())
	sort(angles.begin(),angles.end());
      for (size_t kk = 0; kk < angles.size(); kk++)
	cout << angles[kk] << endl;
      cout << "===" << endl;
      */
      //cout<<C_Count<<" "<<N_Count<<" "<<O_Count<<" "<<F_Count<<" "<<S_Count<<" "<<Cl_Count<<" "<<Br_Count<<" "<<Si_Count<<" "<<U_Count<<" "<<Me_Count<<" "<<R_Count<<" "<<
      //Xx_Count<<" "<<num_rings<<" "<<num_aromatic<<" "<<molecule_statistics.fragments<<" "<<Num_Rings<<" "<<Num_HashBonds<<" "<<Num_WedgeBonds<<" "<<Num_DoubleBonds<<" "<<PositiveCharge<<endl;
      int x[] = {C_Count,N_Count,O_Count,F_Count,S_Count,Cl_Count,Br_Count,Si_Count,U_Count,Me_Count,Li_Count,R_Count, Xx_Count,num_rings,num_aromatic,molecule_statistics.fragments,Num_Rings[5],Num_Rings[6],Num_HashBonds,Num_WedgeBonds,Num_DoubleBonds,PositiveCharge,n_letters,mol.NumAtoms(),mol.NumHvyAtoms(),mol.NumBonds()};
      int n = sizeof(x)/sizeof(int);
      *confidence = confidence_function(x,n);

      if (confidence_parameters)
	{
          std::stringstream cpss;
	  for (int  i = 0; i < n-1; i++)
	    cpss << x[i]<<",";
	  cpss << x[n-1];
	  //       0           1             2              3              4            5               6             7              8             9             10             11
	  //	  cpss<<C_Count<<","<<N_Count<<","<<O_Count<<","<<F_Count<<","<<S_Count<<","<<Cl_Count<<","<< Br_Count<<","<<Si_Count<<","<<U_Count<<","<<Me_Count<<","<<Li_Count<<","<<R_Count<<","<<
	    //  12           13             14                   15                                16                 17                  18                   19                    20
	  //   Xx_Count<<","<<num_rings<<","<<num_aromatic<<","<<molecule_statistics.fragments<<","<<Num_Rings[5]<<","<<Num_Rings[6]<<","<<Num_HashBonds<<","<<Num_WedgeBonds<<","<<Num_DoubleBonds
	    //            21                  22            23                   24                       25
	  //   <<","<<PositiveCharge<<","<<n_letters<<","<<mol.NumAtoms()<<","<<mol.NumHvyAtoms()<<","<<mol.NumBonds();
	  *confidence_parameters = cpss.str();
	}
    }

  if (generate_2D_coordinates)
    {
      if (verbose && !super_atoms.empty())
        std::cout << "Generating 2D coordinates for atoms " << super_atoms << " and bonds " << super_bonds << std::endl;

      for (unsigned int i = 0; i < super_atoms.size(); i++)
        {
          groupRedraw(&mol, super_bonds[i], super_atoms[i], true);
        }

      if (!brackets.empty())
	{
	  int group_id = 1;
	  //const size_t BUFF_SIZE = 32768;
	  char buff[BUFF_SIZE];
	  OBSetData *group = new OBSetData;
	  group->SetAttribute("_SGroupBrackets");
	  group->SetOrigin(userInput);
	  for (size_t i = 0; i < brackets.size(); i++)
	    {
	      double x1 = double(brackets[i].box.x1) * scale;
	      double y1 = -double(brackets[i].box.y1) * scale;
	      double x2 = double(brackets[i].box.x2) * scale;
	      double y2 = -double(brackets[i].box.y2) * scale;

              std::set<OBAtom*> visited;
              std::set<OBAtom*> nbrs;
	      OBAtom* start(NULL);
	      OBAtom* finish(NULL);
	      int bond1 = -1;
	      int bond2 = -1;
	      OBBondIterator bond_iter;
	      for (OBBond *b = mol.BeginBond(bond_iter); b; b = mol.NextBond(bond_iter))
		{
		  OBAtom *a1 = b->GetBeginAtom();
		  OBAtom *a2 = b->GetEndAtom();
		  if (a1->y() > y2 && a1->y() < y1 && a2->y() > y2 && a2->y() < y1)
		    {
		      if (a1->x() < x1 && a2->x() > x1)
			{
			  start = a1;
			  visited.insert(a1);
			  nbrs.insert(a2);
			  bond1 = b->GetIdx();
			}
		      if (a1->x() > x1 && a2->x() < x1)
			{
			  start = a2;
			  visited.insert(a2);
			  nbrs.insert(a1);
			  bond1 = b->GetIdx();
			}
		      if (a1->x() < x2 && a2->x() > x2)
			{
			  finish = a2;
			  visited.insert(a2);
			  nbrs.insert(a1);
			  bond2 = b->GetIdx();
			}
		      if (a1->x() > x2 && a2->x() < x2)
			{
			  finish = a1;
			  visited.insert(a1);
			  nbrs.insert(a2);
			  bond2 = b->GetIdx();
			}
		    }
		}
	      if (bond1 < 0 || bond2 < 0 || bond1 == bond2)
		continue;
	      OBBond *b1 = mol.GetBond(bond1);
	      OBBond *b2 = mol.GetBond(bond2);

	      OBPairData *i1 = new OBPairData;
	      i1->SetAttribute("_SGroup");
	      i1->SetOrigin(userInput);
	      snprintf(buff, BUFF_SIZE, "%d", group_id);
	      i1->SetValue(buff);
	      b1->SetData(i1);

	      OBPairData *i2 = new OBPairData;
	      i2->SetAttribute("_SGroup");
	      i2->SetOrigin(userInput);
	      snprintf(buff, BUFF_SIZE, "%d", group_id);
	      i2->SetValue(buff);
	      b2->SetData(i2);

	      while (!nbrs.empty())
		{
		  OBAtom *a = *nbrs.begin();
		  OBBondIterator bond_iter2;
		  for (OBAtom *b = a->BeginNbrAtom(bond_iter2); b; b = a->NextNbrAtom(bond_iter2))
		    {
		      if (visited.find(b) == visited.end())
			nbrs.insert(b);
		    }
		  visited.insert(a);
		  nbrs.erase(a);
		}
	      visited.erase(start);
	      visited.erase(finish);
	      for (std::set<OBAtom*>::iterator it = visited.begin(); it != visited.end(); ++it)
		{
		  OBAtom *a = *it;
		  OBPairData *i3 = new OBPairData;
		  i3->SetAttribute("_SGroup");
		  i3->SetOrigin(userInput);
		  snprintf(buff, BUFF_SIZE, "%d", group_id);
		  i3->SetValue(buff);
		  a->SetData(i3);
		}

	      OBSetData *coords = new OBSetData;
	      coords->SetAttribute("Coordinates");
	      coords->SetOrigin(userInput);

	      OBPairData *dx1 = new OBPairData;
	      dx1->SetAttribute("x1");
	      dx1->SetOrigin(userInput);
	      snprintf(buff, BUFF_SIZE, "%10.4f", x1);
	      dx1->SetValue(buff);
	      coords->AddData(dx1);

	      OBPairData *dy1 = new OBPairData;
	      dy1->SetAttribute("y1");
	      dy1->SetOrigin(userInput);
	      snprintf(buff, BUFF_SIZE, "%10.4f", y1);
	      dy1->SetValue(buff);
	      coords->AddData(dy1);

	      OBPairData *dx2 = new OBPairData;
	      dx2->SetAttribute("x2");
	      dx2->SetOrigin(userInput);
	      snprintf(buff, BUFF_SIZE, "%10.4f", x2);
	      dx2->SetValue(buff);
	      coords->AddData(dx2);

	      OBPairData *dy2 = new OBPairData;
	      dy2->SetAttribute("y2");
	      dy2->SetOrigin(userInput);
	      snprintf(buff, BUFF_SIZE, "%10.4f", y2);
	      dy2->SetValue(buff);
	      coords->AddData(dy2);

	      OBPairData *lbl = new OBPairData;
	      lbl->SetAttribute("label");
	      lbl->SetOrigin(userInput);
	      lbl->SetValue(brackets[i].a);
	      coords->AddData(lbl);

	      group->AddData(coords);
	      group_id++;
	    }
	  mol.SetData(group);
	}
    }
}

molecule_statistics_t calculate_molecule_statistics(
    std::vector<atom_t> &atom, const std::vector<bond_t> &bond, int n_bond, double avg_bond_length,
    const std::map<std::string, std::string> &superatom, bool verbose)
{
  molecule_statistics_t molecule_statistics;

  #pragma omp critical
  {
    OBMol mol;
    std::vector <bracket_t> brackets;
    create_molecule(mol, atom, bond, n_bond, avg_bond_length, molecule_statistics, false, NULL, superatom, 0, NULL, false, brackets);
    mol.Clear();
  }

  if (verbose)
    std::cout << "Molecule fragments: " << molecule_statistics.fragments << '.' << std::endl;

  return molecule_statistics;
}

const std::string get_formatted_structure(
    std::vector<atom_t> &atom, const std::vector<bond_t> &bond, int n_bond,
    const std::string &format, const std::string &embedded_format,
    molecule_statistics_t &molecule_statistics,
    double &confidence, bool show_confidence,
    double avg_bond_length, double scaled_avg_bond_length, bool show_avg_bond_length,
    const int * const resolution, const int * const page, const box_t * const surrounding_box,
    const std::map<std::string, std::string> &superatom, int n_letters, bool show_learning,
    int resolution_iteration, bool verbose, const std::vector<bracket_t>& brackets)
{
  std::ostringstream strstr;
  #pragma omp critical
  {
    OBMol mol;
    std::string confidence_parameters;
    create_molecule(mol, atom, bond, n_bond, avg_bond_length, molecule_statistics, format == "sdf" || format == "mol" || format == "sd" || format == "mdl", &confidence, superatom,
		    n_letters, &confidence_parameters, verbose, brackets);

    // Add hydrogens to the entire molecule to fill out implicit valence spots:
    mol.AddHydrogens(true, false); // polarOnly, correctForPh
    // Find all chiral atom centers:
    mol.FindChiralCenters();

    // Clear any indication of 2D "wedge" notation:
    OBBondIterator bond_iter;

    for (OBBond *b = mol.BeginBond(bond_iter); b; b = mol.NextBond(bond_iter))
      {
        if (!b->GetBeginAtom()->IsChiral() && !b->GetEndAtom()->IsChiral())
          {
            //b->UnsetHash();
            b->UnsetWedge();
          }
	if (!b->GetBeginAtom()->IsChiral() && (b->IsWedge() || b->IsHash()))
	  {
	    OBAtom* ba = b->GetBeginAtom();
	    OBAtom* ea = b->GetEndAtom();
	    b->SetBegin(ea);
	    b->SetEnd(ba);
	  }
      }

    // Add single bonds based on atom proximity:
    mol.ConnectTheDots();
    // Copies each disconnected fragment as a separate OBMol:
    //mol.Separate();
    // Deletes all atoms except for the largest contiguous fragment:
    mol.StripSalts(MIN_A_COUNT);

    if (show_confidence)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Confidence_estimate");
        std::ostringstream cs;
        cs << confidence;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (show_learning)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Confidence_parameters");
        std::ostringstream cs;
        cs << confidence_parameters;
        label->SetValue(cs.str());
        mol.SetData(label);

	OBPairData *label1 = new OBPairData;
        label1->SetAttribute("Resolution_iteration");
        std::ostringstream cs1;
        cs1 << resolution_iteration;
        label1->SetValue(cs1.str());
        mol.SetData(label1);
      }

    if (show_avg_bond_length)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Average_bond_length");
        std::ostringstream cs;
        cs << scaled_avg_bond_length;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (resolution)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Resolution");
        std::ostringstream cs;
        cs << *resolution;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (page)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Page");
        std::ostringstream cs;
        cs << *page;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (surrounding_box)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Surrounding_box");
        std::ostringstream cs;
        cs << surrounding_box->x1 << 'x' << surrounding_box->y1 << '-' << surrounding_box->x2 << 'x' << surrounding_box->y2;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (!embedded_format.empty())
      {
        std::string value;

        OBConversion conv;

        conv.SetOutFormat(embedded_format.c_str());

        value = conv.WriteString(&mol);
        trim(value);

        if (!value.empty())
          {
            OBPairData *label = new OBPairData;
            label->SetAttribute(embedded_format);
            label->SetValue(value.c_str());
            mol.SetData(label);
            if (embedded_format == "inchi")
              {
                conv.SetOptions("K", conv.OUTOPTIONS);

                value = conv.WriteString(&mol);
                trim(value);

                label = new OBPairData;
                label->SetAttribute("InChI_key");
                label->SetValue(value.c_str());
                mol.SetData(label);
              }
          }
      }

    OBConversion conv;

    conv.SetOutFormat(format.c_str());
    conv.Read(&mol);

    strstr << conv.WriteString(&mol, true);

    if (format == "smi" || format == "can")
      {
        if (show_avg_bond_length)
          strstr << " " << scaled_avg_bond_length;
        if (resolution)
          strstr << " " << *resolution;
        if (show_confidence)
          strstr << " " << confidence;
	if (show_learning)
          strstr << " " << confidence_parameters<<" "<<resolution_iteration;
        if (page)
          strstr << " " << *page;
        if (surrounding_box)
          strstr << " "<< surrounding_box->x1 << 'x' << surrounding_box->y1 << '-' << surrounding_box->x2 << 'x' << surrounding_box->y2;
      }

    strstr << std::endl;

    mol.Clear();
  }

  if (verbose)
    std::cout << "Structure length: " << strstr.str().length() << ", molecule fragments: " << molecule_statistics.fragments << '.' << std::endl;

  return (strstr.str());
}

