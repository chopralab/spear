#ifndef SPEARMINT_STATICS_HPP
#define SPEARMINT_STATICS_HPP

#include <memory>

namespace Spear {
    class Molecule;
    class ScoringFunction;
    class Grid;
}

extern std::unique_ptr<Spear::Molecule> receptor;
extern std::unique_ptr<Spear::Molecule> ligand;
extern std::unique_ptr<Spear::ScoringFunction> score;
extern std::unique_ptr<Spear::Grid> gridrec;

void set_error(const std::string& error);

size_t get_atom_positions(const Spear::Molecule& mol, float* pos);
size_t set_positions(Spear::Molecule& mol, const float* positions);
size_t get_bonds(Spear::Molecule& mol, size_t* bonds);

#define CHECK_MOLECULE(mol, mesg) if (mol == nullptr){ set_error(mesg); return 0;}

#endif
