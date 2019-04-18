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

#endif
