// License:

#ifndef SPEAR_FORCEFIELD_HPP
#define SPEAR_FORCEFIELD_HPP

namespace OpenMM {
class System;
}

namespace Spear {

class SPEAR_EXPORT Forcefield {

public:

    /// Add force to system. `Simulation` take ownership of the force.
    /// This function must be called after all the particles have been added
    /// to the system.
    virtual void add_forces(const Molecule& mol, OpenMM::System& system) const = 0;

    /// Retrieves the masses of all atoms in a molecule
    virtual std::vector<double> masses(const Molecule& mol) const = 0;
};

}

#endif
