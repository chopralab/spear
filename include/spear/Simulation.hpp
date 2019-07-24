// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_SIMULATION_HPP
#define SPEAR_SIMULATION_HPP

#include "spear/exports.hpp"
#include <memory>
#include <vector>

#include "chemfiles.hpp"

#include "Eigen/Geometry"

// Forward declare OpenMM classes to avoid header dependancies
namespace OpenMM {
class System;
class Context;
class Platform;
class Integrator;
}

namespace Spear {

class NonBondedForcefield;
class BondedForcefield;
class Molecule;

class SPEAR_EXPORT Simulation final {

public:

    static void initialize_plugins();

    Simulation();

    ~Simulation();

    /// Adds a molecule to the system with a given forcefield.
    /// The molecule must stay valid during all operations performed by this
    /// class.
    bool add_molecule(const Molecule& mol, const BondedForcefield& ff);

    void add_non_bonded_force(const NonBondedForcefield& ff);

    void set_periodic_vectors(const chemfiles::UnitCell& cell);

    /// Constrains the positions of a coordinate by setting its mass to zero
    void constrain_particle(size_t idx);

    /// Adds a Langevin integrator to the context instead of verlet
    void add_langevin(double temperature, double friction, double stepsize);

    /// Initializes the context of the simulation and prepares it for further
    /// operations.
    void initialize_context(const std::string& platform = "Reference");

    /// Adds random velocities to the system
    void randomize_velocities(double temperature = 300.0);

    /// Run minimization
    void minimize(double tolerance, std::size_t max_iterations);

    /// Run dynamics
    void dynamic_steps(std::size_t steps);

    /// Get the system's time in pico seconds
    double time();

    /// Get the potential energy of the system
    double potential_energy();

    /// Get the kinetic energy of the system
    double kinetic_energy();

    /// Get the positions of the system
    std::vector<Eigen::Vector3d> positions();

    /// Get the velocities of the system
    std::vector<Eigen::Vector3d> velocities();

    /// Get the forces of the system
    std::vector<Eigen::Vector3d> forces();

private:
    /// We want to store pointers as we want to compare if it is the EXACT SAME
    /// molecule being added twice.
    std::vector<std::reference_wrapper<const Molecule>> molecules_;

    std::unique_ptr<OpenMM::System> system_;
    std::unique_ptr<OpenMM::Context> context_;
    std::unique_ptr<OpenMM::Platform> platform_;
    std::unique_ptr<OpenMM::Integrator> integrator_;

    bool uses_periodic_ = false;
    Eigen::Vector3d periodic1_, periodic2_, periodic3_;
};

}

#endif
