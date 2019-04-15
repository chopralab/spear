// License:

#ifndef SPEAR_SIMULATION_HPP
#define SPEAR_SIMULATION_HPP

#include "spear/exports.hpp"
#include <memory>
#include <vector>

#include "Eigen/Geometry"

// Forward declare OpenMM classes to avoid header dependancies
namespace OpenMM {
class System;
class Context;
class Platform;
class Integrator;
}

namespace Spear {

class Forcefield;
class Molecule;

class SPEAR_EXPORT Simulation {

public:

    Simulation();

    ~Simulation();

    /// Adds a molecule to the system with a given forcefield.
    /// The molecule must stay valid during all operations performed by this
    /// class.
    virtual bool add_molecule(const Molecule& mol, const Forcefield& ff);

    /// Adds a Langevin integrator
    virtual void add_langevin(double temperature, double friction, double stepsize);

    /// Initializes the context of the simulation and prepares it for further
    /// operations.
    virtual void initialize_context();

    /// Run minimization
    virtual void minimize(double tolerance, std::size_t max_iterations);

    /// Run dynamics
    virtual void dynamic_steps(std::size_t steps);

    /// Get the system's time in pico seconds
    virtual double time();

    /// Get the potential energy of the system
    virtual double potential_energy();

    /// Get the kinetic energy of the system
    virtual double kinetic_energy();

    /// Get the positions of the system
    virtual std::vector<Eigen::Vector3d> positions();

    /// Get the velocities of the system
    virtual std::vector<Eigen::Vector3d> velocities();

    /// Get the forces of the system
    virtual std::vector<Eigen::Vector3d> forces();

private:
    /// We want to store pointers as we want to compare if it is the EXACT SAME
    /// molecule being added twice.
    std::vector<const Molecule*> molecules_;

    std::unique_ptr<OpenMM::System> system_;
    std::unique_ptr<OpenMM::Context> context_;
    std::unique_ptr<OpenMM::Platform> platform_;
    std::unique_ptr<OpenMM::Integrator> integrator_;
};

}

#endif
