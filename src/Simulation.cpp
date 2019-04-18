#include "spear/Simulation.hpp"

#include "OpenMM.h"

#include "spear/Forcefield.hpp"
#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"

using namespace Spear;

Simulation::Simulation() : system_(new OpenMM::System) {
}

Simulation::~Simulation() {
    context_ = nullptr;
    integrator_ = nullptr;
    system_ = nullptr;
}

bool Simulation::add_molecule(const Molecule& mol, const Forcefield& ff) {

    auto masses = ff.masses(mol);

    // See if its already been added, and bail if it has
    auto old_mol = std::find(molecules_.begin(), molecules_.end(), &mol);
    if (old_mol != molecules_.end()) {
        return false;
    }

    for (auto mass : masses) {
        system_->addParticle(mass);
    }

    ff.add_forces(mol, *system_);

    // Prevent adding molecules multiple times
    molecules_.push_back(&mol);

    return true;
}

void Simulation::add_langevin(double temperature, double friction, double stepsize) {
    integrator_ = std::make_unique<OpenMM::LangevinIntegrator>(temperature, friction, stepsize);
}

void Simulation::initialize_context() {

    // Do nothing if everything is already initialized
    if (context_) {
        return;
    }

    if (!integrator_) {
        integrator_ = std::make_unique<OpenMM::VerletIntegrator>(0.002);
    }

    context_ = std::make_unique<OpenMM::Context>(*system_, *integrator_);

    size_t position_count = 0;
    for (auto& mol : molecules_) {
        position_count += mol->size();
    }

    std::vector<OpenMM::Vec3> positions;
    positions.reserve(position_count);
    for (auto& mol : molecules_) {
        auto& mol_positions = mol->positions();
        for (auto& pos : mol_positions) {
            positions.emplace_back(pos[0] * OpenMM::NmPerAngstrom,
                                   pos[1] * OpenMM::NmPerAngstrom,
                                   pos[2] * OpenMM::NmPerAngstrom
            );
        }
    }

    context_->setPositions(positions);
}

void Simulation::minimize(double tolerance, std::size_t max_iterations) {
    initialize_context();

    OpenMM::LocalEnergyMinimizer::minimize(*context_, tolerance, max_iterations);
}

void Simulation::dynamic_steps(std::size_t steps) {
    initialize_context();

    integrator_->step(steps);
}

double Simulation::time() {
    initialize_context();

    auto state = context_->getState(0);
    return state.getTime();
}

double Simulation::potential_energy() {
    initialize_context();

    auto state = context_->getState(OpenMM::State::Energy);
    return state.getPotentialEnergy();
}

double Simulation::kinetic_energy() {
    initialize_context();

    auto state = context_->getState(OpenMM::State::Energy);
    return state.getKineticEnergy();
}

static std::vector<Eigen::Vector3d> convert_to_eigen(const std::vector<OpenMM::Vec3>& vec) {
    std::vector<Eigen::Vector3d> ret;
    ret.reserve(vec.size());
    for (auto current : vec) {
        ret.emplace_back(current[0] * OpenMM::AngstromsPerNm,
                         current[1] * OpenMM::AngstromsPerNm,
                         current[2] * OpenMM::AngstromsPerNm
        );
    }

    return ret;
}

std::vector<Eigen::Vector3d> Simulation::positions() {
    initialize_context();

    auto state = context_->getState(OpenMM::State::Positions);
    return convert_to_eigen(state.getPositions());
}

std::vector<Eigen::Vector3d> Simulation::velocities() {
    initialize_context();

    auto state = context_->getState(OpenMM::State::Velocities);
    return convert_to_eigen(state.getVelocities());
}

std::vector<Eigen::Vector3d> Simulation::forces() {
    initialize_context();

    auto state = context_->getState(OpenMM::State::Forces);
    return convert_to_eigen(state.getForces());
}
