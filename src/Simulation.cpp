#include "spear/Simulation.hpp"

#include "OpenMM.h"

#include "spear/Forcefield.hpp"
#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"

using namespace Spear;

void Simulation::initialize_plugins() {
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory()
    );
}

Simulation::Simulation() : system_(new OpenMM::System) {
}

Simulation::~Simulation() {
}

bool Simulation::add_molecule(const Molecule& mol, const BondedForcefield& ff) {

    auto masses = ff.masses(mol);

    // See if its already been added, and bail if it has
    for (auto& old_mol : molecules_) {
        if (&old_mol.get() == &mol) {
            return false;
        }
    }

    for (auto mass : masses) {
        system_->addParticle(mass);
    }

    ff.add_forces(mol, *system_);

    // Prevent adding molecules multiple times
    molecules_.push_back(mol);

    return true;
}

void Simulation::add_non_bonded_force(const NonBondedForcefield& ff) {
    ff.add_forces(molecules_, *system_);
}

void Simulation::set_periodic_vectors(const chemfiles::UnitCell& cell) {
    const auto& matrix = cell.matrix();
    periodic1_ = {matrix[0][0], 0, 0};
    periodic2_ = {matrix[1][0], matrix[1][1], 0};
    periodic3_ = {matrix[2][0], matrix[2][1], matrix[2][2]};

    periodic1_ *= OpenMM::NmPerAngstrom;
    periodic2_ *= OpenMM::NmPerAngstrom;
    periodic3_ *= OpenMM::NmPerAngstrom;

    uses_periodic_ = true;
}

void Simulation::constrain_particle(size_t idx) {
    system_->setParticleMass(idx, 0.0);
}

void Simulation::add_langevin(double temperature, double friction, double stepsize) {
    integrator_ = std::make_unique<OpenMM::LangevinIntegrator>(temperature, friction, stepsize);
}

void Simulation::initialize_context(const std::string& platform) {

    // Do nothing if everything is already initialized
    if (context_) {
        return;
    }

    if (!integrator_) {
        integrator_ = std::make_unique<OpenMM::VerletIntegrator>(0.002);
    }

    if (uses_periodic_) {
        system_->setDefaultPeriodicBoxVectors(
            {periodic1_[0], periodic1_[1], periodic1_[2]},
            {periodic2_[0], periodic2_[1], periodic2_[2]},
            {periodic3_[0], periodic3_[1], periodic3_[2]}
        );
    }

    context_ = std::make_unique<OpenMM::Context>(*system_, *integrator_,
        OpenMM::Platform::getPlatformByName(platform)
    );

    size_t position_count = 0;
    for (auto& mol : molecules_) {
        position_count += mol.get().size();
    }

    std::vector<OpenMM::Vec3> positions;
    positions.reserve(position_count);
    for (auto& mol : molecules_) {
        auto& mol_positions = mol.get().positions();
        for (auto& pos : mol_positions) {
            positions.emplace_back(pos[0] * OpenMM::NmPerAngstrom,
                                   pos[1] * OpenMM::NmPerAngstrom,
                                   pos[2] * OpenMM::NmPerAngstrom
            );
        }
    }

    context_->setPositions(positions);
}

void Simulation::randomize_velocities(double temperature) {
    initialize_context();
    context_->setVelocitiesToTemperature(temperature);
}

void Simulation::minimize(double tolerance, std::size_t max_iterations) {
    initialize_context();

    OpenMM::LocalEnergyMinimizer::minimize(*context_, tolerance, static_cast<int>(max_iterations));
}

void Simulation::dynamic_steps(std::size_t steps) {
    initialize_context();

    integrator_->step(static_cast<int>(steps));
}

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#elif defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4018 )
#endif

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

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning( pop )
#endif
