#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"

#include "chemfiles.hpp"
#include <utility>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("Molecule") {
    SECTION("Molecule::iterator") {
        auto traj = chemfiles::Trajectory("data/palmitic.sdf");
        auto mol = Spear::Molecule(traj.read());

        auto begin = mol.begin();
        auto end = mol.end();

        CHECK(std::distance(begin, end) == 18);
        CHECK(std::distance(end, begin) == -18);
        CHECK(begin != end);

        // Bidirectional iterator
        auto first = *(begin++);
        CHECK(first == 0);
        CHECK(first.name() == "C");
        CHECK(first.degree() == 1);
        CHECK(first.atomic_number() == 6);

        auto third = *(++begin);
        CHECK(third == 2);
        CHECK(third.name() == "C");
        CHECK(third.degree() == 2);
        CHECK(third.atomic_number() == 6);

        auto last_atom_vertex = *(--end);
        CHECK(last_atom_vertex == 17);
        CHECK(last_atom_vertex.name() == "O");
        CHECK(last_atom_vertex.degree() == 1);
        CHECK(last_atom_vertex.atomic_number() == 8);
        --end;
        --end;
        auto carbonyl = *(end--);
        CHECK(carbonyl == 15);
        CHECK(carbonyl.name() == "C");
        CHECK(carbonyl.degree() == 3);
        CHECK(carbonyl.atomic_number() == 6);

        // Random access iterator
        CHECK(begin[2] == 4);

        CHECK(end > begin);
        CHECK(begin < end);
        CHECK(end - begin == 12);
        CHECK(begin - end == -12);
        CHECK(end - 12 == begin);

        std::ptrdiff_t diff = -12;
        CHECK(begin - diff == end);

        std::swap(begin, end);
        CHECK(end + 12 == begin);
        begin += diff;
        CHECK(end >= begin);
        CHECK(begin <= end);
    }

    SECTION("Get bonds in a set of atoms") {
        auto traj = chemfiles::Trajectory("data/tibolone.sdf");
        auto mol = Spear::Molecule(traj.read());

        CHECK(mol.bond(0,1).order() == Spear::Bond::SINGLE);
        CHECK_THROWS(mol.bond(0, 5));

        std::set<size_t> test = {0, 1, 2, 3, 4, 5, 20};
        auto bond = mol.get_bonds_in(test);

        CHECK(bond.size() == 5);
        CHECK(bond[0].source() == 0);
        CHECK(bond[0].target() == 1);

        CHECK(bond[1].source() == 1);
        CHECK(bond[1].target() == 2);
    }
}

TEST_CASE("Ring Finding") {
    SECTION("No rings") {
        auto traj = chemfiles::Trajectory("data/palmitic.sdf");
        auto mol = Spear::Molecule(traj.read());
        auto ring = mol.rings();
        auto sssr = mol.smallest_set_of_smallest_rings();

        CHECK(mol.size() == 18);
        CHECK(ring.size() == 0);
        CHECK(sssr.size() == 0);
    }

    SECTION("Pazopanib's rings") {
        auto traj = chemfiles::Trajectory("data/pazopanib.sdf");
        auto mol = Spear::Molecule(traj.read());
        auto ring = mol.rings();
        auto sssr = mol.smallest_set_of_smallest_rings();

        // Holds the total number of expected rings:
        //    ring_size -> #expected rings
        std::map<size_t,size_t> ring_counts = {{5,1}, {6,3}, {9,1}};

        // Decrease the count when we find a ring of given size
        for (const auto& i : ring) {
            ring_counts[i.size()] -= 1;
        }

        CHECK(ring_counts[5] == 0);
        CHECK(ring_counts[6] == 0);
        CHECK(ring_counts[9] == 0);

        //--------------------------------------------------------------
        // SSSR for pazopanib
        //--------------------------------------------------------------
        auto ring_iter = sssr.begin();
        CHECK(sssr.size() == 4);
        CHECK((ring_iter++)->size() == 5);
        CHECK((ring_iter++)->size() == 6);
        CHECK((ring_iter++)->size() == 6);
    }

    SECTION("Steroid rings") {
        auto traj = chemfiles::Trajectory("data/tibolone.sdf");
        auto mol = Spear::Molecule(traj.read());
        auto ring = mol.rings();
        auto sssr = mol.smallest_set_of_smallest_rings();

        std::map<size_t,size_t> ring_counts = {{5,1}, {6,3}, {9,1}, {10,2},
         {13,1}, {14,1}, {17,1}};

        for (const auto& i : ring) {
            ring_counts[i.size()] -= 1;
        }

        CHECK(ring_counts[5] == 0);
        CHECK(ring_counts[6] == 0);
        CHECK(ring_counts[9] == 0);
        CHECK(ring_counts[10] == 0);
        CHECK(ring_counts[13] == 0);
        CHECK(ring_counts[14] == 0);
        CHECK(ring_counts[17] == 0);

        //--------------------------------------------------------------
        // SSSR for tibolone
        //--------------------------------------------------------------
        CHECK(sssr.size() == 4);
        auto ring_iter = sssr.begin();
        CHECK((ring_iter++)->size() == 5);
        CHECK((ring_iter++)->size() == 6);
        CHECK((ring_iter++)->size() == 6);
        CHECK((ring_iter++)->size() == 6);
    }

    SECTION("Large rings - 1") {
        auto traj = chemfiles::Trajectory("data/large_rings.sdf");
        auto mol = Spear::Molecule(traj.read());
        auto sssr = mol.smallest_set_of_smallest_rings();

        CHECK(sssr.size() == 6);
        auto ring_iter = sssr.begin();
        CHECK((ring_iter++)->size() == 5);
        CHECK((ring_iter++)->size() == 6);
        CHECK((ring_iter++)->size() == 6);
        CHECK((ring_iter++)->size() == 6);
        CHECK((ring_iter++)->size() == 6);
        CHECK((ring_iter++)->size() == 6);
    }

    SECTION("Large rings - 2") {
        auto traj = chemfiles::Trajectory("data/large_rings.sdf");
        auto mol = Spear::Molecule(traj.read_step(1));
        auto sssr = mol.smallest_set_of_smallest_rings();

        CHECK(sssr.size() == 9);
        auto ring_iter = sssr.begin();
        CHECK((ring_iter++)->size() == 3);
        CHECK((ring_iter++)->size() == 3);
        CHECK((ring_iter++)->size() == 3);
        CHECK((ring_iter++)->size() == 3);
        CHECK((ring_iter++)->size() == 4);
        CHECK((ring_iter++)->size() == 4);
        CHECK((ring_iter++)->size() == 4);
        CHECK((ring_iter++)->size() == 4);
        CHECK((ring_iter++)->size() == 4);
    }

    SECTION("Large rings - 3") {
        auto traj = chemfiles::Trajectory("data/large_rings.sdf");
        auto mol = Spear::Molecule(traj.read_step(2));
        auto sssr = mol.smallest_set_of_smallest_rings();

        CHECK(sssr.size() == 9);
        auto ring_iter = sssr.begin();
        CHECK((ring_iter++)->size() == 3);
        CHECK((ring_iter++)->size() == 3);
        CHECK((ring_iter++)->size() == 3);
        CHECK((ring_iter++)->size() == 3);
        CHECK((ring_iter++)->size() == 4);
        CHECK((ring_iter++)->size() == 4);
        CHECK((ring_iter++)->size() == 4);
        CHECK((ring_iter++)->size() == 4);
        CHECK((ring_iter++)->size() == 4);
    }
}

TEST_CASE("Dimensionality") {
    SECTION("Point") {
        chemfiles::Frame f;
        f.add_atom(chemfiles::Atom("C"), {1, 2, 3});
        auto mol = Spear::Molecule(f);
        CHECK(mol.dimensionality() == 0);
    }

    SECTION("Linear") {
        chemfiles::Frame f;
        f.add_atom(chemfiles::Atom("C"), {1, 2, 3});
        f.add_atom(chemfiles::Atom("C"), {2, 3, 4});

        auto mol1 = Spear::Molecule(f);
        CHECK(mol1.dimensionality() == 1);

        f.add_atom(chemfiles::Atom("C"), {3, 4, 5});
        auto mol2 = Spear::Molecule(f);
        CHECK(mol2.dimensionality() == 1);

        f.add_atom(chemfiles::Atom("C"), {5, 6, 7});
        auto mol3 = Spear::Molecule(f);
        CHECK(mol3.dimensionality() == 1);

        f.add_atom(chemfiles::Atom("C"), {6, 8, 8});
        auto mol4 = Spear::Molecule(f);
        CHECK(mol4.dimensionality() == 2); // No longer linear

        f.add_atom(chemfiles::Atom("C"), {6, 7, 8});
        auto mol5 = Spear::Molecule(f);
        CHECK(mol5.dimensionality() == 2); // Still planar

        f.add_atom(chemfiles::Atom("C"), {7, 8, 8});
        auto mol6 = Spear::Molecule(f);
        CHECK(mol6.dimensionality() == 3); // Still planar
    }

    SECTION("Planar") {
        chemfiles::Frame f1;
        f1.add_atom(chemfiles::Atom("C"), {0, 0, 0});
        f1.add_atom(chemfiles::Atom("C"), {1, 0, 0});
        f1.add_atom(chemfiles::Atom("C"), {0, 1, 0});

        // L-shaped
        auto mol1 = Spear::Molecule(f1);
        CHECK(mol1.dimensionality() == 2);

        // T-shaped
        f1.add_atom(chemfiles::Atom("C"), {-1, 0, 0});
        auto mol2 = Spear::Molecule(f1);
        CHECK(mol2.dimensionality() == 2);

        // No longer planer
        f1.add_atom(chemfiles::Atom("C"), {0, 0, 1});
        auto mol3 = Spear::Molecule(f1);
        CHECK(mol3.dimensionality() == 3);        
    }
}

TEST_CASE("Default AtomType") {
    auto traj2 = chemfiles::Trajectory("data/tibolone.sdf");
    auto mol2 = Spear::Molecule(traj2.read());

    auto default_atom_type = mol2.get_default_atomtype();
    CHECK(default_atom_type->hybridization(0) == Spear::Hybridization::SP3);
    CHECK(default_atom_type->hybridization(3) == Spear::Hybridization::SP2);
    CHECK(default_atom_type->hybridization(4) == Spear::Hybridization::SP2);
    CHECK(default_atom_type->hybridization(7) == Spear::Hybridization::SP2);
    CHECK(default_atom_type->hybridization(8) == Spear::Hybridization::SP2);
    CHECK(default_atom_type->hybridization(20) == Spear::Hybridization::SP);
    CHECK(default_atom_type->hybridization(21) == Spear::Hybridization::SP);
}
