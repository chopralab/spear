#include "spear/Molecule.hpp"
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
        CHECK_THROWS(*end);

        // Bidirectional iterator
        auto first = *(begin++);
        CHECK(first == 0);
        CHECK(first.name() == "C");
        CHECK(first.neighbor_count() == 1);
        CHECK(first.atomic_number() == 6);

        auto third = *(++begin);
        CHECK(third == 2);
        CHECK(third.name() == "C");
        CHECK(third.neighbor_count() == 2);
        CHECK(third.atomic_number() == 6);

        auto last_atom_vertex = *(--end);
        CHECK(last_atom_vertex == 17);
        CHECK(last_atom_vertex.name() == "O");
        CHECK(last_atom_vertex.neighbor_count() == 1);
        CHECK(last_atom_vertex.atomic_number() == 8);
        --end;
        --end;
        auto carbonyl = *(end--);
        CHECK(carbonyl == 15);
        CHECK(carbonyl.name() == "C");
        CHECK(carbonyl.neighbor_count() == 3);
        CHECK(carbonyl.atomic_number() == 6);

        // Random access iterator
        CHECK_THROWS(begin[1000]);
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

    SECTION("Ring Finding") {
        auto traj0 = chemfiles::Trajectory("data/palmitic.sdf");
        auto mol0 = Spear::Molecule(traj0.read());
        auto ring0 = mol0.rings();

        CHECK(mol0.size() == 18);
        CHECK(ring0.size() == 0);

        auto traj = chemfiles::Trajectory("data/pazopanib.sdf");
        auto mol = Spear::Molecule(traj.read());
        auto rings = mol.rings();

        // Holds the total number of expected rings:
        //    ring_size -> #expected rings
        std::map<size_t,size_t> ring_counts = {{5,1}, {6,3}, {9,1}};

        // Decrease the count when we find a ring of given size
        for (const auto& i : rings) {
            ring_counts[i.size()] -= 1;
        }

        CHECK(ring_counts[5] == 0);
        CHECK(ring_counts[6] == 0);

        auto traj2 = chemfiles::Trajectory("data/tibolone.sdf");
        auto mol2 = Spear::Molecule(traj2.read());
        auto rings2 = mol2.rings();

        std::map<size_t,size_t> ring_counts2 = {{5,1}, {6,3}, {9,1}, {10,2},
         {13,1}, {14,1}, {17,1}};

        for (const auto& i : rings2) {
            ring_counts2[i.size()] -= 1;
        }

        CHECK(ring_counts2[5] == 0);
        CHECK(ring_counts2[6] == 0);
        CHECK(ring_counts2[9] == 0);
        CHECK(ring_counts2[10] == 0);
        CHECK(ring_counts2[13] == 0);
        CHECK(ring_counts2[14] == 0);
        CHECK(ring_counts2[17] == 0);
    }

    SECTION("Get bonds in a set of atoms") {
        auto traj = chemfiles::Trajectory("data/tibolone.sdf");
        auto mol = Spear::Molecule(traj.read());

        std::set<size_t> test = {0, 1, 2, 3, 4, 5, 20};
        auto bond = mol.get_bonds_in(test);

        CHECK(bond.size() == 5);
        CHECK(boost::source(bond[0], mol.graph()) == 0);
        CHECK(boost::target(bond[0], mol.graph()) == 1);

        CHECK(boost::source(bond[1], mol.graph()) == 1);
        CHECK(boost::target(bond[1], mol.graph()) == 2);
    }
}
