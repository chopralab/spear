#include "spearmint.h"
#include "spearmint_statics.hpp"

#include "chemfiles/Trajectory.hpp"
#include "chemfiles/Frame.hpp"

#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/Graph_impl.hpp"
#include "spear/Grid.hpp"

std::map<std::string, char> residue_codes({
    {"GLU", 'E'}, {"ASP", 'D'}, {"ALA", 'A'}, {"VAL", 'V'}, {"PHE", 'F'}, {"LEU", 'L'},
    {"TRP", 'W'}, {"MET", 'M'}, {"ILE", 'I'}, {"PRO", 'P'}, {"GLY", 'G'}, {"SER", 'S'},
    {"TYR", 'Y'}, {"CYS", 'C'}, {"ASN", 'N'}, {"THR", 'T'}, {"GLN", 'Q'}, {"ARG", 'R'},
    {"LYS", 'K'}, {"HIS", 'H'}
});

size_t spear_receptor_atom_count() {
    CHECK_MOLECULE(ligand, "You must run initialize_receptor first");
    return receptor->size();
}

size_t spear_receptor_atoms(float* pos) {
    CHECK_MOLECULE(ligand, "You must run initialize_receptor first");
    return get_atom_positions(*receptor, pos);
}

size_t spear_receptor_atom_details(char* cids, size_t* resi,
                                   char* resn, size_t* elements) {
    CHECK_MOLECULE(ligand, "You must run initialize_receptor first");
    try {

        for (auto residue : receptor->topology().residues()) {
            auto cid = residue.get<chemfiles::Property::STRING>("chainid").value_or(" ")[0];
            auto resn_iter = residue_codes.find(residue.name());
            auto res_name = resn_iter == residue_codes.end() ? ' ' : resn_iter->second;
            auto id = static_cast<bool>(residue.id()) ? (*residue.id()) : 0;
            for (auto atom : residue) {
                cids[atom] = cid;
                resi[atom] = id;
                resn[atom] = res_name;
                elements[atom] = static_cast<size_t>((*receptor)[atom].atomic_number());
            }
        }

        return 1;
    } catch (std::exception& e) {
        set_error(std::string("Error creating receptor atom arrays: ") + e.what());
        return 0;
    }
}

size_t spear_receptor_bond_count() {
    CHECK_MOLECULE(ligand, "You must run initialize_receptor first");
    return receptor->topology().bonds().size();
}

size_t spear_receptor_bonds(size_t* bonds) {
    CHECK_MOLECULE(ligand, "You must run initialize_receptor first");
    return get_bonds(*receptor, bonds);
}

size_t spear_receptor_set_positions(const float* positions) {
    CHECK_MOLECULE(ligand, "You must run initialize_receptor first");
    auto result = set_positions(*receptor, positions);
    gridrec = std::make_unique<Spear::Grid>(receptor->positions());
    return result;
}