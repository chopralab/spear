#include "spearmint.h"
#include "spearmint_statics.hpp"

#include "chemfiles/Trajectory.hpp"
#include "chemfiles/Frame.hpp"

#include "spear/Molecule.hpp"
#include "spear/Grid.hpp"

std::map<std::string, char> residue_codes({
    {"GLU", 'E'}, {"ASP", 'D'}, {"ALA", 'A'}, {"VAL", 'V'}, {"PHE", 'F'}, {"LEU", 'L'},
    {"TRP", 'W'}, {"MET", 'M'}, {"ILE", 'I'}, {"PRO", 'P'}, {"GLY", 'G'}, {"SER", 'S'},
    {"TYR", 'Y'}, {"CYS", 'C'}, {"ASN", 'N'}, {"THR", 'T'}, {"GLN", 'Q'}, {"ARG", 'R'},
    {"LYS", 'K'}, {"HIS", 'H'}
});

uint64_t spear_receptor_atom_count() {
    CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
    return receptor->size();
}

uint64_t spear_receptor_atoms(float* pos) {
    CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
    return get_atom_positions(*receptor, pos);
}

#define GET_RESIDUE(res, atom, receptor)					\
	if (atom >= receptor->size()) {							\
		set_error("Atom index out of bonds.");				\
		return 0;											\
	}														\
	auto res = receptor->topology().residue_for_atom(atom); \
	if (!res) {												\
		set_error("Atom not in a residue.");				\
		return 0;											\
	}

uint8_t spear_receptor_atom_chain(uint64_t atom) {
	CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
	GET_RESIDUE(res, atom, receptor);

	return res->get<chemfiles::Property::STRING>("chainid").value_or(" ")[0];
}

uint64_t spear_receptor_atom_resi(uint64_t atom) {
	CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
	GET_RESIDUE(res, atom, receptor);

	return static_cast<bool>(res->id())? (*res->id()): 0;
}

uint8_t spear_receptor_atom_resn(uint64_t atom) {
	CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
	GET_RESIDUE(res, atom, receptor);

	auto resn_iter = residue_codes.find(res->name());
	return resn_iter == residue_codes.end() ? ' ' : resn_iter->second;
}

uint64_t spear_receptor_atom_element(uint64_t atom) {
	CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
	if (atom >= receptor->size()) {
		set_error("Atom index out of bonds.");
		return 0;
	}
	return static_cast<size_t>((*receptor)[atom].atomic_number());
}

uint64_t spear_receptor_atom_details(uint8_t* cids, uint64_t* resi,
									 uint8_t* resn, uint64_t* elements) {
    CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
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

uint64_t spear_receptor_bond_count() {
    CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
    return receptor->topology().bonds().size();
}

uint64_t spear_receptor_bonds(uint64_t* bonds) {
    CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
    return get_bonds(*receptor, bonds);
}

uint64_t spear_receptor_bonds_in(uint64_t* atoms, uint64_t atoms_size, uint64_t* atoms_out) {
	CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
	std::set<size_t> atoms_to_check(atoms, atoms + atoms_size);

	auto res = receptor->get_bonds_in(atoms_to_check);

	for (size_t i = 0; i < res.size(); ++i) {
		atoms_out[3 * i + 0] = res[i].source();
		atoms_out[3 * i + 1] = res[i].target();
		atoms_out[3 * i + 2] = res[i].order();
	}

	return res.size();
}

uint64_t spear_receptor_set_positions(const float* positions) {
    CHECK_MOLECULE(receptor, "You must run initialize_receptor first");
    auto result = set_positions(*receptor, positions);
    gridrec = std::make_unique<Spear::Grid>(receptor->positions());
    return result;
}
