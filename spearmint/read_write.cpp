#include <memory>

#include "spearmint.h"
#include "spearmint_statics.hpp"

#include "chemfiles/Trajectory.hpp"
#include "chemfiles/Frame.hpp"

#include "spear/Molecule.hpp"
#include "spear/Grid.hpp"
#include "spear/atomtypes/IDATM.hpp"

/******************************************************************************
 * Stolen from Lemon :)
 *****************************************************************************/

const std::set<std::string> common_residues{
    {"CSD"}, {"PCA"}, {"DLE"}, {"KCX"}, {"CAS"}, {"CSO"}, {"PTR"},
    {"CME"}, {"TPO"}, {"SEP"}, {"MLY"}, {"HYP"}, {"MSE"},
    {"CYS"}, {"TRP"}, {"MET"}, {"HIS"}, {"TYR"}, {"GLN"}, {"PHE"},
    {"ASN"}, {"PRO"}, {"ARG"}, {"THR"}, {"ASP"}, {"ILE"}, {"LYS"},
    {"SER"}, {"GLU"}, {"VAL"}, {"GLY"}, {"ALA"}, {"LEU"}, {"EDO"},
    {"FAD"}, {"FMN"}, {"NAD"}, {"NAP"}, {"CLA"}, {"HEM"}, {"HEA"}, {"HEB"},
    {"HEC"}, {"ADP"}, {"ATP"}, {"GDP"}, {"GTP"}, {"UNL"}, {"CIT"}, {"FLC"},
    {"BE7"}, {"MHA"}, {"DHD"}, {"B3P"}, {"BTB"}, {"NHE"}, {"GOL"}, {"DTP"},
    {"SAM"}, {"SIA"}, {"ICT"}, {"EPE"}, {"MES"},
    {"PG6"}, {"PE7"}, {"PG5"}, {"PEU"}, {"PGE"}, {"PIG"}, {"PE8"}, {"PE4"},
    {"P33"}, {"C8E"}, {"OTE"}, {"XPE"}, {"N8E"}, {"DR6"}, {"PEG"}, {"2PE"},
    {"P6G"}, {"1PE"}, {"SPM"}, {"SPK"}, {"SPD"}, {"1PG"}, {"PG4"}, {"MYR"},
    {"OLA"}, {"OLB"}, {"OLC"}, {"PLM"}, {"PEE"}, {"LHG"}, {"MC3"}, {"PAM"},
    {"PRO"}, {"HYP"}, {"PCA"}, {"A"}, {"T"}, {"G"}, {"C"},
};

template <typename Container>
inline void separate_residues(const chemfiles::Frame& input,
                              const Container& accepted_residues,
                              chemfiles::Frame& new_frame) {

    const auto& residues = input.topology().residues();
    const auto& positions = input.positions();
    const auto& old_bonds = input.topology().bonds();
    const auto& bond_ord = input.topology().bond_orders();

    std::unordered_map<size_t, size_t> old_to_new;
    std::unordered_set<size_t> accepted_atoms;
    for (auto res_id : accepted_residues) {
        const auto& res = residues[res_id];

        auto res_new = chemfiles::Residue(res.name(), *(res.id()));

        for (size_t res_atom : res) {
            new_frame.add_atom(input[res_atom], positions[res_atom]);
            res_new.add_atom(new_frame.size() - 1);
            old_to_new.insert({res_atom, new_frame.size() - 1});
            accepted_atoms.insert(res_atom);
        }

        res_new.set("chainid", res.get("chainid")->as_string());
        new_frame.add_residue(std::move(res_new));
    }

    for (size_t bond_idx = 0; bond_idx < old_bonds.size(); ++bond_idx) {
        if (accepted_atoms.count(old_bonds[bond_idx][0]) &&
            accepted_atoms.count(old_bonds[bond_idx][1])) {

            new_frame.add_bond(old_to_new[old_bonds[bond_idx][0]],
                               old_to_new[old_bonds[bond_idx][1]],
                               bond_ord[bond_idx]);
        }
    }
}

inline void protein_and_ligand(const chemfiles::Frame& input, size_t ligand_id,
                               chemfiles::Frame& protein,
                               chemfiles::Frame& ligand) {
    const auto& topo = input.topology();
    const auto& residues = topo.residues();

    std::list<size_t> accepted_residues;
    for (size_t res_id = 0; res_id < residues.size(); ++res_id) {
        if (res_id == ligand_id) {
            continue;
        }
        accepted_residues.push_back(res_id);
    }

    separate_residues(input, accepted_residues, protein);
    separate_residues(input, std::list<size_t>({ligand_id}), ligand);
}

size_t spear_initialize_complex(const char* filename) {
    bool recall_sf = false;
    if (receptor != nullptr) {
        recall_sf = true;
        receptor.reset();
    }

    if (ligand != nullptr) {
        recall_sf = true;
        ligand.reset();
    }

    try {

        chemfiles::Trajectory traj(filename, 'r');

        auto frame = traj.read();
        size_t current = 0;
        bool found = false;
        for (auto res : frame.topology().residues()) {
            if (res.size() < 10) {
                current++;
                continue;
            }

            auto is_common = common_residues.find(res.name());
            if (is_common != common_residues.end()) {
                current++;
                continue;
            }
            found = true;
            break;
        }

        if (!found) {
            set_error("Error in loading complex: everything is protein.");
            receptor = std::make_unique<Spear::Molecule>(frame);
            gridrec  = std::make_unique<Spear::Grid>(receptor->positions());
            ligand   = nullptr;
            return 0;
        } else {
            chemfiles::Frame protein, ligand_frame;
            protein_and_ligand(frame, current, protein, ligand_frame);
            receptor = std::make_unique<Spear::Molecule>(protein);
            gridrec = std::make_unique<Spear::Grid>(receptor->positions());
            ligand = std::make_unique<Spear::Molecule>(ligand_frame);
        }

        if (recall_sf && receptor) {
            auto atomtype_name = receptor->add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);
            receptor->atomtype(atomtype_name);
        }

        if (recall_sf && ligand) {
            auto atomtype_name = ligand->add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);
            ligand->atomtype(atomtype_name);
        }

        return 1;
    } catch (std::exception& e) {
        set_error(std::string("Error in loading complex: ") + e.what());
        return 0;
    }
}

size_t spear_write_complex(const char* filename) {
    if (ligand == nullptr || receptor == nullptr) {
        set_error(
            "You must run initialize_ligand and initialize_receptor before writing a complex"
        );
        return 0;
    }

    try {

        chemfiles::Frame complex;
        complex.reserve(ligand->size() + receptor->size());
        complex.resize(receptor->size());
        complex.set_topology(receptor->topology());

        for (size_t i = 0; i < complex.topology().residues().size(); ++i) {
            // ugh
            auto& res = complex.topology().residues()[i];
            const_cast<chemfiles::Residue&>(res).set("is_standard_pdb", true);
        }

        auto out_positions = complex.positions();
        const auto& rec_positions = receptor->positions();
        for (size_t i = 0; i < rec_positions.size(); ++i) {
            auto& atom_pos = rec_positions[i];
            out_positions[i] = {atom_pos[0], atom_pos[1], atom_pos[2]};
        }

        chemfiles::Residue ligres("LIG", 9999);
        for (auto liga : *ligand) {
            auto& atom_pos = liga.position();
            complex.add_atom(chemfiles::Atom(liga.name(),
                                             Spear::Element::Name[liga.atomic_number()] ),
                {atom_pos[0], atom_pos[1], atom_pos[2]}
            );
            ligres.add_atom(complex.size() - 1);
        }
        complex.add_residue(ligres);

        auto ligand_start = receptor->size();
        for (auto bond : ligand->topology().bonds()) {
            complex.add_bond(ligand_start + bond[0], ligand_start + bond[1]);
        }

        chemfiles::Trajectory traj(filename, 'w');
        traj.write(complex);

        return 1;
    } catch (std::exception& e) {
        set_error(std::string("Error in saving complex: ") + e.what());
        return 0;
    }
}

size_t spear_initialize_receptor(const char* filename) {
    if (receptor != nullptr) {
        receptor.reset();
    }

    try {
        chemfiles::Trajectory traj(filename, 'r');

        receptor = std::make_unique<Spear::Molecule>(traj.read());
        gridrec = std::make_unique<Spear::Grid>(receptor->positions());

        return 1;
    } catch (std::exception& e) {
        set_error(std::string("Error in loading receptor ") + e.what());
        return 0;
    }
}

size_t spear_initialize_ligand(const char* filename) {
    if (ligand != nullptr) {
        ligand.reset();
    }

    try {
        chemfiles::Trajectory traj(filename, 'r');
        ligand = std::make_unique<Spear::Molecule>(traj.read());
        return 1;
    } catch (std::exception& e) {
        set_error(std::string("Error in loading ligand ") + e.what());
        return 0;
    }

}
