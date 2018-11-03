// License:

#ifndef SPEAR_IDATM_HPP
#define SPEAR_IDATM_HPP

#include "spear/AtomType.hpp"

namespace Spear {

class SPEAR_EXPORT IDATM : public AtomType {
public:

    /// algorithm based on E.C. Meng / R.A. Lewis paper
    /// "Determination of Molecular Topology and Atomic Hybridization
    /// States from Heavy Atom Coordinates", J. Comp. Chem., v12#7, 891-898
    /// and on example code from idatm.f implementation by E.C. Meng

    /// differences: No boron types.  Double-bonded Npls are split off
    ///   as N2.  Sox split into Sxd (sulfoxide), and Son (sulfone).
    ///   Carbons in aromatic rings are type Car.  Aromatic oxygens are Oar.
    /// still missing types: C1-,O1+,O1

    std::vector<size_t> type_atoms_3d(const Molecule& mol) override;

    std::vector<size_t> type_atoms_order(const Molecule& mol) override;

private:
    /// infallible pass:  type hydrogens / deuteriums and compute number of
    /// heavy atoms connected to each atom.
    /// also applies templated residues and marks them are mapped
    void infallible_(const Molecule& mol);

    /// valence pass: elements d valences > 1
    /// also fills the 'redo' array for hard to determine types
    ///
    /// valence 4
    ///  C must be sp3 (C3)
    ///  N must be part of an N-oxide (Nox) or a quaternary
    ///    amine (N3+)
    ///  P must be part of a phosphate (Pac), a P-oxide (Pox)
    ///    or a quaternary phosphine (P3+)
    ///  S must be part of a sulfate, sulfonate or sulfamate
    ///    (Sac), or sulfone (Son)
    ///
    /// valence 3
    /// calculate the three bond angles and average them;
    /// since hydrogens may be missing, cannot count on valence
    /// to determine the hybridization state.  Average bond angle
    /// assists in discriminating hybridization
    ///	C may be sp3 (C3), sp2 (C2), or part of a carboxylate
    ///		(Cac)
    ///	N may be sp3 (N3), sp2, or planar (as in amides and
    ///		aniline deriviatives), or part of a nitro
    ///		group (Ntr)
    ///	S may be, depending on oxidation state, sulfoxide (Sxd)
    ///		or S3+
    ///
    /// valence 2
    /// calculate the bond angle and assign a tentative atom
    /// type accordingly (a single angle is often not a good
    /// indicator of type).  Mark these atoms for further
    /// analysis by putting a non-zero value for them in the
    /// 'redo' array.
    ///	C may be sp3 (C3), sp2 (C2), or sp (C1)
    ///	N may be sp3 (N3), sp2 or planar (Npl), or sp (N1)
    ///	O and S are sp3 (O3 and S3, respectively)
    std::vector<size_t> valence_(const Molecule& mol);

    /// terminal pass: determine types of valence 1 atoms.  These were typed by
    /// element only in previous pass, but can be typed more accurately
    /// now that the atoms they are bonded to have been typed.
    /// Bond lengths are used in this pass to perform this assignment.
    void terminal_(const Molecule& mol, std::vector<size_t>& redo_);

    /// Re-examine all atoms with non-zero 'redo' values and
    ///   retype them if necessary
    void redo_(const Molecule& mol, const std::vector<size_t>& redo_);

    /// change isolated sp2 carbons to sp3 since it is
    /// impossible for an atom to be sp2 hybrizided if all its
    /// neighbors are sp3 hybridized.  In addition, a carbon atom cannot
    /// be double bonded to a carboxylate carbon, phosphate phosphorus,
    /// sulfate sulfur, sulfone sulfur, sulfoxide sulfur, or sp1 carbon.
    /// Addition not in original idatm: Nox also
    void fix_C2_(const Molecule& mol);

    /// 1) make decisions about the charge states of nitrogens.  If an
    ///    sp3 nitrogen is bonded to sp3 carbons and/or hydrogens and/or
    ///    deuteriums only, assume that it is positively charged (the pKa
    ///    of its conjugate acid is probably high enough that the
    ///    protonated form predominates at physiological pH).  If an sp2
    ///    carbon is bonded to three planar nitrogens, it may be part of
    ///    a guanidinium group.  Make the nitrogens positively charged
    ///    (Ng+) if guanidine or similar structures can be ruled out (if
    ///    'noplus' is false).
    /// 2) make carboxyl oxygens negatively charged even if the proton is
    ///    present (the pKa of the carboxyl group is probably low enough
    ///    that the unprotonated form predominates at physiological pH).
    void charges_(const Molecule& mol);

    /// Assign aromaticity
    ///
    /// 1) Check that all the atoms of the ring are planar types
    /// 2) Check bond lengths around the ring; see if they are
    ///    consistent with aromatic bond lengths
    void aromatic_(const Molecule& mol);

    /// Split off heavy-atom-valence-Npls that have no hydrogens as type N2.
    /// Discrimination criteria is the average bond length of the two
    /// heavy-atom bonds (shorter implies more double-bond character,
    /// thereby no hydrogen).
    void pass8_(const Molecule& mol);

    /// "pass 9": Assign aromatic types 
    void pass9_(const Molecule& mol);

    /// "pass 10": change O2- to O3- for sulfates, phosphates, N-oxide, S3- for
    /// thiophosphate and other terminal atoms now that we have more types.
    void pass10_(const Molecule& mol);

    /// "pass 11": Additional Special groups
    void pass11_(const Molecule& mol);

    std::vector<size_t> atom_types_;
    std::multimap<size_t, size_t> aromatic_ring_sizes_;
    
    /// number of heavy atoms bonded
    std::vector<size_t> heavys_;

    /// Have we typed the atom?
    std::vector<bool> mapped_;
};

template<> std::string atomtype_name<IDATM>();

template<> std::string atomtype_name_for_id<IDATM>(size_t id);

template<> size_t atomtype_id_for_name<IDATM>(std::string name);

}

#endif
