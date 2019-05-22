using System;
using System.Runtime.InteropServices;
using System.Text;

public static class SpearMintInterface
{

    // ----------------------- ERROR HANDLING -----------------------------

    [DllImport("spearmint")]
    public static extern string spear_get_error();

    // ---------------------- COMPLEX FUNCTIONS ---------------------------

    [DllImport("spearmint")]
    public static extern ulong spear_initialize_complex(string filename);

    [DllImport("spearmint")]
    public static extern ulong spear_write_complex(string filename);


    // --------------------- RECEPTOR FUNCTIONS ---------------------------
    [DllImport("spearmint")]
    public static extern ulong spear_initialize_receptor(string filename);

    [DllImport("spearmint")]
    public static extern ulong spear_receptor_atom_count();

    [DllImport("spearmint")]
    public static extern ulong spear_receptor_atoms([Out] float[] pos);

    [DllImport("spearmint")]
    public static extern ulong spear_receptor_atom_details([Out] char[] cids,
                                                           [Out] ulong[] resi,
                                                           [Out] char[] resn,
                                                           [Out] ulong[] elements);

    [DllImport("spearmint")]
    public static extern ulong spear_receptor_set_positions([In] float[] positions);

    [DllImport("spearmint")]
    public static extern ulong spear_receptor_bond_count();

    [DllImport("spearmint")]
    public static extern ulong spear_receptor_bonds([In] ulong[] bonds);

    // --------------------- LIGAND FUNCTIONS -----------------------------
    [DllImport("spearmint")]
    public static extern ulong spear_initialize_ligand(string filename);

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_atom_count();

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_atoms([Out] float[] pos);

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_atom_details([Out] ulong[] elements);

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_set_positions([In] float[] positions);

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_bond_count();

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_bonds([In] ulong[] bonds);

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_neighbors(ulong atom_idx, [Out] ulong[] neighbors);

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_is_adjacent(ulong atom1, ulong atom2);

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_add_bond(ulong atom1, ulong atom2);

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_remove_bond(ulong atom1, ulong atom2);

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_remove_hydrogens();

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_add_atom([In] ulong element,
            [In] float x, [In] float y, [In] float z
    );

    [DllImport("spearmint")]
    public static extern ulong spear_ligand_add_atom_to( [In] ulong atom, [In] ulong element,
            [In,Out] ref float x, [In,Out] ref float y, [In,Out] ref float z
    );

    // ---------------------- OTHER FUNCTIONS -----------------------------
    [DllImport("spearmint")]
    public static extern ulong spear_initialize_scoring(string obj_dir);

    [DllImport("spearmint")]
    public static extern float spear_calculate_score();

    public static void Main()
    {
        try
        {
            /// Read the file '3qox_pocket.pdb' from disk. Takes a single string argument.
            if (spear_initialize_receptor("3qox_pocket.pdb") == 0)
            {
                Console.WriteLine(spear_get_error());
                return;
            }

            /// Read the file '3qox_ligand.pdb' from disk. Supported formats are PDB and MOL2
            if (spear_initialize_ligand("3qox_ligand.sdf") == 0)
            {
                Console.WriteLine(spear_get_error());
            }

            /// Initialize the scoring function with the objective function given.
            if (spear_initialize_scoring("share/") == 0)
            {
                Console.WriteLine(spear_get_error());
            }

            // calculate_score caclulates the score of the ligand - receptor complex. Takes no arguments
            Console.WriteLine("Score is : " + spear_calculate_score());

            // Initialize atom data here.
            // ligand_atom_count takes no arguments and returns the number of atoms in the ligand
            ulong lig_atom_count = spear_ligand_atom_count();
            float[] positions = new float[lig_atom_count * 3];

            // spear_ligand_atoms fills a float[] with positionss about the ligand atoms. The array must be
            // at least 3 times the number of atoms in the ligand.
            spear_ligand_atoms(positions);

            for (ulong i = 0; i < lig_atom_count; ++i)
            {
                // Change the X position of atom 'i' by 5 angstroms
                positions[i * 3 + 0] = positions[i * 3 + 0] + 5.0f;

                // Change the Y position of atom 'i' by 4 angstroms
                positions[i * 3 + 1] = positions[i * 3 + 1] + 4.0f;

                // Change the Z position of atom 'i' by 3 angstroms
                positions[i * 3 + 2] = positions[i * 3 + 2] + 3.0f;
            }

            ulong[] elements = new ulong[lig_atom_count];
            spear_ligand_atom_details(elements);
            Console.Write("Elements are ");
            for (ulong i = 0; i < lig_atom_count; ++i)
            {
                Console.Write(" " + elements[i]);
            }
            Console.WriteLine();

            // Change the position of the ligand atoms to the postions stored in the float[]
            if (spear_ligand_set_positions(positions) == 0)
            {
                Console.WriteLine(spear_get_error());
            }

            // Nothing else needs to be done to recalculate the score.
            Console.WriteLine("Score is : " + spear_calculate_score());

            // The same as above, but with the receptor coordinates
            ulong rec_atom_count = spear_receptor_atom_count();
            float[] positions_rec = new float[rec_atom_count * 3];
            spear_receptor_atoms(positions_rec);

            for (ulong j = 0; j < rec_atom_count; ++j)
            {
                positions_rec[j * 3 + 0] = positions_rec[j * 3 + 0] + 5.0f;
                positions_rec[j * 3 + 1] = positions_rec[j * 3 + 1] + 4.0f;
                positions_rec[j * 3 + 2] = positions_rec[j * 3 + 2] + 3.0f;
            }

            spear_receptor_set_positions(positions_rec);

            // Should now be the as when we ran this the first time!
            Console.WriteLine("Score is : " + spear_calculate_score());

            // Obtain the number of bonds and create an array of three times its length
            ulong lig_bond_count = spear_ligand_bond_count();
            ulong[] bond = new ulong[lig_bond_count * 3];

            // Fills the array with bond information.
            // The first index gives the first atom index of the bond
            // The next gives the index for the second atom
            // The third gives bond information
            //     1 - single bond
            //     2 - double bond
            //     4 - triple bond
            //   255 - Aromatic  bond
            if (spear_ligand_bonds(bond) == 0)
            {
                Console.WriteLine(spear_get_error());
            }

            // Print out all ligand bonds!
            for (ulong k = 0; k < lig_bond_count; ++k)
            {
                Console.WriteLine("BOND: " + bond[k * 3 + 0] +
                                  " - " + bond[k * 3 + 1] +
                                  " : " + bond[k * 3 + 2]);
            }

            ulong[] neighs = new ulong[5];
            if (spear_ligand_neighbors(1, neighs) == 0)
            {
                Console.WriteLine(spear_get_error());
            }

            Console.WriteLine("Atom 1 neighbors: " +
                                    neighs[0] + " " +
                                    neighs[1] + " " +
                                    neighs[2] + " " +
                                    neighs[3]
            );

            ulong receptor_atom_count = spear_receptor_atom_count();
            char[] cids = new char[receptor_atom_count];
            ulong[] resi = new ulong[receptor_atom_count];
            char[] resn = new char[receptor_atom_count];
            ulong[] relements = new ulong[receptor_atom_count];
            spear_receptor_atom_details(cids, resi, resn, relements);

            ulong current_res = 0;
            Console.Write("Residues are ");
            for (ulong i = 0; i < receptor_atom_count; ++i)
            {
                if (resi[i] != current_res)
                {
                    Console.WriteLine();
                    Console.Write(" " + cids[i] + " " + resn[i] + " " + resi[i] + " ");
                    current_res = resi[i];
                }
                Console.Write(relements[i] + " ");
            }
            Console.WriteLine();

            /**************************************************************
             * Modification time!!!!!!!!!!
             **************************************************************/

            // Attempt (and fail to) add an atom to the ligand because the atom
            // is saturated with hydrogens (see ligand file)!
            float x = 0.0f, y = 0.0f, z = 0.0f;
            if (spear_ligand_add_atom_to(0, 6, ref x, ref y, ref z) == 0)
            {
                Console.WriteLine("Expected error: (" + spear_get_error() + ")");
            }

            // Now we can do it after removing hydrogens! (do it twice because its a nitrogen)
            // Note the decrease in score, this is due to collisions formed with these new atoms
            spear_ligand_remove_hydrogens();
            ulong new_atom1 = spear_ligand_add_atom_to(0, 6, ref x, ref y, ref z);
            Console.WriteLine("Added point at {0} {1} {2}", x, y, z);

            Console.WriteLine("Score is : " + spear_calculate_score());

            ulong new_atom2 = spear_ligand_add_atom_to(0, 6, ref x, ref y, ref z);
            Console.WriteLine("Added point at {0} {1} {2}", x, y, z);

            if (new_atom1 == 0 || new_atom2 == 0)
            {
                Console.WriteLine(spear_get_error() + "\n");
            }

            Console.WriteLine("Score is : " + spear_calculate_score());

            if (spear_ligand_is_adjacent(0, new_atom1) != 0)
            {
                Console.WriteLine("Good, they are connected!");
            }

            // Add a very strained bond!
            if (spear_ligand_add_bond(new_atom1, new_atom2) == 0)
            {
                Console.WriteLine(spear_get_error());
            }
            spear_write_complex("strained.pdb");

            // Try to add the same bond twice, but we can't because it already exists
            if (spear_ligand_add_bond(new_atom1, new_atom2) == 0)
            {
                Console.WriteLine("Expected error: (" + spear_get_error() + ")");
            }

            // Try to add a bond to the nitrogen, but we can't because its saturated
            if (spear_ligand_add_bond(0, 8) == 0)
            {
                Console.WriteLine("Expected error: (" + spear_get_error() + ")");
            }

            // Let's remove the original bond and add ours instead
            if (spear_ligand_remove_bond(0, 1) == 0)
            {
                Console.WriteLine(spear_get_error());
            }

            if (spear_ligand_add_bond(0, 8) == 0)
            {
                Console.WriteLine(spear_get_error());
            }

            if (spear_ligand_add_atom(6, 47.77969f, 54.91803f, 3.669073f) == 0)
            {
                Console.WriteLine(spear_get_error());
            }

            Console.WriteLine("Score after addition of random point is : " + spear_calculate_score());

            // save everything to disk
            spear_write_complex("out.pdb");

            spear_initialize_complex("1aaq.pdb");
            Console.WriteLine("Score after loading non-ligand protein: " + spear_calculate_score());

            receptor_atom_count = spear_receptor_atom_count();
            cids = new char[receptor_atom_count];
            resi = new ulong[receptor_atom_count];
            resn = new char[receptor_atom_count];
            relements = new ulong[receptor_atom_count];
            spear_receptor_atom_details(cids, resi, resn, relements);

            current_res = 0;
            Console.WriteLine("Residues are ");
            for (ulong i = 0; i < receptor_atom_count; ++i)
            {
                if (cids[i] == '\0') {
                    Console.WriteLine("ERROR");
                }
                if (resi[i] != current_res)
                {
                    Console.WriteLine();
                    Console.Write(" " + cids[i] + " " + resn[i] + " " + resi[i] + " ");
                    current_res = resi[i];
                }
                Console.Write(relements[i] + " ");
            }
            Console.WriteLine();

            if (spear_initialize_complex("5tz6.pdb") == 0)
            {
                Console.WriteLine("Expected error: (" + spear_get_error() + ")");
            }

            spear_initialize_complex("5zt6.pdb");
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }
    }
}
