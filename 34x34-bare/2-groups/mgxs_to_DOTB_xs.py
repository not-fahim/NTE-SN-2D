################################################################
## This script reads multi-group cross sections from an OpenMC##
################################################################

import h5py
import numpy as np

cross_section_domain='material'  # could either cell or material

HDF5_FILE_PATH = 'mgxs-2g.h5'
OUTPUT_FILE_PATH = '34x34-bare-2g.xs'


#    Based on the h5dump, have cells or materials "1", "2", and "3".
material_map = {

    '1': 'fuel',
    '2': 'gad',
    '3': "cool"
}


def write_material_block(h5_file, output_file, domain_id, material_name):
    """
    Reads all cross sections for a given cell, formats them in scientific
    notation, and writes them to the output file.
    """
    
    print(f"Processing material '{material_name}' from cell '{domain_id}'...")

    # Material header
    output_file.write(f"mat {material_name}\n")
    
    # Total cross section (sigma_t)
    # try:
    #     total_xs = h5_file[f'/{cross_section_domain}/{domain_id}/total/average'][:]
    #     total_str = ' '.join([f'{x:.8e}' for x in total_xs])
    #     output_file.write(f"tot {total_str}\n")
    # except KeyError:
    #     print(f"  Warning: 'total' cross section not found for cell {domain_id}.")

    # Nu-fission cross section (nu_sigma_f)
    try:
        nuf_xs = h5_file[f'/{cross_section_domain}/{domain_id}/nu-fission/average'][:]
        nuf_str = ' '.join([f'{x:.8e}' for x in nuf_xs])
        output_file.write(f"nuf {nuf_str}\n")
    except KeyError:
        print(f"  Warning: 'nu-fission' cross section not found for cell {domain_id}.")
        
    # Fission spectrum (chi)
    try:
        chi = h5_file[f'/{cross_section_domain}/{domain_id}/chi/average'][:]
        chi_str = ' '.join([f'{x:.8e}' for x in chi])
        output_file.write(f"chi {chi_str}\n")
    except KeyError:
        print(f"  Warning: 'chi' xs not found for cell {domain_id}.")
    
    # absorption cross section (sigma_a)
    try:
        abs = h5_file[f'/{cross_section_domain}/{domain_id}/absorption/average'][:]
        abs_str = ' '.join([f'{x:.8e}' for x in abs])
        output_file.write(f"abs {abs_str}\n")
    except KeyError:
        print(f"  Warning: 'absorption' xs not found for cell {domain_id}.")
    
    # transport cross section (sigma_tr)
    try:
        tra = h5_file[f'/{cross_section_domain}/{domain_id}/transport/average'][:]
        tra_str = ' '.join([f'{x:.8e}' for x in tra])
        output_file.write(f"tot {tra_str}\n")
    except KeyError:
        print(f"  Warning: 'transport' xs not found for cell {domain_id}.")

    # --- Read 2D Dataset (scatter matrix) ---
    try:
        # Note the tab space in 'scatter matrix'
        scatter_matrix = h5_file[f'/{cross_section_domain}/{domain_id}/scatter matrix/average'][:]
        output_file.write("% rows are to groups and columns are from groups\n")
        output_file.write("sca\n")
        # Transpose and iterate through columns (which become rows)
        for row in scatter_matrix.T:
            # Format each number in the row with scientific notation (.8e)
            row_str = '\t'.join([f'{x:.8e}' for x in row])
            output_file.write(f"{row_str}\n")
    except KeyError:
         print(f"  Warning: 'scatter matrix' not found for cell {domain_id}.")

    # Add a blank line for separation
    output_file.write("\n")


def main():
    """Main function to read HDF5 and generate the materials.in file."""
    try:
        with h5py.File(HDF5_FILE_PATH, 'r') as h5f:
            with open(OUTPUT_FILE_PATH, 'w') as out_f:
                
                # Write header comments
                out_f.write("% material input file for DOTB\n")
                out_f.write(f"% Generated from {HDF5_FILE_PATH}\n\n")

                for domain_id, name in material_map.items():
                    write_material_block(h5f, out_f, domain_id, name)
        
        print(f"\nSuccessfully created '{OUTPUT_FILE_PATH}' with {len(material_map)} materials.")

    except FileNotFoundError:
        print(f"Error: The HDF5 file '{HDF5_FILE_PATH}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    main()