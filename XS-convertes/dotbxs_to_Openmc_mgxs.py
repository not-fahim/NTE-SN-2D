import argparse
import textwrap

def parse_dotb_file(filepath):
    """
    Parses a DOTB-style multigroup cross section file.

    Args:
        filepath (str): The path to the input text file.

    Returns:
        dict: A dictionary where keys are material names and values are
              dictionaries of their cross-section data.
        int: The number of energy groups detected.
    """
    materials_data = {}
    current_material = None
    parsing_sca = False
    num_groups = 0

    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('%'):
                    continue

                parts = line.split()
                keyword = parts[0].lower()

                if keyword == 'mat':
                    current_material = parts[1]
                    materials_data[current_material] = {}
                    parsing_sca = False
                elif keyword in ['nuf', 'chi', 'abs', 'tot', 'fis']:
                    data = [float(x) for x in parts[1:]]
                    if not num_groups:
                        num_groups = len(data)
                    materials_data[current_material][keyword] = data
                    parsing_sca = False
                elif keyword == 'sca':
                    materials_data[current_material]['sca'] = []
                    parsing_sca = True
                elif parsing_sca:
                    row_data = [float(x) for x in parts]
                    materials_data[current_material]['sca'].append(row_data)
                    if len(materials_data[current_material]['sca']) == num_groups:
                        parsing_sca = False
    except UnicodeDecodeError:
        # Fallback for non-UTF8 files
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('%'):
                    continue
                # ... (rest of the logic is duplicated, a bit redundant but safe)
                parts = line.split()
                keyword = parts[0].lower()
                if keyword == 'mat':
                    current_material = parts[1]
                    materials_data[current_material] = {}
                    parsing_sca = False
                elif keyword in ['nuf', 'chi', 'abs', 'tot', 'fis']:
                    data = [float(x) for x in parts[1:]]
                    if not num_groups:
                        num_groups = len(data)
                    materials_data[current_material][keyword] = data
                    parsing_sca = False
                elif keyword == 'sca':
                    materials_data[current_material]['sca'] = []
                    parsing_sca = True
                elif parsing_sca:
                    row_data = [float(x) for x in parts]
                    materials_data[current_material]['sca'].append(row_data)
                    if len(materials_data[current_material]['sca']) == num_groups:
                        parsing_sca = False


    if not num_groups:
        raise ValueError("Could not determine the number of energy groups from the file.")

    return materials_data, num_groups


def generate_openmc_script(materials_data, num_groups, output_path):
    """
    Generates the Python script for OpenMC based on parsed data.

    Args:
        materials_data (dict): The parsed material data.
        num_groups (int): The number of energy groups.
        output_path (str): The path to write the generated script to.
    """
    script_content = textwrap.dedent(f"""
    import openmc
    import numpy as np

    print("--- Generating OpenMC MGXS Library and Materials ---")

    # Create a {num_groups}-group structure.
    # NOTE: The energy boundaries here are placeholders. You should
    # replace them with the actual energy boundaries for your problem.
    # The boundaries should go from lowest to highest energy.
    group_bounds = np.logspace(-5, 7, {num_groups + 1}) # Example: 0.01 eV to 10 MeV
    groups = openmc.mgxs.EnergyGroups(group_bounds)

    # Initialize the library
    mg_cross_sections_file = openmc.MGXSLibrary(groups)

    # --------------------------------------------------------------------------
    # Add Cross Section Data for Each Material
    # --------------------------------------------------------------------------
    """)

    for name, data in materials_data.items():
        var_name = name.replace('-', '_')
        script_content += textwrap.dedent(f"""
    # --- Material: {name} ---
    {var_name}_xsdata = openmc.XSdata('{name}', groups)
    {var_name}_xsdata.order = 0 # Isotropic scattering
    """)
        if 'tot' in data:
            script_content += f"{var_name}_xsdata.set_total({data['tot']}, temperature=294.)\n"
        if 'abs' in data:
            script_content += f"{var_name}_xsdata.set_absorption({data['abs']}, temperature=294.)\n"
        if 'fis' in data:
            script_content += f"{var_name}_xsdata.set_fission({data['fis']}, temperature=294.)\n"
        if 'nuf' in data:
            script_content += f"{var_name}_xsdata.set_nu_fission({data['nuf']}, temperature=294.)\n"
        if 'chi' in data:
            script_content += f"{var_name}_xsdata.set_chi({data['chi']}, temperature=294.)\n"
        if 'sca' in data:
            sca_matrix_transposed = [list(row) for row in zip(*data['sca'])]
            matrix_str_rows = [', '.join(map(str, row)) for row in sca_matrix_transposed]
            matrix_str = '[' + '], ['.join(matrix_str_rows) + ']'
            script_content += textwrap.dedent(f"""
    # The scattering matrix from the file is S(g_out, g_in).
    # OpenMC's `set_scatter_matrix` expects the matrix indexed by S[g_in][g_out].
    # Therefore, we use the transposed matrix from the original file.
    # The shape must be (G, G, M+1), so we add a third dimension for M=0.
    {var_name}_scatter_matrix = [[{matrix_str}]]
    {var_name}_scatter_matrix = np.array({var_name}_scatter_matrix)
    {var_name}_scatter_matrix = np.rollaxis({var_name}_scatter_matrix, 0, 3)
    {var_name}_xsdata.set_scatter_matrix({var_name}_scatter_matrix, temperature=294.)
    """)
        script_content += f"mg_cross_sections_file.add_xsdata({var_name}_xsdata)\n\n"

    script_content += textwrap.dedent("""
    # --------------------------------------------------------------------------
    # Export Library and Create Materials XML
    # --------------------------------------------------------------------------

    # Export the MGXS library to an HDF5 file
    mgxs_filepath = 'mgxs.h5'
    mg_cross_sections_file.export_to_hdf5(mgxs_filepath)
    print(f"--> MGXS library saved to {{mgxs_filepath}}")

    # Create a Materials XML file that links to the MGXS library
    materials_dict = {}
    """)
    material_names = [f"'{name}'" for name in materials_data.keys()]
    script_content += f"for xs_name in [{', '.join(material_names)}]:\n"
    script_content += textwrap.dedent("""
        mat = openmc.Material(name=xs_name)
        mat.set_density('macro', 1.0)
        mat.add_macroscopic(xs_name)
        materials_dict[xs_name] = mat

    materials_file = openmc.Materials(materials_dict.values())
    materials_file.cross_sections = mgxs_filepath
    
    materials_filepath = 'materials.xml'
    materials_file.export_to_xml(materials_filepath)
    print(f"--> Materials file saved to {{materials_filepath}}")
    """)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(script_content)

def main():
    parser = argparse.ArgumentParser(
        description="Convert a DOTB-style MGXS file to an OpenMC Python script."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input cross section text file."
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path to write the generated OpenMC Python script."
    )
    args = parser.parse_args()

    try:
        materials_data, num_groups = parse_dotb_file(args.input_file)
        generate_openmc_script(materials_data, num_groups, args.output_file)
        print(f"\\nSuccess! Generated OpenMC script at: {args.output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    main()
