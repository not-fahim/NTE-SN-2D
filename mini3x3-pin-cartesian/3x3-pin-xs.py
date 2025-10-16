
import openmc
import numpy as np

print("--- Generating OpenMC MGXS Library and Materials ---")

# Create a 2-group structure.
# NOTE: The energy boundaries here are placeholders. You should
# replace them with the actual energy boundaries for your problem.
# The boundaries should go from lowest to highest energy.
group_bounds = np.logspace(-5, 7, 3) # Example: 0.01 eV to 10 MeV
groups = openmc.mgxs.EnergyGroups(group_bounds)

# Initialize the library
mg_cross_sections_file = openmc.MGXSLibrary(groups)

# --------------------------------------------------------------------------
# Add Cross Section Data for Each Material
# --------------------------------------------------------------------------

# --- Material: fuel ---
fuel_xsdata = openmc.XSdata('fuel', groups)
fuel_xsdata.order = 0 # Isotropic scattering
fuel_xsdata.set_total([0.368077135, 0.651647585], temperature=294.)
fuel_xsdata.set_absorption([0.0335406364, 0.267179995], temperature=294.)
fuel_xsdata.set_nu_fission([0.0226011124, 0.48644535], temperature=294.)
fuel_xsdata.set_chi([1.0, 0.0], temperature=294.)

# The scattering matrix from the file is S(g_out, g_in).
# OpenMC's `set_scatter_matrix` expects the matrix indexed by S[g_in][g_out].
# Therefore, we use the transposed matrix from the original file.
# The shape must be (G, G, M+1), so we add a third dimension for M=0.
fuel_scatter_matrix = [[[0.333749271, 0.000787296057], [0.000390376409, 0.384062313]]]
fuel_scatter_matrix = np.array(fuel_scatter_matrix)
fuel_scatter_matrix = np.rollaxis(fuel_scatter_matrix, 0, 3)
fuel_xsdata.set_scatter_matrix(fuel_scatter_matrix, temperature=294.)
mg_cross_sections_file.add_xsdata(fuel_xsdata)


# --- Material: cool ---
cool_xsdata = openmc.XSdata('cool', groups)
cool_xsdata.order = 0 # Isotropic scattering
cool_xsdata.set_total([0.280097317, 1.22741586], temperature=294.)
cool_xsdata.set_absorption([0.0007629659, 0.0306393552], temperature=294.)
cool_xsdata.set_nu_fission([0.0, 0.0], temperature=294.)
cool_xsdata.set_chi([0.0, 0.0], temperature=294.)

# The scattering matrix from the file is S(g_out, g_in).
# OpenMC's `set_scatter_matrix` expects the matrix indexed by S[g_in][g_out].
# Therefore, we use the transposed matrix from the original file.
# The shape must be (G, G, M+1), so we add a third dimension for M=0.
cool_scatter_matrix = [[[0.245907516, 0.0334078255], [0.000525517043, 1.19625164]]]
cool_scatter_matrix = np.array(cool_scatter_matrix)
cool_scatter_matrix = np.rollaxis(cool_scatter_matrix, 0, 3)
cool_xsdata.set_scatter_matrix(cool_scatter_matrix, temperature=294.)
mg_cross_sections_file.add_xsdata(cool_xsdata)


# --- Material: gad ---
gad_xsdata = openmc.XSdata('gad', groups)
gad_xsdata.order = 0 # Isotropic scattering
gad_xsdata.set_total([0.388826464, 4.85902365], temperature=294.)
gad_xsdata.set_absorption([0.0461089161, 4.44586931], temperature=294.)
gad_xsdata.set_nu_fission([0.0156646537, 0.149578858], temperature=294.)
gad_xsdata.set_chi([1.0, 0.0], temperature=294.)

# The scattering matrix from the file is S(g_out, g_in).
# OpenMC's `set_scatter_matrix` expects the matrix indexed by S[g_in][g_out].
# Therefore, we use the transposed matrix from the original file.
# The shape must be (G, G, M+1), so we add a third dimension for M=0.
gad_scatter_matrix = [[[0.341855158, 0.000840844768], [0.00119297829, 0.412239303]]]
gad_scatter_matrix = np.array(gad_scatter_matrix)
gad_scatter_matrix = np.rollaxis(gad_scatter_matrix, 0, 3)
gad_xsdata.set_scatter_matrix(gad_scatter_matrix, temperature=294.)
mg_cross_sections_file.add_xsdata(gad_xsdata)


# --------------------------------------------------------------------------
# Export Library and Create Materials XML
# --------------------------------------------------------------------------

# Export the MGXS library to an HDF5 file
mgxs_filepath = 'mgxs.h5'
mg_cross_sections_file.export_to_hdf5(mgxs_filepath)
print(f"--> MGXS library saved to {{mgxs_filepath}}")

# Create a Materials XML file that links to the MGXS library
materials_dict = {}
for xs_name in ['fuel', 'cool', 'gad']:

    mat = openmc.Material(name=xs_name)
    mat.set_density('macro', 1.0)
    mat.add_macroscopic(xs_name)
    materials_dict[xs_name] = mat

materials_file = openmc.Materials(materials_dict.values())
materials_file.cross_sections = mgxs_filepath

materials_filepath = 'materials.xml'
materials_file.export_to_xml(materials_filepath)
print(f"--> Materials file saved to {{materials_filepath}}")
