In order to clip LA geometry you will have to execute the three steps:

#------------------------------------------------------------------------------------
## Files "1_generate_NU_vtu_files_from_xdmf_file.py" and "1_generate_U_vtu_files_from_xdmf_file.py" are paraview/python files that will
generate vtu files that will contain velocity or viscosity fields.

Inside those files you will have to write the correct "input_path" (where is located your "nu.xdmf" or "u.xdmf") and "output_path" (where you want to save your vtu files)

input --> path of the nu.xdmf or u.xdmf
output --> nu_XXX.vtu or u_XXX.vtu
Terminal: pvpython 1_generate_NU_vtu_files_from_xdmf_file
#------------------------------------------------------------------------------------
## To generate the clipped geometry it is necessary to use "2_extract_left_atrium_Volume.py"
having in the same directory this other file "vmtkpointselector.py", which is necessary for
selecting the point in the LAA to generate another branch of the center_lines and then clip 
the LAA.

### For running this script it is necessary to activate the "morphman" environment
### It is necessary to comment the line 464: 
"surface = vtk_convert_unstructured_grid_to_polydata(surface)"

This line is located in the fuction "vtk_compute_threshold", 
located in the line 428, in the script "vtk_wrapper.py" in the morphman environment,
located in the folder: morphMan/morphman/common/vtk_wrapper.py.


This line converts (cast) any type of file to vtp, but in this case we need to keep our files
in vtu format!!!

Input --> path to the surface mesh (vtp) uncapped: "C_XXX_remeshed_surface_EXAMPLE.vtp"
Input --> path to the volume meshes (several vtu meshes): "nu_relative_XXX.vtu"
input --> Three coordinates for a point located in the LAA (be careful with giving the correct point)

output --> several vtu files of the LAA/LA_body/LA_without_PV_MV located in the same folder
where you have your .vtp amd .vtu files

Terminal: conda activate morphman
(morphman) sergionio@serginio-computers:~/Downloads/test_nu/200_tsteps$ python3 extract_left_atrium_Volume.py --laa 6.708876 54.710644 21.220324 --includes_laa_and_la_body 1
#------------------------------------------------------------------------------------
## For visualizing the results in Paraview, you can upload them by using:

Terminal: paraview --data=something..vtu

The point (.) means that you have several files like this and Paraview will do the magic
uploading everything. So you are skipping the step of uploading one by one.
#------------------------------------------------------------------------------------
## For computing statistics data for example from the viscosity field (scalar field),
you can use the script "3_paraview_compute_quartiles_4.py":

Note --> This script works fine for scalar data like pressure, viscosity or shear rate
It wasn't tested in vectorial data like velocity
input --> Path for all vtu files from the LAA and LA_body
output -->Two text files with the values of the min, max, 25%, 50%, 75% of the data and the
average value
Terminal: pvpython 3_paraview_compute_quartiles_4.py



