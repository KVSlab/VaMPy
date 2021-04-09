The following is a modification to the Aneurysm-workflow of KVSLab. The changes incorporated pertain to the file I/O structure and post-processing. These modifications results in a significant storage requirements For example, a test simulation resulted in the storage requirement dropping to 2 GB from 40 GB earlier.

******Modifications To File I/O************

These changes were made to the Artery.py file.

1 - The first change was converting the output timeseries data from XDMFile format to HDF5File format. The changes are reflected from lines 335 onwards. The code for writing the data in HDF5File format is as follows:

folder = "results_folder" #folder is a local variable

if MPI.rank(MPI.comm_world) == 0:
    counter = 1
    to_check = path.join(folder, "data", "%d")
    while path.isdir(to_check % counter):
        counter += 1

    if counter > 1:
        counter -= 1
    if not path.exists(path.join(to_check % counter, "hdf5")):
        makedirs(path.join(to_check % counter, "hdf5"))
else:
    counter = 0

counter = int(MPI.max(MPI.comm_world, counter))

common_path1 = path.join(folder, "data", str(counter), "hdf5")
mesh_path0 = path.join(common_path1, "mesh2.h5")
p_path = path.join(common_path1, "p.h5")
u0_path = path.join(common_path1, "u0.h5")
u1_path = path.join(common_path1, "u1.h5")
u2_path = path.join(common_path1, "u2.h5")

2 - The *.h5 files are written to "results_folder/data/hdf5" folder. The resulting folder structure is highlighted by the black dashed rectangle, whereas the old strucure is encompassed within the red dashed rectangle (see folder_strucure.png). Line 355 of the code indicates the main results_folder to which the data are written.

3 -  The current modification writes the data for the second cardiac cycle. This is done to save on storage, since, often the post-processing analysis are done for/from the second cycle (onwards). The modification is achieved through the 'tstep' variable. Lines 359 onwards reflect the modifications made to the code.

if tstep == 10000:
    viz_mesh = HDF5File(MPI.comm_world, mesh_path0, "w")
    viz_mesh.write(mesh,"mesh")
    viz_mesh.close()
    viz_p = HDF5File(MPI.comm_world, p_path, "w")
    viz_p.write(p_,"/pressure",tstep)
    viz_p.close()
    viz_u0 = HDF5File(MPI.comm_world, u0_path, "w")
    viz_u0.write(u_[0],"/velocity",tstep)
    viz_u0.close()
    viz_u1 = HDF5File(MPI.comm_world, u1_path, "w")
    viz_u1.write(u_[1],"/velocity", tstep)
    viz_u1.close()
    viz_u2 = HDF5File(MPI.comm_world, u2_path, "w")
    viz_u2.write(u_[2],"/velocity", tstep)
    viz_u2.close()
elif tstep>10000 and tstep % store_data == 0:
    viz_p = HDF5File(MPI.comm_world, p_path, "a")
    viz_p.write(p_,"/pressure",tstep)
    viz_p.close()
    viz_u0 = HDF5File(MPI.comm_world, u0_path, "a")
    viz_u0.write(u_[0],"/velocity", tstep)
    viz_u0.close()
    viz_u1 = HDF5File(MPI.comm_world, u1_path, "a")
    viz_u1.write(u_[1],"/velocity", tstep)
    viz_u1.close()
    viz_u2 = HDF5File(MPI.comm_world, u2_path, "a")
    viz_u2.write(u_[2],"/velocity", tstep)
    viz_u2.close()

4 -  It is important to note that the value of the folder variable in the function ‘ NS_parameters.update’ must match with the value in line 335. Only then the resulting folder structure will be same as the one shown in ‘folder_structure.png’.

******Changes To Post-Processing************

1 - The mesh data written in the .h5 format during a parallel run results in the reordering of the degrees of freedom (DOF) of the cells/elements. This is reordering causing interpolation issues during post-processing stage.

2 -  To overcome DOF reordering issue generate the mesh in a serial mode and use this mesh in the post-processing script. The mesh gnerated in the serial run is in XDMFFile format. This is given in the main_results_folder/model_name/data/counter/VTK (see the red dashed rectangle highlighted folder structure in 'folder_structure.png').

3 - Copy the serial generated mesh.xdmf and mesh.h5 files to main_results_folder/data/counter/hdf5 folder before executing the post-processing scripts.

4 -  The files are then read in a usual manner in the post-processing script. Refer to the attached compute_wss_modified.py file for further information.
