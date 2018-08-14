from fenics import compile_extension_module, MPI, HDF5File, mpi_comm_world, dolfin_version
from hashlib import sha1
from os import path

class HDF5Link:
    """Helper class for creating links in HDF5-files."""
    cpp_link_module = None
    def __init__(self):
        cpp_link_code = '''
        #include <hdf5.h>
        void link_dataset(const MPI_Comm comm,
                          const std::string hdf5_filename,
                          const std::string link_from,
                          const std::string link_to, bool use_mpiio)
        {
            hid_t hdf5_file_id = HDF5Interface::open_file(comm, hdf5_filename, "a", use_mpiio);
            herr_t status = H5Lcreate_hard(hdf5_file_id, link_from.c_str(), H5L_SAME_LOC,
                                link_to.c_str(), H5P_DEFAULT, H5P_DEFAULT);
            dolfin_assert(status != HDF5_FAIL);

            HDF5Interface::close_file(hdf5_file_id);
        }
        '''

        self.cpp_link_module = compile_extension_module(cpp_link_code, additional_system_headers=["dolfin/io/HDF5Interface.h"])

    def link(self, hdf5filename, link_from, link_to):
        "Create link in hdf5file."
        use_mpiio = MPI.size(mpi_comm_world()) > 1
        self.cpp_link_module.link_dataset(mpi_comm_world(), hdf5filename, link_from, link_to, use_mpiio)


def save_hdf5(fullname, field_name, data, timestep, hdf5_link):
    # Create "good enough" hash. This is done to avoid data corruption when restarted from
    # different number of processes, different distribution or different function space
    local_hash = sha1()
    local_hash.update(str(data.function_space().mesh().num_cells()))
    local_hash.update(str(data.function_space().ufl_element()))
    local_hash.update(str(data.function_space().dim()))
    local_hash.update(str(MPI.size(mpi_comm_world())))

    # Global hash (same on all processes), 10 digits long
    global_hash = MPI.sum(mpi_comm_world(), int(local_hash.hexdigest(), 16))
    global_hash = str(int(global_hash%1e10)).zfill(10)

    # Open HDF5File
    if not path.isfile(fullname):
        datafile = HDF5File(mpi_comm_world(), fullname, 'w')
    else:
        datafile = HDF5File(mpi_comm_world(), fullname, 'a')


    # Write to hash-dataset if not yet done
    if not datafile.has_dataset(global_hash) or not datafile.has_dataset(global_hash+"/"+field_name):
        datafile.write(data, str(global_hash)+"/"+field_name)


    if not datafile.has_dataset("Mesh"):
        datafile.write(data.function_space().mesh(), "Mesh")

    # Write vector to file
    datafile.write(data.vector(), field_name+str(timestep)+"/vector")

    # HDF5File.close is broken in 1.4, but fixed in dev.
    if dolfin_version() != "1.4.0":
        datafile.close()
    del datafile

    # Link information about function space from hash-dataset
    hdf5filename = str(global_hash)+"/"+field_name+"/%s"
    field_name_current = "%s%s" % (field_name, str(timestep)) +"/%s"
    for l in ["x_cell_dofs", "cell_dofs", "cells"]:
        hdf5_link(fullname, hdf5filename % l, field_name_current % l)
