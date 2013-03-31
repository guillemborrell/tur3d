#include "io.h"

herr_t dump_file(const char *filename, double *field, grid_t phys){
    hsize_t dims[3] = {phys.nz, phys.ny, phys.nx};
    herr_t h5err;
    hid_t file = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    hid_t dset, dspace;

    dspace = H5Screate_simple(3, dims, NULL);
    dset = H5Dcreate(file,"value",H5T_IEEE_F32LE,dspace,
                     H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    h5err = H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,field);
    h5err = H5Dclose(dset);
    h5err = H5Sclose(dspace);
    h5err = H5Fclose(file);
}
