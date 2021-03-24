// #pragma once

#include "data_structure_mod.h"
#include <fstream>
#include <vector>
#include <iostream>
#include "H5Cpp.h"
// #include "H5Cpp.h"
// #include "petsc_data_structure_mod.h"
using namespace H5;
using namespace std;
const H5std_string FILE_NAME("/opt/grids/legacy/hdf5/point_9600.h5");
const H5std_string DATASET_NAME("/1/local");
const int RANK_OUT = 2;
// H5File file(FILE_NAME, H5F_ACC_RDWR);
// DataSet dataset = file.openDataSet(DATASET_NAME);

int read_data()
{
   /*
    * Output buffer initialization.
    */
   /*
    * Try block to detect exceptions raised by any of the calls inside it
    */
   try
   {
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      Exception::dontPrint();
      /*
       * Open the specified file and the specified dataset in the file.
       */
      H5File file(FILE_NAME, H5F_ACC_RDONLY);
      Group group(file.openGroup("/1"));
      DataSet dataset = file.openDataSet(DATASET_NAME);
      //  Attribute attribute = file.openAttribute(DATASET_NAME);
      Attribute attribute = Attribute(group.openAttribute("local"));
      H5Aread(attribute.getId(), H5T_NATIVE_LLONG, &local_points);
      cout << local_points << endl;
      //    //ATTRIBUTE
      //    Attribute att(myObject->openAttribute("local"));

      //    /*
      //     * Get the class of the datatype that is used by the dataset.
      //     */
      H5T_class_t type_class = dataset.getTypeClass();
      //   /*
      //  * Get class of datatype and print message if it's an integer.
      //  */
      cout << typeid(type_class).name() << endl;
      /*
       * Get dataspace of the dataset.
       */
      DataSpace dataspace = dataset.getSpace();
      /*
       * Get the number of dimensions in the dataspace.
       */
      int rank = dataspace.getSimpleExtentNdims();
      /*
    //    * Get the dimension size of each dimension in the dataspace and
    //    * display them.
    //    */
      hsize_t dims_out[2];
      int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
      cout << "rank " << rank << ", dimensions " << (unsigned long)(dims_out[0]) << " x " << (unsigned long)(dims_out[1]) << endl;

      //   /*
      //    * Define hyperslab in the dataset; implicitly giving strike and
      //    * block NULL.
      //    */
      //   hsize_t      offset[2];   // hyperslab offset in the file
      //   hsize_t      count[2];    // size of the hyperslab in the file
      //   offset[0] = 1;
      //   offset[1] = 2;
      //   count[0]  = NX_SUB;
      //   count[1]  = NY_SUB;
      //   dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );
      //   /*
      //    * Define the memory dataspace.
      //    */
      //   hsize_t     dimsm[3];              /* memory space dimensions */
      //   dimsm[0] = NX;
      //   dimsm[1] = NY;
      //   dimsm[2] = NZ ;
      hsize_t dims[2] = {local_points, 30};
      DataSpace memspace(RANK_OUT, dims);
      double *data_out = new double[local_points * 30];
      dataset.read(data_out, PredType::NATIVE_DOUBLE, memspace, dataspace);
      for (int i = 0; i < local_points; i++)
      {
         point.x[i + 1] = data_out[30 * i + 1];
         point.y[i + 1] = data_out[30 * i + 2];
         point.nx[i + 1] = data_out[30 * i + 3];
         point.ny[i + 1] = data_out[30 * i + 4];
         point.min_dist[i + 1] = data_out[30 * i + 5];
         point.left[i + 1] = (int)data_out[30 * i + 6];
         point.right[i + 1] = (int)data_out[30 * i + 7];
         point.flag_1[i + 1] = (int)data_out[30 * i + 9];
         point.flag_2[i + 1] = (int)data_out[30 * i + 10];
         point.nbhs[i + 1] = (int)data_out[30 * i + 11];
         for (int r = 1; r <= point.nbhs[i + 1]; r++)
         {
            point.conn[i + 1][r] = (int)data_out[30 * i + 11 + r];
         }
      }
      return 1;
      //   /*
      //    * Define memory hyperslab.
      //    */
      //   hsize_t      offset_out[3];   // hyperslab offset in memory
      //   hsize_t      count_out[3];    // size of the hyperslab in memory
      //   offset_out[0] = 3;
      //   offset_out[1] = 0;
      //   offset_out[2] = 0;
      //   count_out[0]  = NX_SUB;
      //   count_out[1]  = NY_SUB;
      //   count_out[2]  = 1;
      //   memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );
      //   /*
      //    * Read data from hyperslab in the file into the hyperslab in
      //    * memory and display the data.
      //    */
      //   dataset.read( data_out, PredType::NATIVE_INT, memspace, dataspace );
      //   for (j = 0; j < NX; j++)
      //   {
      // for (i = 0; i < NY; i++)
      //    cout << data_out[j][i][0] << " ";
      // cout << endl;
      //   }
      //   /*
      //    * 0 0 0 0 0 0 0
      //    * 0 0 0 0 0 0 0
      //    * 0 0 0 0 0 0 0
      //    * 3 4 5 6 0 0 0
      //    * 4 5 6 7 0 0 0
      //    * 5 6 7 8 0 0 0
      //    * 0 0 0 0 0 0 0
      //    */
   } // end of try block
     // catch failure caused by the H5File operations
   catch (FileIException error)
   {
      error.printErrorStack();
      return -1;
   }
   // catch failure caused by the DataSet operations
   catch (DataSetIException error)
   {
      error.printErrorStack();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch (DataSpaceIException error)
   {
      error.printErrorStack();
      return -1;
   }
   // catch failure caused by the DataSpace operations
   catch (DataTypeIException error)
   {
      error.printErrorStack();
      return -1;
   }
   return 0; // successfully terminated
}
