# bem3_emf_b2
This is the three-dimensional electromagnetic field analysis program for arbitrary objects irradiated by arbitrary beams. 
This is the extension of 'bem3_emf_b1' using iterative solution. 
This program can analyze multiple scattering between objects with less memory than 'bem3_emf_b1'.
Intel Math Kernel Library is required. 
Gmsh is used for create a mesh data for object. 
The electromagnetic field analysis program "multi_fbeam" is used for analyze incident field.


## Usage of example code  
1. type 'make' command to compile  
   The executable d3b2_create_matrix, d3b2_bv_solver, example1.out, example2.out are created. 
   The executable 'd3b2_create_matrix' is the solver of boundary integral equations, it outputs coefficient matrices and its inverse matrices. 
   The executable 'd3b2_bv_solver' is the sovler for boundary value, it analyzes multiple scattering using iterative solution if defined multipe objests.
   The example1.out is the executable of source code example1.c, it shows a simplest example of usage. 
   The example2.out is the execubable of source code example2.c, it shows a example of electromagnetic field intensity analysis.

2. type './d3b2_create_matrix' with arguments of medium datafile name, mesh datafile name, output object file name, rotation and translation settings ( optional ).  
   For example, './d3b2_create_matrix medium_data.txt cone_m1.msh cone_m1_md_xpi.obj 1 0 0 3.14159265 0 0 0'. 
   In this case, the rotation axis is x-axis, the rotation angle is 3.14159265 ( Rodrigues' rotation formula is used ), the translation vector is zero-vector. 
   The medium_data.txt is the sample of medium datafile, one medium is defined in it. 
   The domain numbers are assigned to the medium from 1 in order. 
   The cone_m1.msh is the sample of mesh datafile, it is a cone object. 
   It was created by Gmsh geometry file cone_m1.geo in mesh_sample folder. 
   The cone_m1_image.pdf is the visualization result of the cone_m1.msh ( using Gmsh ). 
   The ipw.txt is the sample of incident field datafile, a plane-wave is defined in it. Please refer to the "multi_fbeam" for detail.
   The d3b2_create_matrix outputs object datafile, coefficient matrices and its inverve matrices ( binary file with the .cmat file extention ) with the specified datafile name. 
   
3. type './d3b2_bv_solver' with arguments of object setting datafile name, output datafile name.  
   For example, './d3b2_bv_solver object_settings.txt ex.dat'. 
   The object_settings.txt is the sample of object setting datafile, five cone objects defined in it. 
   The objects are set by using the outputed object file by d3b2_create_matrix and a additional translation vector. 
   The setting of reset_flag_of_incident_field in object setting datafile can reset the incident field. 
   It can change the incident field, except the wavelength and the refractive index.  
   
4. type './example1.out' with a argument of datafile name outputed by d3b2_bv_solver.
   For example, './example1.out ex.dat'. 
   This executable calculates electromagnetic field, radiation force and torque.  
   
5. type './example2.out' with a arguemnt of datafile name outputed by d3b2_bv_solver.  
   For example, './example2.out ex.dat'. 
   This executable calculates electromagnetic field intensity distributions, outputs them to text files.
   The I_example2_logcb.pdf is the visualization result of intensity distributions, created by the Gnuplot script gscript_example2_logcb.plt.

Please see 'd3b2_src/bem3_emf_b2.h' for detail of functions. 
The main parts of the code are parallelized by using OpenMP. 
The number of threads is controlled by the environment variable 'OMP_NUM_THREADS'.

## Verification  
The verification results are in folder verification. 
These are the analysis results of two cone objects by using the two different method.
The one is using the iterative solution, the other is using the combined object ( in the folder combined_object ). 


## Reference
1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)
2. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)
3. The electromagnetic field analysis program [multi_fbeam](https://github.com/akohta/multi_fbeam/) 
4. The electromagnetic field analysis program [bem3_emf_b1](https://github.com/akohta/bem3_emf_b1/)
