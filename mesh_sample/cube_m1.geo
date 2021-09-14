// ------------------------------
nx=12; // sampling number of x
ny=12; // sampling number of y
nz=5; // sampling number of z
wx=0.5; // width of x-direction
wy=0.5; // width of y-direction
wz=0.2; // width of z-direction
// ------------------------------

// base point 
Point( 1)={+0.5*wx,+0.5*wy,+0.5*wz};
Point( 2)={+0.5*wx,+0.5*wy,-0.5*wz};
// base line
Line( 1)={1,2};
// transfinite interpolation
Transfinite Line{1}=nz;

o0[]=Extrude {-wx,0,0}{Line{1}; Layers{nx}; Recombine;};
o1[]=Extrude {0,-wy,0}{Line{-1}; Layers{ny}; Recombine;};
o2[]=Extrude {0,-wy,0}{Line{o0[0]}; Layers{ny}; Recombine;};
o3[]=Extrude {0,-wy,0}{Line{-o0[2]}; Layers{ny}; Recombine;};
o4[]=Extrude {0,-wy,0}{Line{-o0[3]}; Layers{ny}; Recombine;};
o5[]=Extrude {-wx,0,0}{Line{o1[0]}; Layers{nx}; Recombine;};

Physical Surface( 1)={ o0[1], o1[1], o2[1], o3[1], o4[1], o5[1]}; // Domain 1
Physical Surface(99)={-o0[1],-o1[1],-o2[1],-o3[1],-o4[1],-o5[1]}; // Domain 0 (opened)


