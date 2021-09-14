// ------------------------------
nx=12; // sampling number of x
ny=12; // sampling number of y
nz= 5; // sampling number of z
wx =0.5; // width of x-direction
wy =0.5; // width of y-direction
wz1=0.2; // width of z-direction
wz2=0.2; // width of z-direction
tx=0.0;    // 
ty=0.0;    // (tx,ty,tz) 
tz=0.001; // translation vector
// ------------------------------

// base point 
Point( 1)={+0.5*wx+tx,+0.5*wy+ty,+wz1+tz};
Point( 2)={+0.5*wx+tx,+0.5*wy+ty,+0.0+tz};
Point( 3)={+0.5*wx+tx,+0.5*wy+ty,-wz2+tz};
// base line
Line( 1)={1,2};
Line( 2)={2,3};
// transfinite interpolation
Transfinite Line{1,2}=nz;

o0[]=Extrude{-wx,0,0}{Line{1}; Layers{nx}; Recombine;};
o1[]=Extrude{0,-wy,0}{Line{-1}; Layers{ny}; Recombine;};
o2[]=Extrude{0,-wy,0}{Line{o0[0]}; Layers{ny}; Recombine;};
o3[]=Extrude{-wx,0,0}{Line{o1[0]}; Layers{nx}; Recombine;};
o4[]=Extrude{0,-wy,0}{Line{-o0[3]}; Layers{ny}; Recombine;};
o5[]=Extrude{0,-wy,0}{Line{-o0[2]}; Layers{ny}; Recombine;};

p0[]=Extrude{-wx,0,0}{Line{ 2}; Layers{nx}; Recombine;};
p1[]=Extrude{0,-wy,0}{Line{-2}; Layers{ny}; Recombine;};
p2[]=Extrude{0,-wy,0}{Line{p0[0]}; Layers{ny}; Recombine;};
p3[]=Extrude{-wx,0,0}{Line{p1[0]}; Layers{nx}; Recombine;};
p4[]=Extrude{0,-wy,0}{Line{-p0[2]}; Layers{ny}; Recombine;};

Physical Surface( 1)={ o0[1], o1[1], o2[1], o3[1], o4[1], o5[1]}; // Domain 1
Physical Surface( 2)={ p0[1], p1[1], p2[1], p3[1], p4[1],-o5[1]}; // Domain 2 
Physical Surface(99)={-o0[1],-o1[1],-o2[1],-o3[1],-o4[1],-p0[1],-p1[1],-p2[1],-p3[1],-p4[1]}; // Domain 0 (opened)



