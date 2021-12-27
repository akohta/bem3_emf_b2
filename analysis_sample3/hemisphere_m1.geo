//------------------------------------------------------------
n =    14;  // division number of 1/4 circle 
r = 0.500;  // radius of sphere 
sx= 0.000;  // x-component of translation vector
sy= 0.000;  // y-component of translation vector
sz=-0.01;  // z-component of translation vector
//------------------------------------------------------------

pp2=Pi/2.0;
l  =2.0/Sqrt(3.0)*r;
lp2=0.5*l;
s2l=0.5*Sqrt(2.0)*l;
sp2=r/Sqrt(2.0);

// basic point and line
Point( 1)={0+sx,0+sy,0+sz}; // origin
Point( 2)={r+sx,0+sy,0+sz}; // x0
Point( 3)={0+sx,r+sy,0+sz}; // y0
Point( 4)={s2l+sx,  0+sy, lp2+sz}; // xp
Point( 5)={  0+sx,s2l+sy, lp2+sz}; // yp
Point( 6)={s2l+sx,  0+sy,-lp2+sz}; // xm
Point( 7)={  0+sx,s2l+sy,-lp2+sz}; // ym
Point( 8)={s2l+sx,  0+sy,   0+sz}; // xs
Point( 9)={  0+sx,s2l+sy,   0+sz}; // ys
Point(10)={-sp2+sx,-sp2+sy, 0+sz};// shift of center   

Circle(1)={2,1,3}; // x0 to y0 
Circle(2)={3,1,5}; // y0 to yp
Circle(3)={5,1,4}; // yp to xp
Circle(4)={6,1,7}; // xm to ym
Circle(5)={7,1,3}; // ym to y0
Circle(6)={8,10,9}; // xs to ys
Line(7)={9,3}; // ys to y0

// upper sphere ( z > 0 )
u0[]={1,2,3};
u1[]=Rotate{{0,0,1},{0+sx,0+sy,0+sz},pp2*1.0} {Duplicata{Line{1,2,3};}};
u2[]=Rotate{{0,0,1},{0+sx,0+sy,0+sz},pp2*2.0} {Duplicata{Line{1,2,3};}};
u3[]=Rotate{{0,0,1},{0+sx,0+sy,0+sz},pp2*3.0} {Duplicata{Line{1,2,3};}};
Line Loop(1)={u0[0],u0[1],u0[2],-u3[1]};     Ruled Surface(1)={1} In Sphere{1};
Line Loop(2)={u1[0],u1[1],u1[2],-u0[1]};     Ruled Surface(2)={2} In Sphere{1};
Line Loop(3)={u2[0],u2[1],u2[2],-u1[1]};     Ruled Surface(3)={3} In Sphere{1};
Line Loop(4)={u3[0],u3[1],u3[2],-u2[1]};     Ruled Surface(4)={4} In Sphere{1};
Line Loop(5)={-u0[2],-u1[2],-u2[2],-u3[2]};  Ruled Surface(5)={5} In Sphere{1};

// separated plane
s0[]={6,7};
s1[]=Rotate{{0,0,1},{0+sx,0+sy,0+sz},pp2*1.0} {Duplicata{Line{6,7};}};
s2[]=Rotate{{0,0,1},{0+sx,0+sy,0+sz},pp2*2.0} {Duplicata{Line{6,7};}};
s3[]=Rotate{{0,0,1},{0+sx,0+sy,0+sz},pp2*3.0} {Duplicata{Line{6,7};}};
Line Loop(6)={s0[0],s0[1],-u0[0],-s3[1]};    Plane Surface(6)={6};
Line Loop(7)={s1[0],s1[1],-u1[0],-s0[1]};    Plane Surface(7)={7};
Line Loop(8)={s2[0],s2[1],-u2[0],-s1[1]};    Plane Surface(8)={8};
Line Loop(9)={s3[0],s3[1],-u3[0],-s2[1]};    Plane Surface(9)={9};
Line Loop(10)={-s0[0],-s1[0],-s2[0],-s3[0]}; Plane Surface(10)={10};

// lower sphere ( z < 0 )
l0[]={4,5};
l1[]=Rotate{{0,0,1},{0+sx,0+sy,0+sz},pp2*1.0} {Duplicata{Line{4,5};}};
l2[]=Rotate{{0,0,1},{0+sx,0+sy,0+sz},pp2*2.0} {Duplicata{Line{4,5};}};
l3[]=Rotate{{0,0,1},{0+sx,0+sy,0+sz},pp2*3.0} {Duplicata{Line{4,5};}};
Line Loop(11)={l0[0],l0[1],-u0[0],-l3[1]};   Ruled Surface(11)={11} In Sphere{1};
Line Loop(12)={l1[0],l1[1],-u1[0],-l0[1]};   Ruled Surface(12)={12} In Sphere{1};
Line Loop(13)={l2[0],l2[1],-u2[0],-l1[1]};   Ruled Surface(13)={13} In Sphere{1};
Line Loop(14)={l3[0],l3[1],-u3[0],-l2[1]};   Ruled Surface(14)={14} In Sphere{1};
Line Loop(15)={-l0[0],-l1[0],-l2[0],-l3[0]}; Ruled Surface(15)={15} In Sphere{1};


// definie domain
/*
Physical Surface(1)={1,2,3,4,5,6,7,8,9,10}; // Domain 1 ( upper sphere )
Physical Surface(2)={11,12,13,14,15,-6,-7,-8,-9,-10}; // Domain 2 ( lower sphere )
Physical Surface(99)={-1,-2,-3,-4,-5,-11,-12,-13,-14,-15}; // Domain 0 ( opened region )
*/
// lower hemi-sphere only
Physical Surface( 1)={ 11, 12, 13, 14, 15,-6,-7,-8,-9,-10}; // Domain 1 
Physical Surface(99)={-11,-12,-13,-14,-15, 6, 7, 8, 9, 10}; // Domain 0 ( opened region )

// transfinite interpolation 
Transfinite Line{u0[0],u1[0],u2[0],u3[0],
		   u0[2],u1[2],u2[2],u3[2],
		   s0[0],s1[0],s2[0],s3[0],
		   l0[0],l1[0],l2[0],l3[0]}=n+1;
Transfinite Line{u0[1],u1[1],u2[1],u3[1],
		   l0[1],l1[1],l2[1],l3[1]}=n/2+1;
Transfinite Line{s0[1],s1[1],s2[1],s3[1]}=n/3+1;
Transfinite Surface{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

// combine triangle element 
Recombine Surface{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
