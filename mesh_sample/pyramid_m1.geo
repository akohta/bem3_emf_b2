//------------------------
h=0.25; // hight of pyramid
a=30.0/180.0*Pi; // angle of face 
zs=-0.5*h;

nf=24;
pf=1.1;
nb=12;
//------------------------

th=2.0*Pi/4;
rb=h*Tan(a);

Point(1)={0,0,h+zs};
Point(2)={rb, rb*Tan(0.5*th),zs};
Point(3)={rb,-rb*Tan(0.5*th),zs};

Line(1)={1,3};
Line(2)={3,2};

s0[]={1,2};
s1[]=Rotate{{0,0,1},{0,0,0},th*1}{Duplicata{Line{1,2};}};
s2[]=Rotate{{0,0,1},{0,0,0},th*2}{Duplicata{Line{1,2};}};
s3[]=Rotate{{0,0,1},{0,0,0},th*3}{Duplicata{Line{1,2};}};

Transfinite Line{s0[0],s1[0],s2[0],s3[0]} = nf Using Progression pf;
Transfinite Line{s0[1],s1[1],s2[1],s3[1]} =nb;

Line Loop(1)={s0[0],s0[1],-s1[0]}; Surface(1)={1};
Line Loop(2)={s1[0],s1[1],-s2[0]}; Surface(2)={2};
Line Loop(3)={s2[0],s2[1],-s3[0]}; Surface(3)={3};
Line Loop(4)={s3[0],s3[1],-s0[0]}; Surface(4)={4};
Line Loop(5)={-s0[1],-s1[1],-s2[1],-s3[1]}; Surface(5)={5};

Transfinite Surface{1,2,3,4,5};
Recombine Surface{1,2,3,4,5};

Physical Surface(1)={1,2,3,4,5};
Physical Surface(99)={-1,-2,-3,-4,-5};

