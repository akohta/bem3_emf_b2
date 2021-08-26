//------------------------------------------
n=5; // basic sampling number 

h=0.25;          // height of cone
a=30.0/180.0*Pi; // angle of gereratrix line, cone angle (aperture) is 2*a
zs=-0.075-h;     // z direction shift

c=0.75;  // coefficient of inner division point on base plane
m=3;     // coefficient of generatrix line sampling number 
p=1.05;  // progression coefficient of sampling number on generatrix line

xd=0.0;    // duplucate model shift parameter
yd=0.0;    // 
zd=h+0.15; // 
thx=Pi; // rotation angle
//------------------------------------------

pp2=Pi/2.0;
pp4=Pi/4.0;
ta=Tan(a);
rb=h*ta;
ri=c*rb;
ric=ri*Cos(pp4);

// basic point
Point( 1)={0,0,0+zs}; // center of base circle
Point( 2)={0,0,h+zs}; // top of cone
Point( 3)={rb,0,0+zs}; // crosspoint of x-axis and generatrix line
Point( 4)={ri,0,0+zs}; // base circle division point on x-axis
Point( 5)={0,ri,0+zs}; // base circle division point on y-axis
Point( 6)={-ri/Sqrt(2.0),-ri/Sqrt(2.0),0+zs}; // center of inner circle on x-y plane

// basic line
Line(1)={2,3};
Line(2)={3,4};
Circle(3)={4,6,5}; 

// transfinite interpolation to cone surface 
Transfinite Line{1}=n*(m/Cos(a))+1 Using Progression p;

// cone surface except base plane
o0[]=Extrude {{0,0,1},{0,0,0}, pp2}{Line{1}; Layers{pp2*n}; Recombine;};
o1[]=Extrude {{0,0,1},{0,0,0}, pp2}{Line{o0[0]}; Layers{pp2*n}; Recombine;};
o2[]=Extrude {{0,0,1},{0,0,0}, pp2}{Line{o1[0]}; Layers{pp2*n}; Recombine;};
o3[]=Extrude {{0,0,1},{0,0,0}, pp2}{Line{o2[0]}; Layers{pp2*n}; Recombine;};
// duplicate
l1[]=Translate {xd,yd,zd} { Duplicata{ Line{1}; } };
Transfinite Line{l1[0]}=n*(m/Cos(a))+1 Using Progression p;
od0[]=Extrude {{0,0,1},{0,0,0}, pp2}{Line{l1[0]}; Layers{pp2*n}; Recombine;};
od1[]=Extrude {{0,0,1},{0,0,0}, pp2}{Line{od0[0]}; Layers{pp2*n}; Recombine;};
od2[]=Extrude {{0,0,1},{0,0,0}, pp2}{Line{od1[0]}; Layers{pp2*n}; Recombine;};
od3[]=Extrude {{0,0,1},{0,0,0}, pp2}{Line{od2[0]}; Layers{pp2*n}; Recombine;};

// base plane
b0[]={2,3};
b1[]=Rotate{{0,0,1},{0,0,0}, pp2*1.0} {Duplicata{Line{2,3};}};
b2[]=Rotate{{0,0,1},{0,0,0}, pp2*2.0} {Duplicata{Line{2,3};}};
b3[]=Rotate{{0,0,1},{0,0,0}, pp2*3.0} {Duplicata{Line{2,3};}};
b4[]={2,3}; // copy of b0[]
sb1=newreg; Line Loop(sb1)={b0[0],b0[1],-b1[0],-o0[2]}; Surface(sb1)={sb1};
sb2=newreg; Line Loop(sb2)={b1[0],b1[1],-b2[0],-o1[2]}; Surface(sb2)={sb2};
sb3=newreg; Line Loop(sb3)={b2[0],b2[1],-b3[0],-o2[2]}; Surface(sb3)={sb3};
sb4=newreg; Line Loop(sb4)={b3[0],b3[1],-b4[0],-o3[2]}; Surface(sb4)={sb4};
sb5=newreg; Line Loop(sb5)={-b0[1],-b1[1],-b2[1],-b3[1]}; Surface(sb5)={sb5};
// transfinite interpolation to base plane
Transfinite Line{b0[0],b1[0],b2[0],b3[0]}=n*c;
Transfinite Line{b0[1],b1[1],b2[1],b3[1]}=pp2*n+1;
Transfinite Surface{sb1,sb2,sb3,sb4,sb5};
Recombine Surface{sb1,sb2,sb3,sb4,sb5};
// duplicate
l2[]=Translate {xd,yd,zd} { Duplicata{ Line{2}; } };
l3[]=Translate {xd,yd,zd} { Duplicata{ Line{3}; } };
bd0[]={l2[0],l3[0]};
bd1[]=Rotate{{0,0,1},{0,0,0}, pp2*1.0} {Duplicata{Line{l2[0],l3[0]};}};
bd2[]=Rotate{{0,0,1},{0,0,0}, pp2*2.0} {Duplicata{Line{l2[0],l3[0]};}};
bd3[]=Rotate{{0,0,1},{0,0,0}, pp2*3.0} {Duplicata{Line{l2[0],l3[0]};}};
bd4[]={l2[0],l3[0]}; // copy of bd0[]
sdb1=newreg; Line Loop(sdb1)={ bd0[0], bd0[1],-bd1[0],-od0[2]}; Surface(sdb1)={sdb1};
sdb2=newreg; Line Loop(sdb2)={ bd1[0], bd1[1],-bd2[0],-od1[2]}; Surface(sdb2)={sdb2};
sdb3=newreg; Line Loop(sdb3)={ bd2[0], bd2[1],-bd3[0],-od2[2]}; Surface(sdb3)={sdb3};
sdb4=newreg; Line Loop(sdb4)={ bd3[0], bd3[1],-bd4[0],-od3[2]}; Surface(sdb4)={sdb4};
sdb5=newreg; Line Loop(sdb5)={-bd0[1],-bd1[1],-bd2[1],-bd3[1]}; Surface(sdb5)={sdb5};
Transfinite Line{bd0[0],bd1[0],bd2[0],bd3[0]}=n*c;
Transfinite Line{bd0[1],bd1[1],bd2[1],bd3[1]}=pp2*n+1;
Transfinite Surface{sdb1,sdb2,sdb3,sdb4,sdb5};
Recombine Surface{sdb1,sdb2,sdb3,sdb4,sdb5};


// define domain
Physical Surface(1)={o0[1],o1[1],o2[1],o3[1],
                        sb1,sb2,sb3,sb4,sb5}; // Domain 1
Physical Surface(2)={od0[1],od1[1],od2[1],od3[1],
                        sdb1,sdb2,sdb3,sdb4,sdb5}; // Domain 2
Physical Surface(99)={-o0[1],-o1[1],-o2[1],-o3[1],-sb1,-sb2,-sb3,-sb4,-sb5,
						-od0[1],-od1[1],-od2[1],-od3[1],-sdb1,-sdb2,-sdb3,-sdb4,-sdb5}; // Domain 0 (opend)
 

