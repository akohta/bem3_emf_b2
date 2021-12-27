// -------------
h=0.25; // hight of cone ( base plane to top )
rt=0.005; // radius of top (sphere shape)
a=25.0/180.0*Pi; // angle of generatrix line
xs=0.0; // shift center
ys=0.0;
zs=-0.5*h;

nt=6; // division number of top sphere radius direction
pt=1.1; // progression number for top sphere raqdius direction
nr=6; // division number of circumference direction (per Pi/4)
ns=26; // division number of cone surface
ps=1.1; // progression number for cone surface
c=0.65; // coefficient inner division point on base plane
nb=4; // division number of base plane radial element
pb=1.1; // progression number for base plane radial element
// -------------

// constant
pp2=Pi/2.0;
ht=h+rt*(1.0/Sin(a)-1.0);
hc=h-rt;
rb=ht*Tan(a);
ri=c*rb;
// point data
Point(1)={xs,ys,zs};
Point(2)={xs+rb,ys,zs};
Point(3)={xs,ys,zs+h-rt};
Point(4)={xs,ys,zs+h};
Point(5)={xs+rt*Cos(a),ys,zs+h-rt+rt*Sin(a)};
Point(6)={ri+xs,0.0+ys,0.0+zs}; // base circle division point on x-axis
Point(7)={0.0+xs,ri+ys,0.0+zs}; // base circle division point on y-axis
Point(8)={-ri/Sqrt(2.0)+xs,-ri/Sqrt(2.0)+ys,0+zs}; // center of inner circle on x-y plane
// line data
Circle(1)={4,3,5};
Line(2)={5,2};
Line(3)={2,6};
Circle(4)={6,8,7};

// ---- cone surface ----
// transfinite settings for cone surface
Transfinite Line{1} = nt Using Progression pt;
Transfinite Line{2} = ns Using Progression ps;
// extrude
o0[]=Extrude{{0,0,1},{xs,ys,zs},pp2}{Line{1,2}; Layers{nr}; Recombine;};
o1[]=Extrude{{0,0,1},{xs,ys,zs},pp2}{Line{o0[0],o0[3]}; Layers{nr}; Recombine;};
o2[]=Extrude{{0,0,1},{xs,ys,zs},pp2}{Line{o1[0],o1[3]}; Layers{nr}; Recombine;};
o3[]=Extrude{{0,0,1},{xs,ys,zs},pp2}{Line{o2[0],o2[3]}; Layers{nr}; Recombine;};

// ---- base plane ----
// line 
b0[]={3,4};
b1[]=Rotate{{0,0,1},{xs,ys,zs}, pp2*1.0} {Duplicata{Line{b0[]};}};
b2[]=Rotate{{0,0,1},{xs,ys,zs}, pp2*2.0} {Duplicata{Line{b0[]};}};
b3[]=Rotate{{0,0,1},{xs,ys,zs}, pp2*3.0} {Duplicata{Line{b0[]};}};
b4[]=b0[];
// surface
sb1=newreg; Line Loop(sb1)={b0[0],b0[1],-b1[0],-o0[5]}; Surface(sb1)={sb1};
sb2=newreg; Line Loop(sb2)={b1[0],b1[1],-b2[0],-o1[5]}; Surface(sb2)={sb2};
sb3=newreg; Line Loop(sb3)={b2[0],b2[1],-b3[0],-o2[5]}; Surface(sb3)={sb3};
sb4=newreg; Line Loop(sb4)={b3[0],b3[1],-b4[0],-o3[5]}; Surface(sb4)={sb4};
sb5=newreg; Line Loop(sb5)={-b0[1],-b1[1],-b2[1],-b3[1]}; Surface(sb5)={sb5};
// transfinite
Transfinite Line{b0[0],b1[0],b2[0],b3[0]}=nb Using Progression pb;
Transfinite Line{b0[1],b1[1],b2[1],b3[1]}=nr+1;
Transfinite Surface{sb1,sb2,sb3,sb4,sb5};
Recombine Surface{sb1,sb2,sb3,sb4,sb5};

Physical Surface(1)={o0[1],o1[1],o2[1],o3[1],o0[4],o1[4],o2[4],o3[4],sb1,sb2,sb3,sb4,sb5}; // Domain 1
Physical Surface(99)={-o0[1],-o1[1],-o2[1],-o3[1],-o0[4],-o1[4],-o2[4],-o3[4],-sb1,-sb2,-sb3,-sb4,-sb5}; // Domain 0

