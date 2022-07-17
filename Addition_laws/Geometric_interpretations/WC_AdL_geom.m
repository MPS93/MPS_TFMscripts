//*** Script intented to study and try to reproduce the geometric interpretation of the add. law. ***//

load "../../functions.m";
//load Ew_AddL_geom 

printf "\n === GEOMETRIC INTERPRETATION OF ADDITION LAW FOR WEIERSTRASS CURVES === \n";

// Example 1 of handpicked curve, for a walkthrough of the process 
BF := Rationals();
R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
R2<X,Y,Z> := ProjectiveSpace(R1,2);
O := [0,1,0];
lst := [27, -23321/132, 98, -1440071/1089, -2401];
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6:= lst[5];

P1 := [144/25, -147503/1375,1]; P2 := [-26/33, -1318/33, 1];
X1 := P1[1]; Y1 := P1[2]; Z1 := P1[3];
X2 := P2[1]; Y2 := P2[2]; Z2 := P2[3];

l := (Z1*Y2 - Y1*Z2)/(X2*Z1-X1*Z2); n := (Y1*X2 - Y2*X1)/(X2*Z1-X1*Z2);

C := Curve(R2, (Y^2)*Z + a1*X*Y*Z + a3*Y*Z^2 - X^3-a2*(X^2)*Z-a4*X*Z^2 - a6*Z^3);
L := Curve(R2,Y - l*X - n*Z);

pts := Points(Intersection(C,L)); R_lst := pts diff Seqset([R2!P1,R2!P2]);

R := R_lst[1]; Xr := R[1]; Zr := R[3];
Lp := Curve(R2,X*Zr - Xr*Z);
pts2 := Points(Intersection(C,Lp)); Rp_lst := pts2 diff Seqset([R2!R,R2!O]);
Rp := Rp_lst[1];




printf "\n Let's see how the geometric interpretation works for some numerical examples: \n";

printf "\n Lets consider the elliptic curve with parameters a1 = %o, a2 = %o, a3 = %o, a4 = %o, a6 = %o over the rationals. So the following curve: \n\n \t%o \n",a1,a2,a3,a4,a6,C;

printf "\n We also have two points on the curve P = (%o:%o:%o) and Q = (%o:%o:%o). Let's also consider the point at infinity O = (%o:%o:%o). \n",P1[1],P1[2],P1[3],P2[1],P2[2],P2[3],O[1],O[2],O[3];

printf "\n Now we have all we need to build the line passing through P and Q: \n\n \t%o \n", L;

printf "\n The points in the intersection of the curve and the line are the following: \n %o \n", pts;

printf "\n And so we can see that R is (%o:%o:%o). We can do some easy comprobations: is this point in the curve? \n\t%o \n And is it in the line? \n\t%o \n",R[1],R[2],R[3], R in C, R in L;

printf "\n Time to build the line passing through R and O: \n\n \t%o \n", Lp;

printf "\n The points in the intersection of the curve and this line are the following: \n %o \n", pts2;

printf "\n And so we can see that Rp = P3 is (%o:%o:%o). \n",Rp[1],Rp[2],Rp[3];

printf "\n Thus according to the geometric interpretation, the addition point is P3 = (%o:%o:%o). And indeed it is as we have checked with the addition law. Last checks: is this point in the curve? \n\t%o \n And is it in Lp? \n\t%o \n",Rp[1],Rp[2],Rp[3], Rp in C, Rp in Lp;



// Example 2 of a handpicked curve
BF := FiniteField(7); 
R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
R2<X,Y,Z> := ProjectiveSpace(R1,2);
O := [0,1,0];
lst := [ BF | 2, 0, 0, 4, 5 ];
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6:= lst[5];

P1 := [ 2, 0, 1 ];  P2:= [ 3, 2, 1 ]; mP2 := [3,6,1];

PP := Ew_AddL_geom(BF,lst,P1,P1);
QmQ := Ew_AddL_geom(BF,lst,P2,mP2);
PO := Ew_AddL_geom(BF,lst,P1,O);
PQ := Ew_AddL_geom(BF,lst,P1,P2);

printf "\n\n Now we have 2 more examples: \n --> First one with parameters a1 = %o, a2 = %o, a3 = %o, a4 = %o, a6 = %o over the GF(7). And the two points on the curve: P = (%o:%o:%o) and Q = (%o:%o:%o).\n",a1,a2,a3,a4,a6,P1[1],P1[2],P1[3],P2[1],P2[2],P2[3];

printf "\n Here is the result of doubling a point:  P + P = (%o:%o:%o). \n", PP[1],PP[2],PP[3];
printf "\n Difference of points:  Q + (-Q) = (%o:%o:%o). \n", QmQ[1],QmQ[2],QmQ[3];
printf "\n Adding a point to O:  P + O = (%o:%o:%o). \n", PO[1],PO[2],PO[3];
printf "\n And adding two different points:  P + Q = (%o:%o:%o). \n", PQ[1],PQ[2],PQ[3];



// Example 3 from safecurves
BF := FiniteField(109454571331697278617670725030735128145969349647868738157201323556196022393859); 
R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
R2<X,Y,Z> := ProjectiveSpace(R1,2);
O := [0,1,0];
a := -3; 
b := 107744541122042688792155207242782455150382764043089114141096634497567301547839;
lst := [ BF | a,b ];
a4 := lst[1]; a6:= lst[2];

P := [ BF|82638672503301278923015998535776227331280144783487139112686874194432446389503,43992510890276411535679659957604584722077886330284298232193264058442323471611, 1 ];

PP := Ew_AddL_geom(BF,lst,P,P);
PO := Ew_AddL_geom(BF,lst,P,O);
PQ := Ew_AddL_geom(BF,lst,P,PP);

printf "\n\n Lastly, an example from safecurves: a curve given in short Weierstrass form over a finite field of prime order. The parameters are  a = a4 = %o, b = a6 = %o. And a point on the curve is P = (%o:%o:%o).\n",a4,a6,P[1],P[2],P[3];

printf "\n Here is the result of doubling a point:  P + P = (%o:%o:%o). \n", PP[1],PP[2],PP[3];
printf "\n Adding a point to O:  P + O = (%o:%o:%o). \n", PO[1],PO[2],PO[3];
printf "\n And adding two different points:  P + 2P = (%o:%o:%o). \n", PQ[1],PQ[2],PQ[3];


quit;
