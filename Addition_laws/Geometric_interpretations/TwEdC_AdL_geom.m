//*** Script intented to study and try to reproduce the geometric interpretation of the add. law. ***//

load "../../functions.m";
//load Ead_AddL_geom

printf "\n === GEOMETRIC INTERPRETATION OF ADDITION LAW FOR TWISTED EDWARDS CURVES === \n";

//** Symbolic part **//

BF := Rationals();
R1<a,d> := FunctionField(BF,2);
R2<X1,Y1,Z1,X2,Y2,Z2> := PolynomialRing(R1,6);
R3<X,Y,Z> := PolynomialRing(R2,3);


// The curve for P and Q
f_crv1 := (a*(X1^2) + (Y1^2))*Z1^2 - Z1^4 - d*(X1^2)*(Y1^2);
f_crv2 := (a*(X2^2) + (Y2^2))*Z2^2 - Z2^4 - d*(X2^2)*(Y2^2);
f_crv := (a*(X^2) + (Y^2))*Z^2 - Z^4 - d*(X^2)*(Y^2);

// Since P and Q are points on the curve, f_crv1 and f_crv2 need to be zero.
// So we need to define the ideal generated by the curve to be able to compute things modulo the curve.
I_crv := ideal<R2|f_crv1,f_crv2>; QI_crv := R2/I_crv;


// The conic for doubling a point 
cz_PP := X1*Z1*(Z1-Y1);
cxy_PP := d*(X1^2)*Y1-Z1^3;
cxz_PP := Z1*(Z1*Y1 - a*X1^2);

f_cnc_PP := cz_PP*(Z^2 + Y*Z)+cxy_PP*X*Y+cxz_PP*X*Z;

// The conic for adding to Op
cz_POp := -X1; cxy_POp := Z1; cxz_POp := Z1;

f_cnc_POp := cz_POp*(Z^2 + Y*Z)+cxy_POp*X*Y+cxz_POp*X*Z;

// The conic otherwise
cz_PQ := X1*X2*(Y1*Z2-Y2*Z1);
cxy_PQ := Z1*Z2*(X1*Z2-X2*Z1 + X1*Y2 - X2*Y1);
cxz_PQ := X2*Y2*(Z1^2)-X1*Y1*(Z2^2)+Y1*Y2*(X2*Z1-X1*Z2);

f_cnc_PQ := cz_PQ*(Z^2 + Y*Z)+cxy_PQ*X*Y+cxz_PQ*X*Z;


// Addition law for this projective space and for each of the 3 cases

// General addition law
A := Z1*Z2; B := A^2; C := X1*X2; D := Y1*Y2; E := d*C*D; F := B-E; G := B+E;
X3 := A*F*((X1+Y1)*(X2+Y2)-C-D);
Y3 := A*G*(D-a*C);
Z3 := F*G;

// Point doubling
X3_PP := Evaluate(X3,[X1,Y1,Z1,X1,Y1,Z1]);
Y3_PP := Evaluate(Y3,[X1,Y1,Z1,X1,Y1,Z1]);
Z3_PP := Evaluate(Z3,[X1,Y1,Z1,X1,Y1,Z1]);

// Adding to Op
X3_POp := Evaluate(X3,[X1,Y1,Z1,0,-1,1]);
Y3_POp := Evaluate(Y3,[X1,Y1,Z1,0,-1,1]);
Z3_POp := Evaluate(Z3,[X1,Y1,Z1,0,-1,1]);

//--> Check if Op, O1, O2, P1, P2 and -P3 are in the curve and in the conic
printf "\n Let's do some symbolic checks first: \n";
printf "\n For starters, we have the elliptic curve in projective form: \n\t%o \n",f_crv;
printf "\n Here is the addition law in the projective space P2: \n\tX3 = %o \n\tY3 = %o \n\tZ3 = %o \n", X3,Y3,Z3;
printf "\n Are Op and the infinity points on the curve? %o \n And P1, P2 are also on the curve, right? %o \n ", (Evaluate(f_crv,[0,-1,1]) in I_crv) and (Evaluate(f_crv,[0,1,0]) in I_crv) and (Evaluate(f_crv,[1,0,0]) in I_crv),(Evaluate(f_crv,[X1,Y1,Z1]) in I_crv) and (Evaluate(f_crv,[X2,Y2,Z2]) in I_crv);
printf "\n Is P3 = (X3 : Y3 : Z3) in the curve? %o \n And is -P3 = (-X3 : Y3 : Z3) also in the curve? %o \n",Evaluate(f_crv,[X3,Y3,Z3]) in I_crv,Evaluate(f_crv,[-X3,Y3,Z3]) in I_crv;

printf "\n The conic has different coefficients for 3 different cases: \n";
printf "\n --> In the case of point doubling, the conic has this form: \n\t %o \n", f_cnc_PP;
printf "\n Are Op and the infinity points on the conic? %o \n And are P1, P2 also on it? %o \n ", (Evaluate(f_cnc_PP,[0,-1,1]) in I_crv) and (Evaluate(f_cnc_PP,[0,1,0]) in I_crv) and (Evaluate(f_cnc_PP,[1,0,0]) in I_crv),Evaluate(f_cnc_PP,[X1,Y1,Z1]) in I_crv;
printf "\n We expect -P3 = (-X3 : Y3 : Z3) to be in the conic. Does it hold? %o \n But in general P3 = (X3 : Y3 : Z3) shouldn't be in the conic. Is it? %o \n",Evaluate(f_cnc_PP,[-X3_PP,Y3_PP,Z3_PP]) in I_crv,Evaluate(f_cnc_PP,[X3_PP,Y3_PP,Z3_PP]) in I_crv;

printf "\n --> In the case of adding a point to Op = (0 : -1 : 1), the conic has this form: \n\t %o \n", f_cnc_POp;
printf "\n Are Op and the infinity points on the conic? %o \n And are P1, P2 also on it? %o \n ", (Evaluate(f_cnc_POp,[0,-1,1]) in I_crv) and (Evaluate(f_cnc_POp,[0,1,0]) in I_crv) and (Evaluate(f_cnc_POp,[1,0,0]) in I_crv), Evaluate(f_cnc_POp,[X1,Y1,Z1]) in I_crv;
printf "\n We expect -P3 = (-X3 : Y3 : Z3) to be in the conic. Does it hold? %o \n But in general P3 = (X3 : Y3 : Z3) shouldn't be in the conic. Is it? %o \n",Evaluate(f_cnc_POp,[-X3_POp,Y3_POp,Z3_POp]) in I_crv,Evaluate(f_cnc_POp,[X3_POp,Y3_POp,Z3_POp]) in I_crv;

printf "\n --> In any other case, the conic has this form: \n\t %o \n", f_cnc_PQ;
printf "\n Are Op and the infinity points on the conic? %o \n And are P1, P2 also on it? %o \n ", (Evaluate(f_cnc_PQ,[0,-1,1]) in I_crv) and (Evaluate(f_cnc_PQ,[0,1,0]) in I_crv) and (Evaluate(f_cnc_PQ,[1,0,0]) in I_crv), (Evaluate(f_cnc_PQ,[X1,Y1,Z1]) in I_crv) and (Evaluate(f_cnc_PQ,[X2,Y2,Z2]) in I_crv);
printf "\n We expect -P3 = (-X3 : Y3 : Z3) to be in the conic. Does it hold? %o \n But in general P3 = (X3 : Y3 : Z3) shouldn't be in the conic. Is it? %o \n",Evaluate(f_cnc_PQ,[-X3,Y3,Z3]) in I_crv,Evaluate(f_cnc_PQ,[X3,Y3,Z3]) in I_crv;

printf "\n Thus we see that the Op, O1, O2, P1, P2 and -P3 are in both the curve and the conic in any case. So the geometric interpretation holds. \n";
printf "Op: %o; O1: %o; O2: %o; \n", Evaluate(f_crv,[0,-1,1]), Evaluate(f_crv,[0,1,0]), Evaluate(f_crv,[1,0,0]);
printf "P3: %o; -P3: %o \n", QI_crv!Evaluate(f_crv,[X3,Y3,Z3]),QI_crv!Evaluate(f_crv,[-X3,Y3,Z3]);
printf "\n";

/* Prints for the thesis for this part: 

printf "\n";
printf "Op: %o; O1: %o; O2: %o; \n", Evaluate(f_crv,[0,-1,1]), Evaluate(f_crv,[1,0,0]), Evaluate(f_crv,[0,1,0]);
printf "P3: %o; -P3: %o; \n", QI_crv!Evaluate(f_crv,[X3,Y3,Z3]),QI_crv!Evaluate(f_crv,[-X3,Y3,Z3]);
printf "\n";

printf "Op: %o; O1: %o; O2: %o; \n", Evaluate(f_cnc_PP,[0,-1,1]), Evaluate(f_cnc_PP,[1,0,0]), Evaluate(f_cnc_PP,[0,1,0]);
printf "P1: %o; -P3: %o; \n", QI_crv!Evaluate(f_cnc_PP,[X1,Y1,Z1]), QI_crv!Evaluate(f_cnc_PP,[-X3_PP,Y3_PP,Z3_PP]);
printf "\n";

printf "Op: %o; O1: %o; O2: %o; \n", Evaluate(f_cnc_POp,[0,-1,1]), Evaluate(f_cnc_POp,[1,0,0]), Evaluate(f_cnc_POp,[0,1,0]);
printf "P1: %o; -P3: %o; \n", QI_crv!Evaluate(f_cnc_POp,[X1,Y1,Z1]), Evaluate(f_cnc_POp,[-X3_POp,Y3_POp,Z3_POp]);
printf "P3: %o; \n",QI_crv!Evaluate(f_cnc_POp,[X3_POp,Y3_POp,Z3_POp]);
printf "\n";

printf "Op: %o; O1: %o; O2: %o; \n", Evaluate(f_cnc_PQ,[0,-1,1]), Evaluate(f_cnc_PQ,[1,0,0]), Evaluate(f_cnc_PQ,[0,1,0]);
printf "P1: %o; P2: %o; -P3: %o; \n", Evaluate(f_cnc_PQ,[X1,Y1,Z1]), Evaluate(f_cnc_PQ,[X2,Y2,Z2]), QI_crv!Evaluate(f_cnc_PQ,[-X3,Y3,Z3]);
printf "\n";
*/



//** Numerical part **//

// Example 1 of handpicked curve, for a walkthrough of the process 
BF := Rationals();
R1<a,d> := FunctionField(BF,2);
R2<X,Y,Z> := ProjectiveSpace(R1,2);
Op := [0,-1,1]; O1 := [1,0,0]; O2 := [0,1,0]; 
a := 17; d := 82; lst := [a,d];

P1 := [1,4/9,1]; P2 := [-72/1393, -1361/1231, 1];
X1 := P1[1]; Y1 := P1[2]; Z1 := P1[3];
X2 := P2[1]; Y2 := P2[2]; Z2 := P2[3];

cz := X1*X2*(Y1*Z2-Y2*Z1);
cxy := Z1*Z2*(X1*Z2-X2*Z1 + X1*Y2 - X2*Y1);
cxz := X2*Y2*(Z1^2)-X1*Y1*(Z2^2)+Y1*Y2*(X2*Z1-X1*Z2);

Crv := Curve(R2, (a*(X^2) + (Y^2))*Z^2 - Z^4 - d*(X^2)*(Y^2));
Cnc := Conic(R2,cz*(Z^2 + Y*Z)+cxy*X*Y+cxz*X*Z);

pts := Points(Intersection(Crv,Cnc)); P3_lst := pts diff Seqset([R2!P1,R2!P2,R2!Op,R2!O1,R2!O2]); 
mP3 := P3_lst[1];P3 := [-mP3[1], mP3[2], mP3[3]];

printf "\n ~~~~~~~~~~~~~~~~~~~~~~ \n";
printf " ~~~~~~~~~~~~~~~~~~~~~~ \n";

printf "\n Let's see how the geometric interpretation works for some numerical examples now: \n";

printf "\n Lets consider the elliptic curve with parameters a = %o, d = %o. So the following curve: \n\n \t%o \n",a,d,Crv;

printf "\n We also have two points on the curve P = (%o:%o:%o) and Q = (%o:%o:%o). Let's also consider the points at infinity O1 = (%o:%o:%o), O2 = (%o:%o:%o) and the point Op = (%o:%o:%o). \n",P1[1],P1[2],P1[3],P2[1],P2[2],P2[3],O1[1],O1[2],O1[3],O2[1],O2[2],O2[3],Op[1],Op[2],Op[3];

printf "\n Now we have all we need to build the conic passing through P,Q,Op,O1 and O2: \n\n \t%o \n", Cnc;

printf "\n The points in the intersection of the curve and the conic are the following: \n %o \n", pts;

printf "\n And so we can see that -P3 is (%o:%o:%o). We can do some easy comprobations: is this point in the curve? \n\t%o \n And is it in the conic? \n\t%o \n",mP3[1],mP3[2],mP3[3], mP3 in Crv, mP3 in Cnc;

printf "\n So according to the geometric interpretation, the addition point is P3 = (%o:%o:%o). And indeed it is as we have checked with the addition law. Last checks: is this point in the curve? \n\t%o \n And is it in the conic? \n\t%o \n",P3[1],P3[2],P3[3], P3 in Crv, P3 in Cnc;


// Example 2 of a handpicked curve
BF := Rationals();
a := 82; d := 17; lst := [a,d];
Op := [0,-1,1];O1 := [1,0,0]; O2 := [0,1,0]; 

P1 := [1,9/4,1]; P2 := [72/1393, 1231/1361, 1];
P3 := [1935265/5286047, 16708329/5977204,1]; mP3 := [-1935265/5286047, 16708329/5977204,1];

PP := Ead_AddL_geom(BF,lst,P1,P1);
POp := Ead_AddL_geom(BF,lst,P1,Op);
PQ := Ead_AddL_geom(BF,lst,P1,P2);

printf "\n\n Now we have 2 more examples: \n --> First one with parameters a = %o, d = %o over the rationals. And the two points on the curve: P = (%o:%o:%o) and Q = (%o:%o:%o).\n",a,d,P1[1],P1[2],P1[3],P2[1],P2[2],P2[3];

printf "\n Here is the result of doubling a point:  P + P = (%o:%o:%o). \n", PP[1],PP[2],PP[3];
printf "\n Adding a point to Op:  P + Op = (%o:%o:%o). \n", POp[1],POp[2],POp[3];
printf "\n And adding two different points:  P + Q = (%o:%o:%o). \n", PQ[1],PQ[2],PQ[3];


// Example 3 from savecurves
p := 2^414 - 17;
BF := FiniteField(p);
R1<a,d> := FunctionField(BF,2);
R2<X,Y,Z> := ProjectiveSpace(R1,2);    
a := 1; d := 3617; lst := [a,d];
Op := [0,-1,1];O1 := [1,0,0]; O2 := [0,1,0];
P := [17319886477121189177719202498822615443556957307604340815256226171904769976866975908866528699294134494857887698432266169206165,
34,1];
Crv := Curve(R2, (a*(X^2) + (Y^2))*Z^2 - Z^4 - d*(X^2)*(Y^2));

PP := R2!Ead_AddL_geom(BF,lst,P,P); 
POp := R2!Ead_AddL_geom(BF,lst,P,Op);
PQ := R2!Ead_AddL_geom(BF,lst,P,PP);

printf "\n --> Second one with parameters a = %o, d = %o over GF(2^414 - 17). And the point on the curve P = (%o:%o:%o).\n",a,d,P[1],P[2],P[3];

printf "\n Here is the result of doubling a point: \n \tP + P = (%o:%o:%o). \n \tIs is in the curve? %o \n", PP[1],PP[2],PP[3], PP in Crv;
printf "\n Adding a point to Op: \n \tP + Op = (%o:%o:%o). \n \tIs is in the curve? %o \n", POp[1],POp[2],POp[3], POp in Crv;
printf "\n And adding two different points: \n \tP + 2P = (%o:%o:%o). \n \tIs is in the curve? %o \n", PQ[1],PQ[2],PQ[3], PQ in Crv;


quit;
