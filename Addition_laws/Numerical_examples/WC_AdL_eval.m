//*** Script to evaluate the complete set of addition laws, studying different cases ***//

load "../../functions.m";
//load Ew_AddL_eval


printf "\n\n === TRYING THE ADDITION LAW SYSTEM WITH SOME EXAMPLES === \n";

//--> set up to try examples <--//
BF := Rationals();
R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
R2<X1,Y1,Z1, X2,Y2,Z2> := PolynomialRing(R1,6); 
P1 := [X1,Y1,Z1]; Q1 := [X2,Y2,Z2];
 
//~ evaluating the AddL in general ~//
printf "\n The evaluation function works for general addition: \n";
Ew_AddL_eval(BF,[a1,a2,a3,a4,a6],P1,Q1);
printf "\n If provided with an exceptional point, it yields the dual addition law: \n";
Ew_AddL_eval(BF,[a1,a2,a3,a4,a6],P1,P1);

//~ evaluating the AddL in some numeric cases ~//

// example of handpicked curve
BF := Rationals();
R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
R2<x,y> := AffineSpace(R1,2);
//lst := [0, 0, 1, -7, 6];
lst := [19, 2, -13/7, -5, -6];
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6 := lst[5];
C := Curve(R2, y^2+a1*x*y+a3*y-x^3-a2*x^2-a4*x-a6);//C;
//P1 := [ 175912024457 * 278846, -41450244419357361, 278846^3 ]; P2 := [ -151 * 8, -1845, 8^3 ];
//P3 := [-3, -1, 1]; P4 := [406, 8180, 1];  P5 := [4, -7, 1];
P1 := [2, -253/7, 1]; P2 := [-33/49, 4978/343, 1]; P3 := [-3, 412/7, 1];

PQ := Ew_AddL_eval(BF,lst,P1,P2); 
printf "\n We have also tried some numerical examples: \n";
printf "\n Given the curve with params a1 = %o, a2 = %o, a3 = %o, a4 = %o, a6 = %o and the points P = (%o : %o : %o), Q = (%o : %o : %o),  we have that P + Q = (%o : %o : %o).\n",a1,a2,a3,a4,a6,P1[1],P1[2],P1[3],P2[1],P2[2],P2[3],PQ[1],PQ[2],PQ[3];
printf "\n Is P+Q in the curve? %o \n",R2![R1|PQ[1]/PQ[3],PQ[2]/PQ[3]] in C;
printf "\n For the same curve we have some other point additions to test its behaviour with exceptional points: \n";


//Exceptional points
P := P3; PP := Ew_AddL_eval(BF,lst,P,P);
printf "\n 2P = 2[%o,%o,%o] = \n\t[%o,%o,%o] \n Is 2P in the curve? %o \n",P[1],P[2],P[3],PP[1],PP[2],PP[3],R2![R1|PP[1]/PP[3],PP[2]/PP[3]] in C;
//P := [66/169, -5253/2197,1]; Q := [-475/289, 15471/4913,1]; mQ := [-475/289, -20384/4913,1]; 
P := [-33/49, 4978/343, 1]; Q := [8, 27/7, 1]; mQ := [8, -154,1]; 
PmQ := Ew_AddL_eval(BF,lst,P,mQ); PQ := Ew_AddL_eval(BF,lst,P,Q);
printf "\n P - Q = [%o,%o,%o] - [%o,%o,%o] =  [%o,%o,%o] ==> \n\t P + Q = [%o,%o,%o] \n Is P + Q in the curve? %o \n",P[1],P[2],P[3],Q[1],Q[2],Q[3],PmQ[1]/PmQ[3],PmQ[2]/PmQ[3],1,PQ[1],PQ[2],PQ[3],R2![R1|PQ[1]/PQ[3],PQ[2]/PQ[3]] in C;


// example from safecurves

p := 2^384 - 2^128 - 2^96 + 2^32 - 1;
BF := FiniteField(p);
R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
R2<x,y> := AffineSpace(R1,2);

a := -3; 
b := 27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575;
lst := [a,b];
C := Curve(R2, y^2 - x^3 - a*x - b);//C;
P := [26247035095799689268623156744566981891852923491109213387815615900925518854738050089022388053975719786650872476732087,8325710961489029985546751289520108179287853048861315594709205902480503199884419224438643760392947333078086511627871,1]; 
PP := Ew_AddL_eval(BF,lst,P,P); //Q; 

printf "\n Let's try a curve over a finite field F(p), where p = 2^384 - 2^128 - 2^96 + 2^32 - 1, params a = %o, b = %o and the point P = (%o : %o : %o). \n\n Then we have that 2P = (%o : %o : %o).\n Is 2P in the curve? %o",a,b,P[1],P[2],P[3],PP[1],PP[2],PP[3], R2![R1|PP[1]/PP[3],PP[2]/PP[3]] in C;

quit;
