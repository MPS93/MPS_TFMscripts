//*** Script to evaluate the complete set of addition laws, studying different cases ***//

printf "\n\n === TRYING THE ADDITION LAW SYSTEM WITH SOME EXAMPLES === \n";


load "../../functions.m";
//load Ead_AddL_eval, ProjCoord, AffCoord


//--> set up to try examples <--//
BF := Rationals();
R1<a,d> := FunctionField(BF,2);
R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8); PR := R2;
P1 := [[X1,Z1],[Y1,T1]]; Q1 := [[X2,Z2],[Y2,T2]];
P2 := [[0,Z1],[Y1,T1]]; Q2 := [[X2,Z2],[Y2,0]];

Ead1 := Z1^2*T1^2 + d*X1^2*Y1^2 - a*X1^2*T1^2 - Y1^2*Z1^2; Ead2 := Z2^2*T2^2 + d*X2^2*Y2^2 - a*X2^2*T2^2 - Y2^2*Z2^2;
I := ideal< R2 | Ead1, Ead2>;
QQ<X1,Y1,Z1,T1,X2,Y2,Z2,T2> := quo<R2 | I>; 
 
//~ evaluating the AddL in general ~//
printf "\n The evaluation function works for general addition: \n";
Ead_AddL_eval(BF,[a,d],P1,Q1);
printf "\n If provided with an exceptional point, it yields the dual addition law: \n";
Ead_AddL_eval(BF,[a,d],P2,Q2);

//~ evaluating the AddL in some numeric cases ~//

// example of handpicked curve
BF := Rationals();
R1<a,d> := FunctionField(BF,2);
R2<x,y> := AffineSpace(R1,2);
a := 82; d := 17;
C := Curve(R2, 1 + d*x^2*y^2 - a*x^2 - y^2);//C;
o := Origin(R2);
p := R2 ! [1,9/4];
q := R2 ! [1,-9/4];

lst := [a,d]; P:= ProjCoord(PR,p); Q := Ead_AddL_eval(BF,lst,P,P); 
printf "\n We have also tried some numerical examples: \n";
printf "\n Given the curve with params a = %o, d = %o and the point P = ((%o : %o),(%o : %o)), we have that 2P = ((%o : %o),(%o : %o)).\n",a,d,P[1][1],P[1][2],P[2][1],P[2][2],Q[1][1],Q[1][2],Q[2][1],Q[2][2];
printf "\n For the same curve we have some other point additions to test its behaviour with exceptional points: \n";

//Exceptional points
BF := RealField();
P := [[0,1],[1,1]]; Q := [[1,Sqrt(17)],[1,0]]; PQ := Ead_AddL_eval(BF,lst,P,Q);
printf "\n [[%.3o,%.3o],[%.3o,%.3o]] + [[%.3o,%.3o],[%.3o,%.3o]] = \n\t[[%.5o,%.5o],[%.5o,%.5o]] \n",P[1][1],P[1][2],P[2][1],P[2][2],Q[1][1],Q[1][2],Q[2][1],Q[2][2],PQ[1][1],PQ[1][2],PQ[2][1],PQ[2][2];
P := [[1,1],[9/4,1]]; Q := [[1,Sqrt(17*82)],[Sqrt(82),Sqrt(17)*9/4]]; PQ := Ead_AddL_eval(BF,lst,P,Q);
printf "\n [[%.3o,%.3o],[%.3o,%.3o]] + [[%.3o,%.3o],[%.3o,%.3o]] = \n\t[[%.3o,%.3o],[%.3o,%.3o]] \n",P[1][1],P[1][2],P[2][1],P[2][2],Q[1][1],Q[1][2],Q[2][1],Q[2][2],PQ[1][1],PQ[1][2],PQ[2][1],PQ[2][2];
P := [[1,1],[9/4,1]]; Q := [[9/4,Sqrt(82)],[-Sqrt(82),1]]; PQ := Ead_AddL_eval(BF,lst,P,Q);
printf "\n [[%.3o,%.3o],[%.3o,%.3o]] + [[%.3o,%.3o],[%.3o,%.3o]] = \n\t[[%.3o,%.3o],[%.3o,%.3o]] \n",P[1][1],P[1][2],P[2][1],P[2][2],Q[1][1],Q[1][2],Q[2][1],Q[2][2],PQ[1][1],PQ[1][2],PQ[2][1],PQ[2][2];
P := [[1,0],[Sqrt(a),Sqrt(d)]]; Q := [[1,0],[Sqrt(a/d),1]]; PQ := Ead_AddL_eval(BF,lst,P,Q);
printf "\n [[%.3o,%.3o],[%.3o,%.3o]] + [[%.3o,%.3o],[%.3o,%.3o]] = \n\t[[%.3o,%.3o],[%.3o,%.3o]] \n",P[1][1],P[1][2],P[2][1],P[2][2],Q[1][1],Q[1][2],Q[2][1],Q[2][2],PQ[1][1],PQ[1][2],PQ[2][1],PQ[2][2];


// example from safecurves

p := 2^414 - 17;
BF := FiniteField(p);
R1<a,d> := FunctionField(BF,2);
R2<x,y> := AffineSpace(R1,2);
PR<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8); 
a := 1; d := 3617;
C := Curve(R2, 1 + d*x^2*y^2 - a*x^2 - y^2);//C;
o := Origin(R2);
p := R2 ! [17319886477121189177719202498822615443556957307604340815256226171904769976866975908866528699294134494857887698432266169206165,
34];

lst := [a,d]; P:= ProjCoord(PR,p); Q := Ead_AddL_eval(BF,lst,P,P); //Q; 

printf "\n Let's try a curve over a finite field F(p), where p = 2^414 - 17, params a = %o, d = %o and the point P = ((%o : %o),(%o : %o)). \n\n Then we have that 2P = ((%o : %o),(%o : %o)).\n Is 2P in the curve? %o",a,d,P[1][1],P[1][2],P[2][1],P[2][2],Q[1][1],Q[1][2],Q[2][1],Q[2][2],AffCoord(BF,Q) in C;

quit;
