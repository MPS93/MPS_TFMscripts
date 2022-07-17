//*** script to check facts about the complete set of addition laws for Ew ***//

printf"\n === SYMBOLIC CHECKS ABOUT THE ADDITION LAW SYSTEM === \n";

load "../../functions.m";
//load Ew_AddL

//--> set up to do checks <--//
BF := Rationals();
R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
R2<X1,Y1,Z1, X2,Y2,Z2> := PolynomialRing(R1,6);

Ew1 := Y1^2*Z1+a1*X1*Y1*Z1+a3*Y1*Z1^2-X1^3-a2*X1^2*Z1-a4*X1*Z1^2-a6*Z1^3; 
Ew2 := Y2^2*Z2+a1*X2*Y2*Z2+a3*Y2*Z2^2-X2^3-a2*X2^2*Z2-a4*X2*Z2^2-a6*Z2^3;
I := ideal< R2 | Ew1, Ew2>;
QQ<x1,y1,z1,x2,y2,z2> := quo<R2 | I>;

X3_001 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[1][1]); 
Y3_001 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[1][2]);
Z3_001 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[1][3]); 

X3_010 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[2][1]); 
Y3_010 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[2][2]);
Z3_010 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[2][3]); 


//checking that the result of the AddL is in fact on the curve
printf "\nIs (X3_001 : Y3_001 : Z3_001) on the curve? \n";
if QQ ! Evaluate(Ew1,[X3_001,Y3_001,Z3_001,1,1,1]) eq 0 then printf "   Yes \n"; else printf "   No \n"; end if;
printf "Is (X3_010 : Y3_010 : Z3_010) on the curve? \n";
if QQ ! Evaluate(Ew1,[X3_010,Y3_010,Z3_010,1,1,1]) eq 0 then printf "   Yes \n"; else printf "   No \n"; end if;


//checking that the system is well defined:
f_XZ := X3_001*Z3_010; g_XZ := Z3_001*X3_010;
printf "\nDoes it hold that X3*Z'3 = X'3*Z3 \n"; printf "   "; (QQ ! f_XZ) eq (QQ ! g_XZ);

f_YZ := Y3_001* Z3_010; g_YZ := Z3_001*Y3_010; 
printf "Does it hold that Y3*Z'3 = Y'3*Z3 \n"; printf "   "; (QQ ! f_YZ) eq (QQ ! g_YZ);

printf "\nWe will not be proving completeness in all cases, but at least we will try to give a taste of it.\n";
printf "\nFor starters, we can check if P1 = P2 is indeed exceptional for (X3_001 : Y3_001 : Z3_001), i.e., do the three formulae evaluate to 0 in this case? \n %o \n", 
QQ!Evaluate(X3_001,[X1,Y1,Z1,X1,Y1,Z1]) eq 0 and QQ!Evaluate(Y3_001,[X1,Y1,Z1,X1,Y1,Z1]) eq 0 and QQ!Evaluate(Z3_001,[X1,Y1,Z1,X1,Y1,Z1]) eq 0;

printf "\nThings get complicated when trying to check the fact that P1, P2 are exceptional for (X3_010 : Y3_010 : Z3_010) when Y(P1-P2) = 0.";
printf " So we can at least take a look if P1 = P2 is exceptional for this second addition law.";

a1 := 0; a2 := 0; a3 := 0; a4 := 0; 
X3_010 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[2][1]); 
Y3_010 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[2][2]);
Z3_010 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[2][3]); 
X3_PP := Evaluate(X3_010,[X1,Y1,1,X1,Y1,1]);
Y3_PP := Evaluate(Y3_010,[X1,Y1,1,X1,Y1,1]);
Z3_PP := Evaluate(Z3_010,[X1,Y1,1,X1,Y1,1]);


printf " That happens to be also rather complicated for the general case so we will take a look at the case a1 = a2 = a3 = a4 = 0 and Z1 = 1 over the rationals. That is, for the Weierstrass equation Y^2 = X^3 + a6.\n";
printf "\nIn said case, (X3_010 : Y3_010 : Z3_010) is the following: \n X3_010 = %o \n Y3_010 = %o \n Z3_010 = %o \n", X3_PP, Y3_PP, Z3_PP;
printf "\nX3_010 evaluates to 0 in three cases: if X1 = 0, Y1^2 = 9a6,  or Y1 = 0.  Thus let's see how in none of these cases it can happen that Y3_010 = Z3_010 = 0.\n";
printf "\nIf X1=0, then Z3_010 = %o. So for Z3_010 = 0, it needs to hold that Y1^2 = -3a6. But looking at the curve equation, when X1 = 0 then Y1^2 = a6 so if both hold, then -3 = 1, which is a contradiction. \n", Evaluate(Z3_PP,[0,Y1,1,0,Y1,1]);
printf "\nIf Y1^2 = 9a6, then by looking at Z3_010 again it needs to hold that X1^3 = -4a6. But looking at the curve equation, we have that X1^3 = 8a6, which is again a contradiction.\n";
printf "\nLastly, if Y1=0 then Y3_010 = %o. If we want that to be 0, then we need X1^3 = a6/2, but one again that contradicts the curve equation since in this case it is X1^3 = - a6.\n", Evaluate(Y3_PP,[X1,0,1,X1,0,1]);

printf "\nTherefore in this case we can see how P1=P2 is not exceptional for (X3_010 : Y3_010 : Z3_010), so the system is complete.\n";

//Factorization(Evaluate(X3_010,[X1,Y1,1,X1,Y1,1]));
//Factorization(Evaluate(Y3_010,[X1,Y1,1,X1,Y1,1]));
//Factorization(Evaluate(Z3_010,[X1,Y1,1,X1,Y1,1]));


//--> This doesn't help, waaay too long output
//GB := GroebnerBasis(ideal<R2|Ew1,Ew2,X3_PP,Y3_PP,Z3_PP>);

// Numeric example of exceptional pair for 010:
/*
BF := Rationals();
R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
R2<X1,Y1,Z1, X2,Y2,Z2> := PolynomialRing(R1,6);

lst := [0, 0, 1, -7, 6];
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6 := lst[5];

Ew1 := Y1^2*Z1+a1*X1*Y1*Z1+a3*Y1*Z1^2-X1^3-a2*X1^2*Z1-a4*X1*Z1^2-a6*Z1^3; 
Ew2 := Y2^2*Z2+a1*X2*Y2*Z2+a3*Y2*Z2^2-X2^3-a2*X2^2*Z2-a4*X2*Z2^2-a6*Z2^3;
I := ideal< R2 | Ew1, Ew2>;
QQ<x1,y1,z1,x2,y2,z2> := quo<R2 | I>;

X3_001 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[1][1]); 
Y3_001 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[1][2]);
Z3_001 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[1][3]); 

X3_010 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[2][1]); 
Y3_010 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[2][2]);
Z3_010 := (R2 ! Ew_AddL(BF,[a1,a2,a3,a4,a6])[2][3]); 

P := [66/169, -5253/2197,1]; Q := [-475/289, 15471/4913,1];

Evaluate(X3_001, [P[1],P[2],P[3],Q[1],Q[2],Q[3]]);
Evaluate(Y3_001, [P[1],P[2],P[3],Q[1],Q[2],Q[3]]);
Evaluate(Z3_001, [P[1],P[2],P[3],Q[1],Q[2],Q[3]]);

Evaluate(X3_010, [P[1],P[2],P[3],Q[1],Q[2],Q[3]]);
Evaluate(Y3_010, [P[1],P[2],P[3],Q[1],Q[2],Q[3]]);
Evaluate(Z3_010, [P[1],P[2],P[3],Q[1],Q[2],Q[3]]);
*/

quit;
