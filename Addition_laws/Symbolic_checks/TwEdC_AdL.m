//*** script to check facts about the complete set of addition laws for EadC ***//

load "../../functions.m";
//load Ead_AddL

printf"\n === SYMBOLIC CHECKS ABOUT THE ADDITION LAW SYSTEM === \n";

//--> set up to do checks <--//
BF := Rationals();
R1<a,d> := FunctionField(BF,2);
R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8);
P := [[X1,Z1],[Y1,T1]]; Q := [[X2,Z2],[Y2,T2]];

Ead1 := Z1^2*T1^2 + d*X1^2*Y1^2 - a*X1^2*T1^2 - Y1^2*Z1^2; 
Ead2 := Z2^2*T2^2 + d*X2^2*Y2^2 - a*X2^2*T2^2 - Y2^2*Z2^2;
I := ideal< R2 | Ead1, Ead2>;
QQ<X1,Y1,Z1,T1,X2,Y2,Z2,T2> := quo<R2 | I>;

X3_o := (R2 ! Ead_AddL(BF,[a,d])[1][1][1]); Y3_o := (R2 ! Ead_AddL(BF,[a,d])[1][2][1]);
Z3_o := (R2 ! Ead_AddL(BF,[a,d])[1][1][2]); T3_o := (R2 ! Ead_AddL(BF,[a,d])[1][2][2]);

X3_d := (R2 ! Ead_AddL(BF,[a,d])[2][1][1]); Y3_d := (R2 ! Ead_AddL(BF,[a,d])[2][2][1]);
Z3_d := (R2 ! Ead_AddL(BF,[a,d])[2][1][2]); T3_d := (R2 ! Ead_AddL(BF,[a,d])[2][2][2]);

//checking that the system is indeed well defined:
f_XZ := X3_o*Z3_d; g_XZ := Z3_o*X3_d;
printf "\nDoes it hold that X3*Z'3 = X'3*Z3 \n"; printf "   "; (QQ ! f_XZ) eq (QQ ! g_XZ);

f_YT := Y3_o* T3_d; g_YT := T3_o*Y3_d; 
printf "Does it hold that Y3*T'3 = Y'3*T3 \n"; printf "   "; (QQ ! f_YT) eq (QQ ! g_YT);


//checking that the result of the AddL is in fact on the curve
printf "\nIs ((X3_o : Z3_o),(Y3_o : T3_o)) on the curve? \n";
if QQ ! Evaluate(Ead1,[X3_o,Y3_o,Z3_o,T3_o,1,1,1,1]) eq 0 then printf "   Yes \n"; else printf "   No \n"; end if;
printf "Is ((X3_d : Z3_d),(Y3_d : T3_d)) on the curve? \n";
if QQ ! Evaluate(Ead1,[X3_d,Y3_d,Z3_d,T3_d,1,1,1,1]) eq 0 then printf "   Yes \n"; else printf "   No \n"; end if;


//checking the exceptional points is done in TwEdC_AdL_completeness.m

printf "\nLast check left is that the system is indeed complete. This is done in the file TwEdC_AdL_completeness.m \n";

quit;
