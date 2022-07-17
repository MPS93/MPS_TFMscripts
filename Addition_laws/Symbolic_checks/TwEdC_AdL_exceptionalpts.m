//*** Script to solve the system of equations defined by the Groebner basis of each addition law coordinate ***//
//** This way we find all the exceptional points of the addition law **//

load "../../functions.m";
//load Ead_AddL, Ead_AddL_evals, ProjCoord, AffCoord

printf "\n === STUDYING THE EXCEPTIONAL POINTS OF THE SYSTEM OF ADD. LAWS === \n";

//** Symbolic part **//

// Set up to get the formulae to work with them:
BF := Rationals();
R1<a,d,sa,sd> := FunctionField(BF,4);
R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8); 
a := sa^2; d := sd^2;
P := [[X1,Z1],[Y1,T1]]; Q := [[X2,Z2],[Y2,T2]];

Ead1 := Z1^2*T1^2 + d*X1^2*Y1^2 - a*X1^2*T1^2 - Y1^2*Z1^2; 
Ead2 := Z2^2*T2^2 + d*X2^2*Y2^2 - a*X2^2*T2^2 - Y2^2*Z2^2;

X3_o := (R2 ! Ead_AddL(BF,[a,d] : except := true)[1][1][1]); Y3_o := (R2 ! Ead_AddL(BF,[a,d] : except := true)[1][2][1]);
Z3_o := (R2 ! Ead_AddL(BF,[a,d] : except := true)[1][1][2]); T3_o := (R2 ! Ead_AddL(BF,[a,d] : except := true)[1][2][2]);

X3_d := (R2 ! Ead_AddL(BF,[a,d] : except := true)[2][1][1]); Y3_d := (R2 ! Ead_AddL(BF,[a,d] : except := true)[2][2][1]);
Z3_d := (R2 ! Ead_AddL(BF,[a,d] : except := true)[2][1][2]); T3_d := (R2 ! Ead_AddL(BF,[a,d] : except := true)[2][2][2]);

//In order to simplify the systems, we do some variable changes.
RR2<x1,y1,z1,t1, x2,y2,z2,t2> := FunctionField(R1,8); 
RR3<m,n,x,y> := FunctionField(R1,4);
h := hom<RR3->RR2|x1/z1,y1/t1,x2/z2,y2/t2>;
hinv := hom<RR2->RR3| Numerator(m),Numerator(n),Denominator(m),Denominator(n),Numerator(x),Numerator(y),Denominator(x),Denominator(y)>;

hX3_o := hinv(RR2!X3_o);
hZ3_o := hinv(RR2!Z3_o);

hY3_o := hinv(RR2!Y3_o);
hT3_o := hinv(RR2!T3_o);

hX3_d := hinv(RR2!X3_d);
hZ3_d := hinv(RR2!Z3_d);

hY3_d := hinv(RR2!Y3_d);
hT3_d := hinv(RR2!T3_d);


// By hand or with Matlab we solve the systems and get the following results:
solx_o1 := [-1/(n*sd),1/(n*sd)]; soly_o1 := [1/(m*sd),-1/(m*sd)];
solx_o2 := [1/(m*sa*sd),-1/(m*sa*sd)]; soly_o2 := [sa/(n*sd),-sa/(n*sd)];
solx_d1 := [-n/sa,n/sa]; soly_d1 := [m*sa,-m*sa]; 
solx_d2 := [-m,m]; soly_d2 := [-n,n];

// Then we undo the variable change to obtain the solutions:

// here we have done some computations with the coefficients in order to obtain the exceptional points as we have written them in the theory. 
// But they are not actually necessary for the computations.

solxz_o1 := [h(s) : s in solx_o1]; solx2_o1 := [Numerator(s) : s in solxz_o1]; solz2_o1 := [Denominator(s) : s in solxz_o1];
solyt_o1 := [h(s) : s in soly_o1]; soly2_o1 := [Numerator(s) : s in solyt_o1]; solt2_o1 := [Denominator(s) : s in solyt_o1];
solP2_o1 := [[[R2|solx2_o1[i]/Coefficients(solx2_o1[i])[1],solz2_o1[i]/Coefficients(solx2_o1[i])[1]],[R2|soly2_o1[i]/Coefficients(soly2_o1[i])[1],solt2_o1[i]/Coefficients(soly2_o1[i])[1]]] : i in [1..#solx2_o1]]; 

solxz_o2 := [h(s) : s in solx_o2]; solx2_o2 := [Numerator(s) : s in solxz_o2]; solz2_o2 := [Denominator(s) : s in solxz_o2];
solyt_o2 := [h(s) : s in soly_o2]; soly2_o2 := [Numerator(s) : s in solyt_o2]; solt2_o2 := [Denominator(s) : s in solyt_o2];
solP2_o2 := [[[R2|solx2_o2[i]*Numerator(1/Coefficients(solx2_o2[i])[1]),solz2_o2[i]*Numerator(1/Coefficients(solx2_o2[i])[1])],[R2|soly2_o2[i]*Denominator(Coefficients(soly2_o2[i])[1]),solt2_o2[i]*Denominator(Coefficients(soly2_o2[i])[1])]] : i in [1..#solx2_o2]]; 

solxz_d1 := [h(s) : s in solx_d1]; solx2_d1 := [Numerator(s) : s in solxz_d1]; solz2_d1 := [Denominator(s) : s in solxz_d1];
solyt_d1 := [h(s) : s in soly_d1]; soly2_d1 := [Numerator(s) : s in solyt_d1]; solt2_d1 := [Denominator(s) : s in solyt_d1];
solP2_d1 := [[[R2|solx2_d1[i]/Coefficients(solx2_d1[i])[1],solz2_d1[i]/Coefficients(solx2_d1[i])[1]],[R2|soly2_d1[i],solt2_d1[i]]] : i in [1..#solx2_d1]]; 

solxz_d2 := [h(s) : s in solx_d2]; solx2_d2 := [Numerator(s) : s in solxz_d2]; solz2_d2 := [Denominator(s) : s in solxz_d2];
solyt_d2 := [h(s) : s in soly_d2]; soly2_d2 := [Numerator(s) : s in solyt_d2]; solt2_d2 := [Denominator(s) : s in solyt_d2];
solP2_d2 := [[[R2|solx2_d2[i],solz2_d2[i]],[R2|soly2_d2[i],solt2_d2[i]]] : i in [1..#solx2_d2]]; 

P2_o1 := solP2_o1[2]; mP2_o1 := solP2_o1[1]; 
P2_o2 := solP2_o2[1]; mP2_o2 := solP2_o2[2]; 
P2_d1 := solP2_d1[2]; mP2_d1 := solP2_d1[1];
P2_d2 := solP2_d2[2]; mP2_d2 := solP2_d2[1]; 

lst_excp := [P2_o1,mP2_o1,P2_o2,mP2_o2,P2_d1,mP2_d1,P2_d2,mP2_d2];


// Exceptional pairs to print
prP1 := "((X1 : Z1), (Y1 : T1))";
prP2_o1 := "((T1 : Sqrt(d)*Y1),(Z1 : -Sqrt(d)*X1))"; prmP2_o1 := "((T1 : -Sqrt(d)*Y1),(Z1 : Sqrt(d)*X1))";
prP2_o2 := "((Z1 : Sqrt(a*d)*X1),(Sqrt(a)*T1 : Sqrt(d)*Y1))"; prmP2_o2 := "((Z1 : -Sqrt(a*d)*X1),(-Sqrt(a)*T1 : Sqrt(d)*Y1))";
prP2_d1 := "((Y1 : Sqrt(a)*T1),(-Sqrt(a)*X1 : Z1))"; prmP2_d1 := "((Y1 : -Sqrt(a)*T1),(Sqrt(a)*X1 : Z1))";
prP2_d2 := "((X1 : Z1),(Y1 : T1))"; prmP2_d2 := "((-X1 : Z1),(-Y1 : T1))"; 


// Checking that tbe solutions are indeed exceptional points
I := ideal<R2|Ead1,Ead2>; QI := R2/I;

o1_excp := [[QI!Evaluate(X3_o,[X1,Y1,Z1,T1,pt[1][1],pt[2][1],pt[1][2],pt[2][2]]),QI!Evaluate(Z3_o,[X1,Y1,Z1,T1,pt[1][1],pt[2][1],pt[1][2],pt[2][2]])] : pt in lst_excp];
o2_excp := [[QI!Evaluate(Y3_o,[X1,Y1,Z1,T1,pt[1][1],pt[2][1],pt[1][2],pt[2][2]]),QI!Evaluate(T3_o,[X1,Y1,Z1,T1,pt[1][1],pt[2][1],pt[1][2],pt[2][2]])] : pt in lst_excp];
d1_excp := [[QI!Evaluate(X3_d,[X1,Y1,Z1,T1,pt[1][1],pt[2][1],pt[1][2],pt[2][2]]),QI!Evaluate(Z3_d,[X1,Y1,Z1,T1,pt[1][1],pt[2][1],pt[1][2],pt[2][2]])] : pt in lst_excp];
d2_excp := [[QI!Evaluate(Y3_d,[X1,Y1,Z1,T1,pt[1][1],pt[2][1],pt[1][2],pt[2][2]]),QI!Evaluate(T3_d,[X1,Y1,Z1,T1,pt[1][1],pt[2][1],pt[1][2],pt[2][2]])] : pt in lst_excp];

chck_o1_0 := &+Flat(o1_excp[1] cat o1_excp[2]) eq 0; o1_excp := Remove(Remove(o1_excp,1),1); chck_o1_1 := &*Flat(o1_excp) eq 0; //do this modulo the curve?
chck_o2_0 := &+Flat(o2_excp[3] cat o2_excp[4]) eq 0; o2_excp := Remove(Remove(o2_excp,3),3); chck_o2_1 := &*Flat(o2_excp) eq 0;
chck_d1_0 := &+Flat(d1_excp[5] cat d1_excp[6]) eq 0; d1_excp := Remove(Remove(d1_excp,5),5); chck_d1_1 := &*Flat(d1_excp) eq 0;
chck_d2_0 := &+Flat(d2_excp[7] cat d2_excp[8]) eq 0; d2_excp := Remove(Remove(d2_excp,7),7); chck_d2_1 := &*Flat(d2_excp) eq 0;

printf "\n Let's start doing some symbolic computations. \n";
printf "\n To find the exceptional points, we are going to solve each of these systems for P2, assuming we know P1 = %o: \n",prP1;
printf "\n --> For (X3:Z3) of the original law: \n\t%o \t= %o \t= 0 \n\t%o \t= %o \t= 0 \n",X3_o,hX3_o,Z3_o,hZ3_o;
printf "\n \t -> Solutions: \tP2 = %o \n \t\t\tP2 = %o \n",prP2_o1,prmP2_o1;
printf "\n \t -> Are these solutions points on the curve? %o \n",(QI!Evaluate(Ead1,[P2_o1[1][1],P2_o1[2][1],P2_o1[1][2],P2_o1[2][2],1,1,1,1]) eq 0) and (QI!Evaluate(Ead1,[mP2_o1[1][1],mP2_o1[2][1],mP2_o1[1][2],mP2_o1[2][2],1,1,1,1]) eq 0);
printf "\n \t -> Does it indeed hold that X3 = 0, Z3 = 0 for these solutions? %o \n",chck_o1_0;
printf "\n --> For (Y3:T3) of the original law: \n\t%o \t= %o \t= 0 \n\t%o \t= %o \t= 0 \n",Y3_o,hY3_o,T3_o,hT3_o;
printf "\n \t -> Solutions: \tP2 = %o \n \t\t\tP2 = %o \n",prP2_o2,prmP2_o2;
printf "\n \t -> Are these solutions points on the curve? %o \n",(QI!Evaluate(Ead1,[P2_o2[1][1],P2_o2[2][1],P2_o2[1][2],P2_o2[2][2],1,1,1,1]) eq 0) and (QI!Evaluate(Ead1,[mP2_o2[1][1],mP2_o2[2][1],mP2_o2[1][2],mP2_o2[2][2],1,1,1,1]) eq 0);
printf "\n \t -> Does it indeed hold that Y3 = 0, T3 = 0 for these solutions? %o \n",chck_o2_0;
printf "\n --> For (X3:Z3) of the dual law: \n\t%o \t= %o \t= 0 \n\t%o \t= %o \t= 0 \n",X3_d,hX3_d,Z3_d,hZ3_d;
printf "\n \t -> Solutions: \tP2 = %o \n \t\t\tP2 = %o \n",prP2_d1,prmP2_d1;
printf "\n \t -> Are these solutions points on the curve? %o \n",(QI!Evaluate(Ead1,[P2_d1[1][1],P2_d1[2][1],P2_d1[1][2],P2_d1[2][2],1,1,1,1]) eq 0) and (QI!Evaluate(Ead1,[mP2_d1[1][1],mP2_d1[2][1],mP2_d1[1][2],mP2_d1[2][2],1,1,1,1]) eq 0);
printf "\n \t -> Does it indeed hold that X3 = 0, Z3 = 0 for these solutions? %o \n",chck_d1_0;
printf "\n --> For (Y3:T3) of the dual law: \n\t%o \t= %o \t= 0 \n\t%o \t= %o \t= 0 \n",Y3_d,hY3_d,T3_d,hT3_d;
printf "\n \t -> Solutions: \tP2 = %o \n \t\t\tP2 = %o \n",prP2_d2,prmP2_d2;
printf "\n \t -> Are these solutions points on the curve? %o \n",(QI!Evaluate(Ead1,[P2_d2[1][1],P2_d2[2][1],P2_d2[1][2],P2_d2[2][2],1,1,1,1]) eq 0) and (QI!Evaluate(Ead1,[mP2_d2[1][1],mP2_d2[2][1],mP2_d2[1][2],mP2_d2[2][2],1,1,1,1]) eq 0);
printf "\n \t -> Does it indeed hold that Y3 = 0, T3 = 0 for these solutions? %o \n",chck_d2_0;

printf "\n One last check: do any of the pairs X3,Z3 or Y3,T3 evaluate to zero in any exceptional pair aside from their own? %o \n", chck_o1_1 or chck_o2_1 or chck_d1_1 or chck_d2_1;


//** Numerical part **//

printf "\n ~~~~~~~~~~~~~~~~~~~~~~ \n";
printf " ~~~~~~~~~~~~~~~~~~~~~~ \n";

printf "\n Now it's time to try some numerical examples. \n";


// Example to illustrate the completeness of the original law for d non-square 
p := 2^222 - 117;
BF := FiniteField(p);
R1<a,d> := FunctionField(BF,2);
R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8);
a := 1; d := 160102; lst := [a,d];
p := [2705691079882681090389589001251962954446177367541711474502428610129, 28];
P1 := ProjCoord(R2,p); 
Ead1 := Z1^2*T1^2 + d*X1^2*Y1^2 - a*X1^2*T1^2 - Y1^2*Z1^2; 

// Checking if a and/or d are squares in BF
tf_d := IsSquare(d);
tf,sa := IsSquare(a); tf_a := IsCoercible(BF,sa);

// Computing the exceptional points
P2_d1 := [[Y1,sa*T1],[-sa*X1,Z1]]; mP2_d1 := [[Y1,-sa*T1],[sa*X1,Z1]];
P2_d2 := [[X1,Z1],[Y1,T1]]; mP2_d2 := [[-X1,Z1],[-Y1,T1]];

P2_d1_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in P2_d1];
mP2_d1_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in mP2_d1];

P2_d2_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in P2_d2];
mP2_d2_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in mP2_d2];

// Checking that they are on the curve
chck_P2d1 := Evaluate(Ead1,[P2_d1_num[1][1], P2_d1_num[2][1],P2_d1_num[1][2], P2_d1_num[2][2],1,1,1,1]) eq 0;
chck_mP2d1 := Evaluate(Ead1,[mP2_d1_num[1][1], mP2_d1_num[2][1],mP2_d1_num[1][2], mP2_d1_num[2][2],1,1,1,1]) eq 0;
chck_P2d2 := Evaluate(Ead1,[P2_d2_num[1][1], P2_d2_num[2][1],P2_d2_num[1][2], P2_d2_num[2][2],1,1,1,1]) eq 0;
chck_mP2d2 := Evaluate(Ead1,[mP2_d2_num[1][1], mP2_d2_num[2][1],mP2_d2_num[1][2], mP2_d2_num[2][2],1,1,1,1]) eq 0;

// Checking that they are indeed exceptional
d1 := Ead_AddL_evals(BF,lst,P1,P2_d1_num)[2][1];md1 := Ead_AddL_evals(BF,lst,P1,mP2_d1_num)[2][1];
d2 :=  Ead_AddL_evals(BF,lst,P1,P2_d2_num)[2][2];md2 := Ead_AddL_evals(BF,lst,P1,mP2_d2_num)[2][2];

// Printing it all
printf "\n For starters, we have an example of the Edwards curve over the finite field K = GF(2^222 - 117) defined by a = %o and d = %o. \n", a,d;
printf "\n Note that: Is a square in K? %o.  And is d square in K? %o. \n", tf_a,tf_d;
printf "\n So according to the theory, the original addition law should be complete in this case. And if we look at its exceptional points, we can see that they all include a square root of d so this curve indeed doesn't have exceptional points over K for the original addition law. \n";
printf "\n And about the dual addition law: let's consider the point on the curve P1 = ((%o : %o),(%o : %o)).\n",P1[1][1],P1[1][2],P1[2][1],P1[2][2];
printf "\n   So for X3,Z3 the pair P1, P2 will be exceptional iff P2 =  ((%o : %o),(%o : %o)) or P2 = ((%o : %o),(%o : %o)). Does this hold? %o \n", P2_d1_num[1][1],P2_d1_num[1][2],P2_d1_num[2][1],P2_d1_num[2][2],mP2_d1_num[1][1],mP2_d1_num[1][2],mP2_d1_num[2][1],mP2_d1_num[2][2],(&+d1 eq 0) and (&+md1 eq 0); 
printf "     Are both versions of P2 on the curve? %o \n", chck_P2d1 and chck_mP2d1;
printf "\n   Similarly for Y3,T3 the pair P1, P2 will be exceptional iff P2 =  ((%o : %o),(%o : %o)) or P2 = ((%o : %o),(%o : %o)). Does this hold? %o \n",P2_d2_num[1][1],P2_d2_num[1][2],P2_d2_num[2][1],P2_d2_num[2][2],mP2_d2_num[1][1],mP2_d2_num[1][2],mP2_d2_num[2][1],mP2_d2_num[2][2],(&+d2 eq 0) and (&+md2 eq 0);
printf "     Are both versions of P2 on the curve? %o \n", chck_P2d2 and chck_mP2d2;

//DO THIS TO TALK ABOUT PRECISION
BF := RealField(); 
R1<a,d> := PolynomialRing(BF,2);
R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8);
a := 82; d := 17; lst := [a,d];
p := [1,9/4]; 
P1 := ProjCoord(R2,p); 
Ead1 := Z1^2*T1^2 + d*X1^2*Y1^2 - a*X1^2*T1^2 - Y1^2*Z1^2;
sa := Sqrt(a); sd := Sqrt(d);

// Computing the exceptional points
P2_o1 := [[T1,sd*Y1],[Z1,-sd*X1]]; mP2_o1 := [[T1,-sd*Y1],[Z1,sd*X1]];
P2_o2 := [[Z1,sa*sd*X1],[sa*T1,sd*Y1]]; mP2_o2 := [[Z1,-sa*sd*X1],[-sa*T1,sd*Y1]];
P2_d1 := [[Y1,sa*T1],[-sa*X1,Z1]]; mP2_d1 := [[Y1,-sa*T1],[sa*X1,Z1]];
P2_d2 := [[X1,Z1],[Y1,T1]]; mP2_d2 := [[-X1,Z1],[-Y1,T1]];

P2_o1_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in P2_o1];
mP2_o1_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in mP2_o1];
o1 := Ead_AddL_evals(BF,lst,P1,P2_o1_num)[2][1];mo1 := Ead_AddL_evals(BF,lst,P1,mP2_o1_num)[2][1];

P2_o2_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in P2_o2];
mP2_o2_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in mP2_o2];
o2 := Ead_AddL_evals(BF,lst,P1,P2_o2_num)[2][1];mo2 := Ead_AddL_evals(BF,lst,P1,mP2_o2_num)[2][1];

P2_d1_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in P2_d1];
mP2_d1_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in mP2_d1];
d1 := Ead_AddL_evals(BF,lst,P1,P2_d1_num)[2][1];md1 := Ead_AddL_evals(BF,lst,P1,mP2_d1_num)[2][1];

P2_d2_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in P2_d2];
mP2_d2_num := [[Evaluate(c,[P1[1][1],P1[2][1],P1[1][2],P1[2][2],1,1,1,1]) : c in pair] : pair in mP2_d2];
d2 :=  Ead_AddL_evals(BF,lst,P1,P2_d2_num)[2][2];md2 := Ead_AddL_evals(BF,lst,P1,mP2_d2_num)[2][2];


// Checking that the exceptional points are actually on the curve
chck_P2o1 := Evaluate(Ead1,[P2_o1_num[1][1], P2_o1_num[2][1],P2_o1_num[1][2], P2_o1_num[2][2],1,1,1,1]) lt 1E-10;
chck_mP2o1 := Evaluate(Ead1,[mP2_o1_num[1][1], mP2_o1_num[2][1],mP2_o1_num[1][2], mP2_o1_num[2][2],1,1,1,1]) lt 1E-10;
chck_P2o2 := Evaluate(Ead1,[P2_o2_num[1][1], P2_o2_num[2][1],P2_o2_num[1][2], P2_o2_num[2][2],1,1,1,1]) lt 1E-10;
chck_mP2o2 := Evaluate(Ead1,[mP2_o2_num[1][1], mP2_o2_num[2][1],mP2_o2_num[1][2], mP2_o2_num[2][2],1,1,1,1]) lt 1E-10;
chck_P2d1 := Evaluate(Ead1,[P2_d1_num[1][1], P2_d1_num[2][1],P2_d1_num[1][2], P2_d1_num[2][2],1,1,1,1]) lt 1E-10;
chck_mP2d1 := Evaluate(Ead1,[mP2_d1_num[1][1], mP2_d1_num[2][1],mP2_d1_num[1][2], mP2_d1_num[2][2],1,1,1,1]) lt 1E-10;
chck_P2d2 := Evaluate(Ead1,[P2_d2_num[1][1], P2_d2_num[2][1],P2_d2_num[1][2], P2_d2_num[2][2],1,1,1,1]) lt 1E-10;
chck_mP2d2 := Evaluate(Ead1,[mP2_d2_num[1][1], mP2_d2_num[2][1],mP2_d2_num[1][2], mP2_d2_num[2][2],1,1,1,1]) lt 1E-10;

// Checking that the exceptional points are indeed exceptional
o1 := Ead_AddL_evals(BF,lst,P1,P2_o1_num)[1][1]; mo1 := Ead_AddL_evals(BF,lst,P1,mP2_o1_num)[1][1];
o2 := Ead_AddL_evals(BF,lst,P1,P2_o2_num)[1][2]; mo2 := Ead_AddL_evals(BF,lst,P1,mP2_o2_num)[1][2];
d1 := Ead_AddL_evals(BF,lst,P1,P2_d1_num)[2][1]; md1 := Ead_AddL_evals(BF,lst,P1,mP2_d1_num)[2][1];
d2 := Ead_AddL_evals(BF,lst,P1,P2_d2_num)[2][2]; md2 := Ead_AddL_evals(BF,lst,P1,mP2_d2_num)[2][2];

// Printing it all
printf "\n ~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n We can also try an example for a curve over the reals. Let's consider the curve over K = R with parameters a = %o, d = %o. Since we are working over the reals and both a and d are positive integers, they are squares in K. So both addition laws have exceptional points for this curve over K.\n",a,d;
printf "\n The point P1 = ((%o:%o),(%o:%o)) is a point on this curve and with it, we can compute all the exceptional P2's of both addition laws: \n",P1[1][1],P1[1][2],P1[2][1],P1[2][2];

printf "\n --> For (X3:Z3) of the original law:";
printf "\n \t -> P1, Q are exceptional iff Q = P2 = ((%o:%o),(%o:%o)) or Q = mP2 = ((%o:%o),(%o:%o))",P2_o1_num[1][1],P2_o1_num[1][2],P2_o1_num[2][1],P2_o1_num[2][2],mP2_o1_num[1][1],mP2_o1_num[1][2],mP2_o1_num[2][1],mP2_o1_num[2][2];
printf "\n \t -> Are P2 and mP2 on the curve? %o",chck_P2o1 and chck_mP2o1;
printf "\n \t -> Computing (X3:Z3) for P1,P2 yields : \n\t\t(%o:%o)",o1[1],o1[2];
printf "\n \t -> Computing (X3:Z3) for P1,mP2 yields : \n\t\t(%o:%o) \n",mo1[1],mo1[2];
printf "\n --> For (Y3:T3) of the original law:";
printf "\n \t -> P1, Q are exceptional iff Q = P2 = ((%o:%o),(%o:%o)) or Q = mP2 = ((%o:%o),(%o:%o))",P2_o2_num[1][1],P2_o2_num[1][2],P2_o2_num[2][1],P2_o2_num[2][2],mP2_o2_num[1][1],mP2_o2_num[1][2],mP2_o2_num[2][1],mP2_o2_num[2][2];
printf "\n \t -> Are P2 and mP2 on the curve? %o",chck_P2o2 and chck_mP2o2;
printf "\n \t -> Computing (Y3:T3) for P1,P2 yields : \n\t\t(%o:%o)",o2[1],o2[2];
printf "\n \t -> Computing (Y3:T3) for P1,mP2 yields : \n\t\t(%o:%o) \n",mo2[1],mo2[2];
printf "\n --> For (X3:Z3) of the dual law:";
printf "\n \t -> P1, Q are exceptional iff Q = P2 = ((%o:%o),(%o:%o)) or Q = mP2 = ((%o:%o),(%o:%o))",P2_d1_num[1][1],P2_d1_num[1][2],P2_d1_num[2][1],P2_d1_num[2][2],mP2_d1_num[1][1],mP2_d1_num[1][2],mP2_d1_num[2][1],mP2_d1_num[2][2];
printf "\n \t -> Are P2 and mP2 on the curve? %o",chck_P2d1 and chck_mP2d1;
printf "\n \t -> Computing (X3:Z3) for P1,P2 yields : \n\t\t(%o:%o)",d1[1],d1[2];
printf "\n \t -> Computing (X3:Z3) for P1,mP2 yields : \n\t\t(%o:%o) \n",md1[1],md1[2];
printf "\n --> For (Y3:T3) of the dual law:";
printf "\n \t -> P1, Q are exceptional iff Q = P2 = ((%o:%o),(%o:%o)) or Q = mP2 = ((%o:%o),(%o:%o))",P2_d2_num[1][1],P2_d2_num[1][2],P2_d2_num[2][1],P2_d2_num[2][2],mP2_d2_num[1][1],mP2_d2_num[1][2],mP2_d2_num[2][1],mP2_d2_num[2][2];
printf "\n \t -> Are P2 and mP2 on the curve? %o",chck_P2d2 and chck_mP2d2;
printf "\n \t -> Computing (Y3:T3) for P1,P2 yields : \n\t\t(%o:%o)",d2[1],d2[2];
printf "\n \t -> Computing (Y3:T3) for P1,mP2 yields : \n\t\t(%o:%o) \n",md2[1],md2[2];

printf "\n Note that we don't get always 0 but sometimes some very small number. That is due to the fact that we are working with reals and the precision is not 100 percent accurate.";
printf "\n In any case, we can see that numerically things also work as the symbolic part tells us they should.";


quit;
