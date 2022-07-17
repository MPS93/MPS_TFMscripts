//*** Examples and numerical checks on the Ew-EdC  birational equivalence without P4 ***//

load "../functions.m";
//load ProjCoord, AffCoord, Ew2Ed_eval, Ed2Ew_eval, Ead_AddL_eval, Ew_AddL_eval


printf "\n\n === NUM. CHECKS ABOUT THE PRESERVATION OF ADDITION BY THE Ew_EdC BIRAT. EQ. WO/ P4 === \n";

//--> Data to begin the Ew->EdC example <--//
p := 2^448 - 2^224 - 1;
BF := FiniteField(p);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);

lst := [BF| 0, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018287279, 0, 1527402724, 0];
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6 := lst[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
lst := [a1,a2,a3,a4,a6];
p_w := [161519716510134864566516401752889896523031413486070680062553377595691628481495727263643640821928121317797529880409697191345221781815562,717225809973557004946798282930174671419297352528066850822411768352295942696132428073016719856723894754274572935636178936371882792255446];
q_w := [39082, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018287275];
P_w := [p_w[1],p_w[2],1]; Q_w := [q_w[1],q_w[2],1];
PQ_w := Ew_AddL_eval(BF,lst,P_w,Q_w);
pq_w := [PQ_w[1]/PQ_w[3], PQ_w[2]/PQ_w[3]];

p_d := Ew2Ed_eval(BF,lst,p_w);
q_d := Ew2Ed_eval(BF,lst,q_w);
pq_d := Ew2Ed_eval(BF,lst,pq_w);
d := p_d[1][1]; lst_d := [BF!1,d];
p_d := [p_d[2][1],p_d[2][2]];
q_d := [q_d[2][1],q_d[2][2]];
pq_d := [pq_d[2][1],pq_d[2][2]];
Cd := Curve(Rd_var, x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_d := Ead_AddL_eval(BF,lst_d, ProjCoord(Rw,p_d), ProjCoord(Rw,q_d)); r_d :=  AffCoord(BF,r_d);
PP4 := [39082, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018287275];
// Ed2Ew_eval(BF,[d],PP4,r_d);

printf "\n Here we have an example over the finite field GF(2^448 - 2^224 - 1). Let's start by considering this Ew: \n  %o \n",Cw;
printf "\n Two points on this curve are: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n", p_w[1], p_w[2],q_w[1],q_w[2]; 
printf "\n The addition of these two points is: \n\t p_w + q_w = pq_w = (%o, %o) \n", pq_w[1], pq_w[2];
printf "\n Using the Ew2Ed birational equivalence, we have that these points go to: \n\t p_d = (%o,%o) \n\t q_d = (%o, %o) \n\t pq_d = (%o, %o)\n", p_d[1], p_d[2], q_d[1], q_d[2], pq_d[1], pq_d[2];
printf "\n These points are supposed to be points over the EdC with param d = %o. Are they really on this curve? \n\t p_d -> %o, q_d -> %o, pq_d -> %o \n", d, Rd_var!p_d in Cd, Rd_var!q_d in Cd, Rd_var!pq_d in Cd;
printf "\n Now it's time to add p_d and q_d with the addition law for Edwards curves: \n\t p_d + q_d = r_d = (%o, %o) \n", r_d[1], r_d[2];
printf "\n Is r_d equal to pq_d? That is, does the birational equivalence preserve addition in this example? \n\t %o", r_d eq pq_d;

// Ew -> EdC example over the rationals
BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);

lst := [BF| 0, 48/5, 0, 196/25, 0 ];
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6 := lst[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
lst := [a1,a2,a3,a4,a6];
p_w := [1/5, -7/5 ]; q_w := [-14/9, 364/135];
P_w := [p_w[1],p_w[2],1]; Q_w := [q_w[1],q_w[2],1];
PQ_w := Ew_AddL_eval(BF,lst,P_w,Q_w);
pq_w := [PQ_w[1]/PQ_w[3], PQ_w[2]/PQ_w[3]];

p_d := Ew2Ed_eval(BF,lst,p_w);
q_d := Ew2Ed_eval(BF,lst,q_w);
pq_d := Ew2Ed_eval(BF,lst,pq_w);
d := p_d[1][1]; lst_d := [BF!1,d]; 
p_d := [p_d[2][1],p_d[2][2]];
q_d := [q_d[2][1],q_d[2][2]];
pq_d := [pq_d[2][1],pq_d[2][2]];
Cd := Curve(Rd_var, x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_d := Ead_AddL_eval(BF,lst_d, ProjCoord(Rw,p_d), ProjCoord(Rw,q_d)); r_d := AffCoord(BF,r_d);
// Ed2Ew_eval(BF,[d],PP4,r_d);

printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over the rationals. Let's start by considering this Ew: \n  %o \n",Cw;
printf "\n Two points on this curve are: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n", p_w[1], p_w[2],q_w[1],q_w[2]; 
printf "\n The addition of these two points is: \n\t p_w + q_w = pq_w = (%o, %o) \n", pq_w[1], pq_w[2];
printf "\n Using the Ew2Ed birational equivalence, we have that these points go to: \n\t p_d = (%o,%o) \n\t q_d = (%o, %o) \n\t pq_d = (%o, %o)\n", p_d[1], p_d[2], q_d[1], q_d[2], pq_d[1], pq_d[2];
printf "\n These points are supposed to be points over the EdC with param d = %o. Are they really on this curve? \n\t p_d -> %o, q_d -> %o, pq_d -> %o \n", d, Rd_var!p_d in Cd, Rd_var!q_d in Cd, Rd_var!pq_d in Cd;
printf "\n Now it's time to add p_d and q_d with the addition law for Edwards curves: \n\t p_d + q_d = r_d = (%o, %o) \n", r_d[1], r_d[2];
printf "\n Is r_d equal to pq_d? That is, does the birational equivalence preserve addition in this example? \n\t %o", r_d eq pq_d;


//--> Data to begin the Ew->EdC example <--//
BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);

lst := [BF| 27, -23321/132, 98, -1440071/1089, -2401 ]; u2p := 0;
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6 := lst[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
lst := [a1,a2,a3,a4,a6];
p_w := [ 144/25,-201087/1375 ];
q_w := [ 30688841142001/43970426240400, -16489454909765872019419/291568775808617208000];
P_w := [p_w[1],p_w[2],1]; Q_w := [q_w[1],q_w[2],1];
PQ_w := Ew_AddL_eval(BF,lst,P_w,Q_w);
pq_w := [PQ_w[1]/PQ_w[3], PQ_w[2]/PQ_w[3]];

p_d := Ew2Ed_eval(BF,lst,p_w);
q_d := Ew2Ed_eval(BF,lst,q_w);
pq_d := Ew2Ed_eval(BF,lst,pq_w);
d := p_d[1][1]; lst_d := [BF!1,d];
p_d := [p_d[2][1],p_d[2][2]];
q_d := [q_d[2][1],q_d[2][2]];
pq_d := [pq_d[2][1],pq_d[2][2]];
Cd := Curve(Rd_var, x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_d := Ead_AddL_eval(BF,lst_d, ProjCoord(Rw,p_d), ProjCoord(Rw,q_d)); r_d :=  AffCoord(BF,r_d);

printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over the rationals. Let's start by considering this complete Ew: \n  %o \n",Cw;
printf "\n Two points on this curve are: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n", p_w[1], p_w[2],q_w[1],q_w[2]; 
printf "\n The addition of these two points is: \n\t p_w + q_w = pq_w = (%o, %o) \n", pq_w[1], pq_w[2];
printf "\n Using the Ew2Ed birational equivalence, we have that these points go to: \n\t p_d = (%o,%o) \n\t q_d = (%o, %o) \n\t pq_d = (%o, %o)\n", p_d[1], p_d[2], q_d[1], q_d[2], pq_d[1], pq_d[2];
printf "\n These points are supposed to be points over the EdC with param d = %o. Are they really on this curve? \n\t p_d -> %o, q_d -> %o, pq_d -> %o \n", d, Rd_var!p_d in Cd, Rd_var!q_d in Cd, Rd_var!pq_d in Cd;
printf "\n Now it's time to add p_d and q_d with the addition law for Edwards curves: \n\t p_d + q_d = r_d = (%o, %o) \n", r_d[1], r_d[2];
printf "\n Is r_d equal to pq_d? That is, does the birational equivalence preserve addition in this example? \n\t %o", r_d eq pq_d;

//--> Data to begin the EdC->Ew example <--//
p := 2^251 - 9;
BF := FiniteField(p);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);


d := -1174; lst := [d];
Cd := Curve(Rd_var, x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 
p_d := [1582619097725911541954547006453739763381091388846394833492296309729998839514,
3037538013604154504764115728651437646519513534305223422754827055689195992590];
P_d := [[p_d[1],1],[p_d[2],1]];
q_d := [2549750364055938208658959522391484219652589930931770035797162350828388088291,
209746268288649404466725102591950163741896080620509071752107886722442588196];
PQ_d := Ead_AddL_eval(BF,[1,d], ProjCoord(Rw,p_d), ProjCoord(Rw,q_d)); pq_d := AffCoord(BF,PQ_d);

p_w := Ed2Ew_eval(BF,lst,p_d);  
q_w := Ed2Ew_eval(BF,lst,q_d)[2]; q_w := [BF|q_w[1],q_w[2]];
pq_w := Ed2Ew_eval(BF,lst,pq_d)[2];pq_w := [BF|pq_w[1],pq_w[2]];
lst_w := p_w[1]; p_w := [BF|p_w[2][1],p_w[2][2]];
a1 := BF!lst_w[1]; a2 := BF!lst_w[2]; a3 := BF!lst_w[3]; a4 := BF!lst_w[4]; a6 := BF!lst_w[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
P_w := [p_w[1],p_w[2],1]; Q_w := [q_w[1],q_w[2],1];
R_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w);
r_w := [R_w[1]/R_w[3], R_w[2]/R_w[3]];

printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over the finite field GF(2^251 - 9). Let's start by considering this EdC: \n  %o \n",Cd;
printf "\n Two points on this curve are: \n\t p_d = (%o,%o) \n\t q_d = (%o, %o) \n", p_d[1], p_d[2],q_d[1],q_d[2]; 
printf "\n The addition of these two points is: \n\t p_d + q_d = pq_d = (%o, %o) \n", pq_d[1], pq_d[2];
printf "\n Using the Ed2Ew birational equivalence, we have that these points go to: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n\t pq_w = (%o, %o)\n", p_w[1], p_w[2], q_w[1], q_w[2], pq_w[1], pq_w[2];
printf "\n These points are supposed to be points over the Ew with params [%o,%o,%o,%o,%o]. Are they really on this curve? \n\t p_w -> %o, q_w -> %o, pq_w -> %o \n", a1,a2,a3,a4,a6, Rw_var!p_w in Cw, Rw_var!q_w in Cw, Rw_var!pq_w in Cw;
printf "\n Now it's time to add p_w and q_w with the addition law for Weierstrass curves: \n\t p_w + q_w = r_w = (%o, %o) \n", r_w[1], r_w[2];
printf "\n Is r_d equal to pq_d? That is, does the birational equivalence preserve addition in this example? \n\t %o", r_w eq pq_w;


//--> Data to begin the EdC->Ew example over the rationals <--//

BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);


d := 36; lst := [d];
Cd := Curve(Rd_var, x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 
p_d := [-1/8, -3/2];
P_d := [[p_d[1],1],[p_d[2],1]];
q_d := [-16919760/430106689,-429773761/417954239];
PQ_d := Ead_AddL_eval(BF,[1,d], ProjCoord(Rw,p_d), ProjCoord(Rw,q_d)); pq_d := AffCoord(BF,PQ_d);

p_w := Ed2Ew_eval(BF,lst,p_d);  
q_w := Ed2Ew_eval(BF,lst,q_d)[2]; q_w := [BF|q_w[1],q_w[2]];
pq_w := Ed2Ew_eval(BF,lst,pq_d)[2];pq_w := [BF|pq_w[1],pq_w[2]];
lst_w := p_w[1]; p_w := [BF|p_w[2][1],p_w[2][2]];
a1 := BF!lst_w[1]; a2 := BF!lst_w[2]; a3 := BF!lst_w[3]; a4 := BF!lst_w[4]; a6 := BF!lst_w[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
P_w := [p_w[1],p_w[2],1]; Q_w := [q_w[1],q_w[2],1];
R_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w);
r_w := [R_w[1]/R_w[3], R_w[2]/R_w[3]];


printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over the rationals. Let's start by considering this EdC: \n  %o \n",Cd;
printf "\n Two points on this curve are: \n\t p_d = (%o,%o) \n\t q_d = (%o, %o) \n", p_d[1], p_d[2],q_d[1],q_d[2]; 
printf "\n The addition of these two points is: \n\t p_d + q_d = pq_d = (%o, %o) \n", pq_d[1], pq_d[2];
printf "\n Using the Ed2Ew birational equivalence, we have that these points go to: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n\t pq_w = (%o, %o)\n", p_w[1], p_w[2], q_w[1], q_w[2], pq_w[1], pq_w[2];
printf "\n These points are supposed to be points over the Ew with params [%o,%o,%o,%o,%o]. Are they really on this curve? \n\t p_w -> %o, q_w -> %o, pq_w -> %o \n", a1,a2,a3,a4,a6, Rw_var!p_w in Cw, Rw_var!q_w in Cw, Rw_var!pq_w in Cw;
printf "\n Now it's time to add p_w and q_w with the addition law for Weierstrass curves: \n\t p_w + q_w = r_w = (%o, %o) \n", r_w[1], r_w[2];
printf "\n Is r_d equal to pq_d? That is, does the birational equivalence preserve addition in this example? \n\t %o", r_w eq pq_w;

BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);

d := 59/33; lst := [d];
Cd := Curve(Rd_var, x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 
p_d := [-1980/3349, 2701/2051];
q_d := [-73468399352040/98315324391601, -65332207270801/3954524986799];
PQ_d := Ead_AddL_eval(BF,[1,d], ProjCoord(Rw,p_d), ProjCoord(Rw,q_d)); pq_d := AffCoord(BF,PQ_d);
a1 := 27; a3 := 98; u2p := 0;


p_w := Ed2Ew_eval(BF,lst,p_d : PP4 := [0,0], a1 := a1, a3 := a3);  
q_w := Ed2Ew_eval(BF,lst,q_d : PP4 := [0,0], a1 := a1, a3 := a3)[2]; q_w := [BF|q_w[1],q_w[2]];
pq_w := Ed2Ew_eval(BF,lst,pq_d : PP4 := [0,0], a1 := a1, a3 := a3)[2];pq_w := [BF|pq_w[1],pq_w[2]];
lst_w := p_w[1]; p_w := [BF|p_w[2][1],p_w[2][2]];
a1 := BF!lst_w[1]; a2 := BF!lst_w[2]; a3 := BF!lst_w[3]; a4 := BF!lst_w[4]; a6 := BF!lst_w[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
P_w := [p_w[1],p_w[2],1]; Q_w := [q_w[1],q_w[2],1];
R_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w);
r_w := [R_w[1]/R_w[3], R_w[2]/R_w[3]];

printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over the rationals. Let's start by considering this EdC: \n  %o \n",Cd;
printf "\n Two points on this curve are: \n\t p_d = (%o,%o) \n\t q_d = (%o, %o) \n", p_d[1], p_d[2],q_d[1],q_d[2]; 
printf "\n The addition of these two points is: \n\t p_d + q_d = pq_d = (%o, %o) \n", pq_d[1], pq_d[2];
printf "\n Using the Ed2Ew birational equivalence, we have that these points go to: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n\t pq_w = (%o, %o)\n", p_w[1], p_w[2], q_w[1], q_w[2], pq_w[1], pq_w[2];
printf "\n These points are supposed to be points over the Ew with params [%o,%o,%o,%o,%o]. Are they really on this curve? \n\t p_w -> %o, q_w -> %o, pq_w -> %o \n", a1,a2,a3,a4,a6, Rw_var!p_w in Cw, Rw_var!q_w in Cw, Rw_var!pq_w in Cw;
printf "\n Now it's time to add p_w and q_w with the addition law for Weierstrass curves: \n\t p_w + q_w = r_w = (%o, %o) \n", r_w[1], r_w[2];
printf "\n Is r_d equal to pq_d? That is, does the birational equivalence preserve addition in this example? \n\t %o", r_w eq pq_w;

quit;
