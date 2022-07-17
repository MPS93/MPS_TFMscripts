//*** Checking that addition is preserved through the MD birat. eq. ***//

load "../functions.m";
//load ProjCoord, AffCoord, Ew2Ead_eval, Ead2Ew_eval, Ead_AddL_eval, Ew_AddL_eval

printf "\n\n === ADDITION PRESERVATION THROUGH MD === \n\n";

//--> Data to begin the Ew->Ead example <--//
BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);

lst_w :=  [ 0, 137/8, 0, 11025/256, 0];
a1 := lst_w[1]; a2 := lst_w[2]; a3 := lst_w[3]; a4 := lst_w[4]; a6 := lst_w[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
p_w := [-15/4, 165/32]; P4 := p_w; P_w := [p_w[1],p_w[2],1];
q_w := [49/9, -12985/432]; Q_w := [q_w[1],q_w[2],1];
//[-105/16, 105/8], [-105/16, -105/8], [105/16, 1155/32], [105/16, -1155/32] are of order 4
PQ_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w);
pq_w := [PQ_w[1]/PQ_w[3], PQ_w[2]/PQ_w[3]];

p_ad := Ew2Ead_eval(BF,lst_w,p_w);
q_ad := Ew2Ead_eval(BF,lst_w,q_w);
pq_ad := Ew2Ead_eval(BF,lst_w,pq_w);
lst_ad := p_ad[2]; a:= lst_ad[1]; d := lst_ad[2];
p_ad := [p_ad[3][1],p_ad[3][2]];
q_ad := [q_ad[3][1],q_ad[3][2]];
pq_ad := [pq_ad[3][1],pq_ad[3][2]];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(Rw,p_ad), ProjCoord(Rw,q_ad)); r_ad := AffCoord(BF,r_ad);

printf "\n Here we have an example over the rationals. Let's start by considering this Ew: \n  %o \n",Cw;
printf "\n Two points on this curve are: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n", p_w[1], p_w[2],q_w[1],q_w[2]; 
printf "\n The addition of these two points is: \n\t p_w + q_w = pq_w = (%o, %o) \n", pq_w[1], pq_w[2];
printf "\n Using the W2Ead birational equivalence, we have that these points go to: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n\t pq_ad = (%o, %o)\n", p_ad[1], p_ad[2], q_ad[1], q_ad[2], pq_ad[1], pq_ad[2];
printf "\n These points are supposed to be points over the Ead with params a= %o, d = %o. Are they really on this curve? \n\t p_ad -> %o, q_ad -> %o, pq_ad -> %o \n", a, d, Rd_var!p_ad in Cad, Rd_var!q_ad in Cad, Rd_var!pq_ad in Cad;
printf "\n Now it's time to add p_ad and q_ad with the addition law for twisted Edwards curves: \n\t p_ad + q_ad = r_ad = (%o, %o) \n", r_ad[1], r_ad[2];
printf "\n Is r_ad equal to pq_ad? That is, does the birational equivalence preserve addition in this example? \n\t %o", r_ad eq pq_ad;


p_ad := Ew2Ead_eval(BF,lst_w,p_w: a := 25);
q_ad := Ew2Ead_eval(BF,lst_w,q_w: a := 25);
pq_ad := Ew2Ead_eval(BF,lst_w,pq_w: a := 25);
lst_ad := p_ad[2]; a:= lst_ad[1]; d := lst_ad[2];
p_ad := [p_ad[3][1],p_ad[3][2]];
q_ad := [q_ad[3][1],q_ad[3][2]];
pq_ad := [pq_ad[3][1],pq_ad[3][2]];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(Rw,p_ad), ProjCoord(Rw,q_ad)); r_ad := AffCoord(BF,r_ad);

printf "\n\n\n The same Ew can also be birat. equivalent to a Ead with a != 1. So using the same Ew and the same points on it, let's see how the W2Ead birational equivalence works for a = 25.\n";
printf "\n We have that p_w, q_w and pq_w go to: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n\t pq_ad = (%o, %o)\n", p_ad[1], p_ad[2], q_ad[1], q_ad[2], pq_ad[1], pq_ad[2];
printf "\n These points are supposed to be points over the Ead with params a= %o, d = %o. Are they really on this curve? \n\t p_ad -> %o, q_ad -> %o, pq_ad -> %o \n", a, d, Rd_var!p_ad in Cad, Rd_var!q_ad in Cad, Rd_var!pq_ad in Cad;
printf "\n Now it's time to add p_ad and q_ad with the addition law for twisted Edwards curves: \n\t p_ad + q_ad = r_ad = (%o, %o) \n", r_ad[1], r_ad[2];
printf "\n Is r_ad equal to pq_ad? That is, does the birat. eq. preserve addition in this example? \n\t %o", r_ad eq pq_ad;


p_ad := Ew2Ead_eval(BF,lst_w,p_w: a := 5);
BF := p_ad[1]; 
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);
q_ad := Ew2Ead_eval(BF,lst_w,q_w: a := 5);
pq_ad := Ew2Ead_eval(BF,lst_w,pq_w: a := 5);
lst_ad := p_ad[2]; a:= lst_ad[1]; d := lst_ad[2];
p_ad := [BF|p_ad[3][1],p_ad[3][2]];
q_ad := [BF|q_ad[3][1],q_ad[3][2]];
pq_ad := [BF|pq_ad[3][1],pq_ad[3][2]];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(BF,p_ad), ProjCoord(BF,q_ad)); r_ad := AffCoord(BF,r_ad);

printf "\n\n\n What happens if 'a' is not a square on BF? Well, then we have to work on the algebraic closure of BF. So using the same Ew and the same points on it, let's see how the W2Ead birational equivalence works for a = 5.\n";
printf "\n We have that p_w, q_w and pq_w go to: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n\t pq_ad = (%o, %o)\n", p_ad[1], p_ad[2], q_ad[1], q_ad[2], pq_ad[1], pq_ad[2];
printf "\n Here 'r1' represents the square root of 'a'. And the BF has been changed to: \n\t %o \n", BF;
printf "\n These points are supposed to be points over the Ead with params a= %o, d = %o. Are they really on this curve? \n\t p_ad -> %o, q_ad -> %o, pq_ad -> %o \n", a, d, Rd_var!p_ad in Cad, Rd_var!q_ad in Cad, Rd_var!pq_ad in Cad;
printf "\n Now it's time to add p_ad and q_ad with the addition law for twisted Edwards curves: \n\t p_ad + q_ad = r_ad = (%o, %o) \n", r_ad[1], r_ad[2];
printf "\n Is r_ad equal to pq_ad? That is, does the birat. eq. preserve addition in this example? \n\t %o", r_ad eq pq_ad;


//--> Data to begin the FF Ew->Ead example <--//
BF := FiniteField(23);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);

a1 := 0; a2 := 12; a3 := 0; a4 := 16; a6 := 0;
lst_w := [BF|a1,a2,a3,a4,a6];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
p_w := [1 , 12]; P_w := [p_w[1],p_w[2],1];
q_w := [15, 17]; Q_w := [q_w[1],q_w[2],1];
// [19,8], [19,15] are of order 4
PQ_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w);
pq_w := [PQ_w[1]/PQ_w[3], PQ_w[2]/PQ_w[3]];

p_ad := Ew2Ead_eval(BF,lst_w,p_w);
q_ad := Ew2Ead_eval(BF,lst_w,q_w);
pq_ad := Ew2Ead_eval(BF,lst_w,pq_w);
lst_ad := p_ad[2]; a:= lst_ad[1]; d := lst_ad[2];
p_ad := [p_ad[3][1],p_ad[3][2]];
q_ad := [q_ad[3][1],q_ad[3][2]];
pq_ad := [pq_ad[3][1],pq_ad[3][2]];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(Rw,p_ad), ProjCoord(Rw,q_ad)); r_ad := AffCoord(BF,r_ad);

printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over GF(23). Let's start by considering this Ew: \n  %o \n",Cw;
printf "\n Two points on this curve are: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n", p_w[1], p_w[2],q_w[1],q_w[2]; 
printf "\n The addition of these two points is: \n\t p_w + q_w = pq_w = (%o, %o) \n", pq_w[1], pq_w[2];
printf "\n Using the W2Ead birational equivalence, we have that these points go to: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n\t pq_ad = (%o, %o)\n", p_ad[1], p_ad[2], q_ad[1], q_ad[2], pq_ad[1], pq_ad[2];
printf "\n These points are supposed to be points over the Ead with params a= %o, d = %o. Are they really on this curve? \n\t p_ad -> %o, q_ad -> %o, pq_ad -> %o \n", a, d, Rd_var!p_ad in Cad, Rd_var!q_ad in Cad, Rd_var!pq_ad in Cad;
printf "\n Now it's time to add p_ad and q_ad with the addition law for twisted Edwards curves: \n\t p_ad + q_ad = r_ad = (%o, %o) \n", r_ad[1], r_ad[2];
printf "\n Is r_ad equal to pq_ad? That is, does the birational equivalence preserve addition in this example? \n\t %o", r_ad eq pq_ad;


p_ad := Ew2Ead_eval(BF,lst_w,p_w: a := 3);
q_ad := Ew2Ead_eval(BF,lst_w,q_w: a := 3);
pq_ad := Ew2Ead_eval(BF,lst_w,pq_w: a := 3);
lst_ad := p_ad[2]; a:= lst_ad[1]; d := lst_ad[2];
p_ad := [p_ad[3][1],p_ad[3][2]];
q_ad := [q_ad[3][1],q_ad[3][2]];
pq_ad := [pq_ad[3][1],pq_ad[3][2]];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(Rw,p_ad), ProjCoord(Rw,q_ad)); r_ad := AffCoord(BF,r_ad);

printf "\n\n\n The same Ew can also be birat. equivalent to a Ead with a != 1. So using the same Ew and the same points on it, let's see how the W2Ead birational equivalence works for a = 25.\n";
printf "\n We have that p_w, q_w and pq_w go to: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n\t pq_ad = (%o, %o)\n", p_ad[1], p_ad[2], q_ad[1], q_ad[2], pq_ad[1], pq_ad[2];
printf "\n These points are supposed to be points over the Ead with params a= %o, d = %o. Are they really on this curve? \n\t p_ad -> %o, q_ad -> %o, pq_ad -> %o \n", a, d, Rd_var!p_ad in Cad, Rd_var!q_ad in Cad, Rd_var!pq_ad in Cad;
printf "\n Now it's time to add p_ad and q_ad with the addition law for twisted Edwards curves: \n\t p_ad + q_ad = r_ad = (%o, %o) \n", r_ad[1], r_ad[2];
printf "\n Is r_ad equal to pq_ad? That is, does the birat. eq. preserve addition in this example? \n\t %o", r_ad eq pq_ad;


p_ad := Ew2Ead_eval(BF,lst_w,p_w: a := 11);
BF := p_ad[1];
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);
q_ad := Ew2Ead_eval(BF,lst_w,q_w: a := 11);
pq_ad := Ew2Ead_eval(BF,lst_w,pq_w: a := 11);
lst_ad := p_ad[2]; a:= lst_ad[1]; d := lst_ad[2];
p_ad := [BF|p_ad[3][1],p_ad[3][2]];
q_ad := [BF|q_ad[3][1],q_ad[3][2]];
pq_ad := [BF|pq_ad[3][1],pq_ad[3][2]];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(BF,p_ad), ProjCoord(BF,q_ad)); r_ad := AffCoord(BF,r_ad);

printf "\n\n\n What happens if 'a' is not a square on BF? Well, then we have to work on the algebraic closure of BF. So using the same Ew and the same points on it, let's see how the W2Ead birational equivalence works for a = 11.\n";
printf "\n We have that p_w, q_w and pq_w go to: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n\t pq_ad = (%o, %o)\n", p_ad[1], p_ad[2], q_ad[1], q_ad[2], pq_ad[1], pq_ad[2];
printf "\n Here 'r1' represents the square root of 'a'. And the BF has been changed to: \n\t %o \n", BF;
printf "\n These points are supposed to be points over the Ead with params a= %o, d = %o. Are they really on this curve? \n\t p_ad -> %o, q_ad -> %o, pq_ad -> %o \n", a, d, Rd_var!p_ad in Cad, Rd_var!q_ad in Cad, Rd_var!pq_ad in Cad;
printf "\n Now it's time to add p_ad and q_ad with the addition law for twisted Edwards curves: \n\t p_ad + q_ad = r_ad = (%o, %o) \n", r_ad[1], r_ad[2];
printf "\n Is r_ad equal to pq_ad? That is, does the birat. eq. preserve addition in this example? \n\t %o", r_ad eq pq_ad;

//--> Data to begin the Ead->Ew example <--//
p := 2^382 - 105;
BF := FiniteField(p);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);


a := 216547; d := BF!a*(-67254); lst_ad := [a,d];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 
p_ad := [1752167411121869277917751722267788655961702307508145952969064921766738032475012772864787403207074614859447249041651, 17];
P_ad := [[p_ad[1],1],[p_ad[2],1]];
q_ad := [ 2608701730757893826068041123479891646848773283990489975691454083477968234893211605423340254939607215612005254071209, 981565035501027873780508780039633250980186172985292787536388573334951241433846757449675882375747228241007013820068];
PQ_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(BF,p_ad), ProjCoord(BF,q_ad)); pq_ad := AffCoord(BF,PQ_ad);
// p := [BF|67255, 4925250774549309901534880012517951725634967408808180833493536675530715221437151326426783281860614455100818821787429]; is of order 4

p_w := Ead2Ew_eval(BF,lst_ad,p_ad : a1 := 287765, a3 := 97646);  
q_w := Ead2Ew_eval(BF,lst_ad,q_ad : a1 := 287765, a3 := 97646)[3]; q_w := [BF|q_w[1],q_w[2]]; Q_w := [BF|q_w[1],q_w[2],1];
pq_w := Ead2Ew_eval(BF,lst_ad,pq_ad : a1 := 287765, a3 := 97646)[3]; pq_w := [BF|pq_w[1],pq_w[2]];
lst_w := p_w[2]; p_w := [BF|p_w[3][1],p_w[3][2]]; P_w := [BF|p_w[1],p_w[2],1];
a1 := BF!lst_w[1]; a2 := BF!lst_w[2]; a3 := BF!lst_w[3]; a4 := BF!lst_w[4]; a6 := BF!lst_w[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
R_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w); 
r_w := [BF!(R_w[1])/BF!(R_w[3]), BF!(R_w[2])/BF!(R_w[3])];

printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over the finite field GF(2^382 - 105). Let's start by considering this Ead: \n  %o \n",Cad;
printf "\n Two points on this curve are: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n", p_ad[1], p_ad[2],q_ad[1],q_ad[2]; 
printf "\n The addition of these two points is: \n\t p_ad + q_ad = pq_ad = (%o, %o) \n", pq_ad[1], pq_ad[2];
printf "\n Using the E2W birational equivalence, we have that these points go to: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n\t pq_w = (%o, %o)\n", p_w[1], p_w[2], q_w[1], q_w[2], pq_w[1], pq_w[2];
printf "\n These points are supposed to be points over the Ew with params [%o,%o,%o,%o,%o]. Are they really on this curve? \n\t p_w -> %o, q_w -> %o, pq_w -> %o \n", a1,a2,a3,a4,a6, Rw_var!p_w in Cw, Rw_var!q_w in Cw, Rw_var!pq_w in Cw;
printf "\n Now it's time to add p_w and q_w with the addition law for Weierstrass curves: \n\t p_w + q_w = r_w = (%o, %o) \n", r_w[1], r_w[2];
printf "\n Is r_d equal to pq_d? That is, does the birational equivalence preserve addition in this example? \n\t %o", r_w eq pq_w;

quit;
