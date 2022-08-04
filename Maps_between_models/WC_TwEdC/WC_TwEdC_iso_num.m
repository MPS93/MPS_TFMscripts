//*** Checking that addition is preserved through the isogeny and its doubling property ***//

load "../functions.m";
//load ProjCoord, AffCoord, Ead_AddL_eval, Ew2Ead_iso_eval, Ead2Ew_iso_eval, Ew_AddL_eval

printf "\n\n === ADDITION PRESERVATION AND POINT DOUBLING THROUGH ISOGENY === \n\n";

//--> Data to begin the Ew->Ead example <--//
BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);

lst_w :=  [ 0, 137/8, 0, 11025/256, 0];
a1 := lst_w[1]; a2 := lst_w[2]; a3 := lst_w[3]; a4 := lst_w[4]; a6 := lst_w[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
p_w := [-15/4, 165/32]; P4 := [-105/16, 105/8]; P_w := [p_w[1],p_w[2],1];
q_w := [49/9, -12985/432]; Q_w := [q_w[1],q_w[2],1];

PQ_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w);
pq_w := [PQ_w[1]/PQ_w[3], PQ_w[2]/PQ_w[3]];

p_ad := Ew2Ead_iso_eval(BF,lst_w,p_w);
q_ad := Ew2Ead_iso_eval(BF,lst_w,q_w);
pq_ad := Ew2Ead_iso_eval(BF,lst_w,pq_w);
lst_ad := p_ad[2]; a:= lst_ad[1]; d := lst_ad[2];
p_ad := [p_ad[3][1],p_ad[3][2]];
q_ad := [q_ad[3][1],q_ad[3][2]];
pq_ad := [pq_ad[3][1],pq_ad[3][2]];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(Rw,p_ad), ProjCoord(Rw,q_ad)); r_ad := AffCoord(BF,r_ad);

out_w := Ead2Ew_iso_eval(BF,lst_ad,p_ad); pp_w := out_w[3];
out_w := Ead2Ew_iso_eval(BF,lst_ad,q_ad); qq_w := out_w[3];

PP_w := Ew_AddL_eval(BF,lst_w,P_w,P_w);
QQ_w := Ew_AddL_eval(BF,lst_w,Q_w,Q_w);


printf "\n Here we have an example over the rationals. Let's start by considering this Ew: \n  %o \n",Cw;
printf "\n Two points on this curve are: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n", p_w[1], p_w[2],q_w[1],q_w[2]; 
printf "\n The addition of these two points is: \n\t p_w + q_w = pq_w = (%o, %o) \n", pq_w[1], pq_w[2];
printf "\n Using the W2Ead isogeny, we have that these points go to: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n\t pq_ad = (%o, %o)\n", p_ad[1], p_ad[2], q_ad[1], q_ad[2], pq_ad[1], pq_ad[2];
printf "\n These points are supposed to be points over the Ead with params a = %o, d = %o. Are they really on this curve? \n\t p_ad -> %o, q_ad -> %o, pq_ad -> %o \n", a, d, Rd_var!p_ad in Cad, Rd_var!q_ad in Cad, Rd_var!pq_ad in Cad;
printf "\n Now it's time to add p_ad and q_ad with the addition law for twisted Edwards curves: \n\t p_ad + q_ad = r_ad = (%o, %o) \n", r_ad[1], r_ad[2];
printf "\n Is r_ad equal to pq_ad? That is, does the isogeny preserve addition in this example? \n\t %o \n", r_ad eq pq_ad;
printf "\n Applying the other isogeny now, so the TWE2W, to p_ad and q_ad we get the following points: \n\t pp_w = (%o,%o) \n\t qq_w = (%o,%o) \n", pp_w[1], pp_w[2],qq_w[1], qq_w[2]; 
printf "\n Are the output parameters the same of the original Ew? %o \n", out_w[1] eq lst_w; 
printf "\n Are pp_w, qq_w really on our Ew? \n\t pp_w -> %o, qq_w -> %o \n",Rw_var!pp_w in Cw,Rw_var!qq_w in Cw;
printf "\n Does it hold that pp_w = 2p_w? %o ", pp_w eq [PP_w[1]/PP_w[3],PP_w[2]/PP_w[3]];
printf "\n Does it hold that qq_w = 2q_w? %o \n", qq_w eq [QQ_w[1]/QQ_w[3],QQ_w[2]/QQ_w[3]];


// Data to begin the FF Ew->Ead example
BF := FiniteField(23);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);

a1 := 0; a2 := 10; a3 := 0; a4 := 19; a6 := 0;
lst_w := [BF|a1,a2,a3,a4,a6];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
p_w := [8 , 4]; P_w := [p_w[1],p_w[2],1];
q_w := [9, 13]; Q_w := [q_w[1],q_w[2],1];
PQ_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w);
pq_w := [PQ_w[1]/PQ_w[3], PQ_w[2]/PQ_w[3]];
p_u1 := [7,0]; p_u2 := [6,0];  u_lst := [0,p_u1[1], p_u2[1]];

p_ad := Ew2Ead_iso_eval(BF,lst_w,p_w : u_lst := u_lst);
q_ad := Ew2Ead_iso_eval(BF,lst_w,q_w : u_lst := u_lst);
pq_ad := Ew2Ead_iso_eval(BF,lst_w,pq_w : u_lst := u_lst);
lst_ad := p_ad[2]; a:= lst_ad[1]; d := lst_ad[2];
p_ad := [p_ad[3][1],p_ad[3][2]];
q_ad := [q_ad[3][1],q_ad[3][2]];
pq_ad := [pq_ad[3][1],pq_ad[3][2]];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 

r_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(Rw,p_ad), ProjCoord(Rw,q_ad)); r_ad := AffCoord(BF,r_ad);

out_w := Ead2Ew_iso_eval(BF,lst_ad,p_ad); pp_w := out_w[3];
out_w := Ead2Ew_iso_eval(BF,lst_ad,q_ad); qq_w := out_w[3];

PP_w := Ew_AddL_eval(BF,lst_w,P_w,P_w);
QQ_w := Ew_AddL_eval(BF,lst_w,Q_w,Q_w);

printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over GF(23). Let's start by considering this Ew: \n  %o \n",Cw;
printf "\n Two points on this curve are: \n\t p_w = (%o,%o) \n\t q_w = (%o,%o) \n", p_w[1], p_w[2],q_w[1],q_w[2]; 
printf "\n The addition of these two points is: \n\t p_w + q_w = pq_w = (%o,%o) \n", pq_w[1], pq_w[2];
printf "\n Using the W2Ead isogeny, we have that these points go to: \n\t p_ad = (%o,%o) \n\t q_ad = (%o,%o) \n\t pq_ad = (%o,%o)\n", p_ad[1], p_ad[2], q_ad[1], q_ad[2], pq_ad[1], pq_ad[2];
printf "\n These points are supposed to be points over the Ead with params a = %o, d = %o. Are they really on this curve? \n\t p_ad -> %o, q_ad -> %o, pq_ad -> %o \n", a, d, Rd_var!p_ad in Cad, Rd_var!q_ad in Cad, Rd_var!pq_ad in Cad;
printf "\n Now it's time to add p_ad and q_ad with the addition law for twisted Edwards curves: \n\t p_ad + q_ad = r_ad = (%o, %o) \n", r_ad[1], r_ad[2];
printf "\n Is r_ad equal to pq_ad? That is, does the isogeny preserve addition in this example? \n\t %o \n", r_ad eq pq_ad;
printf "\n Applying the other isogeny now, so the TWE2W, to p_ad and q_ad we get the following points: \n\t pp_w = (%o,%o) \n\t qq_w = (%o,%o) \n", pp_w[1], pp_w[2],qq_w[1], qq_w[2]; 
printf "\n Are the output parameters the same of the original Ew? %o \n", out_w[1] eq lst_w; 
printf "\n Are pp_w, qq_w really on our Ew? \n\t pp_w -> %o, qq_w -> %o \n",Rw_var!pp_w in Cw,Rw_var!qq_w in Cw;
printf "\n Does it hold that pp_w = 2p_w? %o ", pp_w eq [PP_w[1]/PP_w[3],PP_w[2]/PP_w[3]];
printf "\n Does it hold that qq_w = 2q_w? %o \n", qq_w eq [QQ_w[1]/QQ_w[3],QQ_w[2]/QQ_w[3]];


//--> Data to begin the Ead->Ew example <--//
BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
Rw<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(BF,8);


a := 82; d := 17; lst_ad := [a,d];
Cad := Curve(Rd_var, a*x^2 +  y^2 - 1 - d*(x^2)*(y^2)); 
p_ad := [1,9/4];
q_ad := [-72/1393, -1231/1361 ];
PQ_ad := Ead_AddL_eval(BF,lst_ad, ProjCoord(BF,p_ad), ProjCoord(BF,q_ad)); pq_ad := AffCoord(BF,PQ_ad);

p_ad_db := AffCoord(BF,Ead_AddL_eval(BF,[a,d],[[p_ad[1],1],[p_ad[2],1]],[[p_ad[1],1],[p_ad[2],1]]));
q_ad_db := AffCoord(BF,Ead_AddL_eval(BF,[a,d],[[q_ad[1],1],[q_ad[2],1]],[[q_ad[1],1],[q_ad[2],1]]));

p_w := Ead2Ew_iso_eval(BF,lst_ad,p_ad);
q_w := Ead2Ew_iso_eval(BF,lst_ad,q_ad)[3]; q_w := [BF|q_w[1],q_w[2]]; Q_w := [q_w[1],q_w[2],1];
pq_w := Ead2Ew_iso_eval(BF,lst_ad,pq_ad)[3]; pq_w := [BF|pq_w[1],pq_w[2]];
lst_w := p_w[1];u_lst := p_w[2]; p_w := [BF|p_w[3][1],p_w[3][2]]; P_w := [p_w[1],p_w[2],1];
a1 := BF!lst_w[1]; a2 := BF!lst_w[2]; a3 := BF!lst_w[3]; a4 := BF!lst_w[4]; a6 := BF!lst_w[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
R_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w);
r_w := [R_w[1]/R_w[3], R_w[2]/R_w[3]];

out_ad := Ew2Ead_iso_eval(BF,lst_w,p_w: u_lst := u_lst); pp_ad := out_ad[3];
out_ad := Ew2Ead_iso_eval(BF,lst_w,q_w: u_lst := u_lst); qq_ad := out_ad[3];

printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over the rationals. Let's start by considering this Ead: \n  %o \n",Cad;
printf "\n Two points on this curve are: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n", p_ad[1], p_ad[2],q_ad[1],q_ad[2]; 
printf "\n The addition of these two points is: \n\t p_ad + q_ad = pq_ad = (%o, %o) \n", pq_ad[1], pq_ad[2];
printf "\n Using the E2W isogeny, we have that these points go to: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n\t pq_w = (%o, %o)\n", p_w[1], p_w[2], q_w[1], q_w[2], pq_w[1], pq_w[2];
printf "\n These points are supposed to be points over the Ew with params [%o,%o,%o,%o,%o]. Are they really on this curve? \n\t p_w -> %o, q_w -> %o, pq_w -> %o \n", a1,a2,a3,a4,a6, Rw_var!p_w in Cw, Rw_var!q_w in Cw, Rw_var!pq_w in Cw;
printf "\n Now it's time to add p_w and q_w with the addition law for Weierstrass curves: \n\t p_w + q_w = r_w = (%o, %o) \n", r_w[1], r_w[2];
printf "\n Is r_d equal to pq_d? That is, does the isogeny preserve addition in this example? \n\t %o \n", [r_w[1],r_w[2]] eq pq_w;
printf "\n Applying the other isogeny now, so the W2TWE, to p_w and q_w we get the following points: \n\t pp_ad = (%o,%o) \n\t qq_ad = (%o,%o) \n", pp_ad[1], pp_ad[2],qq_ad[1], qq_ad[2]; 
printf "\n Are the output parameters the same of the original Ead? %o \n", out_ad[2] eq lst_ad; 
printf "\n Are pp_ad, qq_ad really on our Ead? \n\t pp_ad -> %o, qq_ad -> %o \n",pp_ad in Cad,qq_ad in Cad;
printf "\n Does it hold that pp_ad = 2p_ad? %o ", pp_ad eq p_ad_db;
printf "\n Does it hold that qq_ad = 2q_ad? %o \n", qq_ad eq q_ad_db;

//--> Data to begin the FF Ead->Ew example <--//
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

p_ad_db := AffCoord(BF,Ead_AddL_eval(BF,[a,d],[[p_ad[1],1],[p_ad[2],1]],[[p_ad[1],1],[p_ad[2],1]]));
q_ad_db := AffCoord(BF,Ead_AddL_eval(BF,[a,d],[[q_ad[1],1],[q_ad[2],1]],[[q_ad[1],1],[q_ad[2],1]]));

p_w := Ead2Ew_iso_eval(BF,lst_ad,p_ad);
q_w := Ead2Ew_iso_eval(BF,lst_ad,q_ad)[3]; q_w := [BF|q_w[1],q_w[2]]; Q_w := [q_w[1],q_w[2],1];
pq_w := Ead2Ew_iso_eval(BF,lst_ad,pq_ad)[3]; pq_w := [BF|pq_w[1],pq_w[2]];
lst_w := p_w[1];u_lst := p_w[2]; p_w := [BF|p_w[3][1],p_w[3][2]]; P_w := [p_w[1],p_w[2],1];
a1 := BF!lst_w[1]; a2 := BF!lst_w[2]; a3 := BF!lst_w[3]; a4 := BF!lst_w[4]; a6 := BF!lst_w[5];
Cw := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); //Cw;
R_w := Ew_AddL_eval(BF,lst_w,P_w,Q_w);
r_w := [R_w[1]/R_w[3], R_w[2]/R_w[3]];

out_ad := Ew2Ead_iso_eval(BF,lst_w,p_w: u_lst := u_lst); pp_ad := out_ad[3];
out_ad := Ew2Ead_iso_eval(BF,lst_w,q_w: u_lst := u_lst); qq_ad := out_ad[3];

printf "\n\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
printf "\n Here we have an example over the finite field GF(2^382 - 105). Let's start by considering this Ead: \n  %o \n",Cad;
printf "\n Two points on this curve are: \n\t p_ad = (%o,%o) \n\t q_ad = (%o, %o) \n", p_ad[1], p_ad[2],q_ad[1],q_ad[2]; 
printf "\n The addition of these two points is: \n\t p_ad + q_ad = pq_ad = (%o, %o) \n", pq_ad[1], pq_ad[2];
printf "\n Using the E2W isogeny, we have that these points go to: \n\t p_w = (%o,%o) \n\t q_w = (%o, %o) \n\t pq_w = (%o, %o)\n", p_w[1], p_w[2], q_w[1], q_w[2], pq_w[1], pq_w[2];
printf "\n These points are supposed to be points over the Ew with params [%o,%o,%o,%o,%o]. Are they really on this curve? \n\t p_w -> %o, q_w -> %o, pq_w -> %o \n", a1,a2,a3,a4,a6, Rw_var!p_w in Cw, Rw_var!q_w in Cw, Rw_var!pq_w in Cw;
printf "\n Now it's time to add p_w and q_w with the addition law for Weierstrass curves: \n\t p_w + q_w = r_w = (%o, %o) \n", r_w[1], r_w[2];
printf "\n Is r_d equal to pq_d? That is, does the isogeny preserve addition in this example? \n\t %o \n", [r_w[1],r_w[2]] eq pq_w;
printf "\n Applying the other isogeny now, so the W2TWE, to p_w and q_w we get the following points: \n\t pp_ad = (%o,%o) \n\t qq_ad = (%o,%o) \n", pp_ad[1], pp_ad[2],qq_ad[1], qq_ad[2]; 
printf "\n Are the output parameters the same of the original Ead? %o \n", out_ad[2] eq lst_ad; 
printf "\n Are pp_ad, qq_ad really on our Ead? \n\t pp_ad -> %o, qq_ad -> %o \n",pp_ad in Cad,qq_ad in Cad;
printf "\n Does it hold that pp_ad = 2p_ad? %o ", pp_ad eq p_ad_db;
printf "\n Does it hold that qq_ad = 2q_ad? %o \n", qq_ad eq q_ad_db;

quit;
