//*** Implementing the birat. eq. between Edwards and Weierstrass ***//
//-> Here we will use the BL  birat. eq. <-//

load "../functions.m";
//load Ew2Ed_P4_symb, Ed2Ew_P4_symb, Ew2Ed_P4_eval, Ed2Ew_P4_eval


printf "\n\n BL BIRATIONAL EQUIVALENCE BETWEEN Ed AND Ew \n\n";

//--> Symbolic examples <--//

BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := FunctionField(R,2);
Rd_var<x,y> := FunctionField(R,2);

xx := Evaluate(Ew2Ed_P4_symb(BF,[0,a2,0,a4,a6],[up,vp])[2][1],[u,v]);
yy := Evaluate(Ew2Ed_P4_symb(BF,[0,a2,0,a4,a6],[up,vp])[2][2],[u,v]);


print "\n Change of coordinates from Ew to Ed: \n";
" x := ", xx; "\n y := ", yy;

uu := Evaluate(Ed2Ew_P4_symb(BF,[d],[up,vp])[2][1],[x,y]);
vv := Evaluate(Ed2Ew_P4_symb(BF,[d],[up,vp])[2][2],[x,y]); 

print "\n Change of coordinates from Ed to Ew: \n";
" u := ", uu; "\n v := ", vv;

print "\n Are these two maps inverse of each other? \n";
"  ", u eq Evaluate(uu,[xx,yy]) and v eq Evaluate(vv,[xx,yy]) and x eq Evaluate(xx,[uu,vv]) and y eq Evaluate(yy,[uu,vv]);

"\n===============\n";


//--> Numerical examples <--//

BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
a1 := 0; a2 := 10; a3 := 0; a4 := 9; a6 := 0; 
C := Curve(Rw_var, v^2 - u^3 - a2*u^2 - a4*u - a6); C;
lst := [a1,a2,a3,a4,a6];
P4 := [-3, -6]; p:= [-3, 6];

"\n The given Ew above and its point p_w := ", Rw_var!p, " correspond to the following Ed and p_ad \n";

ch := Ew2Ed_P4_eval(BF,lst,P4,p);

d := ch[1][1];
C := Curve(Rd_var, 1 + d*x^2*y^2 - x^2 - y^2); C;
p_ad := Rd_var![ch[2][1],ch[2][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in Ed? \n\t", tf;

"\n===============\n";

BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
d:= 4; lst := [d];
C := Curve(Rd_var, 1 + d*x^2*y^2 - x^2 - y^2); C;	
P4 := [-3,-6];
p := [-1,0];

"\n The given Ed above and its point p_ad := ", Rd_var!p, " correspond to the following Ew and p_w \n";

ch := Ed2Ew_P4_eval(BF,lst,P4,p);
a1 := ch[1][1]; a2 := ch[1][2]; a3 := ch[1][3]; a4 := ch[1][4]; a6 := ch[1][5]; 
C := Curve(Rw_var, v^2 - u^3 - a2*u^2 - a4*u - a6); C;
p_w := Rw_var![ch[2][1],ch[2][2]]; "\n p_w := ", p_w;
tf,pt := p_w in C; "\n Is p_w in Ew? \n\t", tf;

"\n===============\n";

p := 2^448 - 2^224 - 1;
BF := FiniteField(p);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
d := -39081; lst := [d];
C := Curve(Rd_var, 1 + d*x^2*y^2 - x^2 - y^2); C;
p := [117812161263436946737282484343310064665180535357016373416879082147939404277809514858788439644911793978499419995990477371552926308078495, 19];
P4 := [39082, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018287275];

"\n The given Ed above and its point p_ad := ", Rd_var!p, " correspond to the following Ew and p_w \n";

ch := Ed2Ew_P4_eval(BF,lst,P4,p);	

a1 := ch[1][1]; a2 := ch[1][2]; a3 := ch[1][3]; a4 := ch[1][4]; a6 := ch[1][5]; 
C := Curve(Rw_var, v^2 - u^3 - a2*u^2 - a4*u - a6); C;
p_w := Rw_var![BF!ch[2][1],BF!ch[2][2]]; "\n p_w := ", p_w;
tf,pt := p_w in C; "\n Is p_w in Ew? \n\t", tf;


"\n===============\n";


BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
lst :=  [ 126, -3959, -6367, 401130, -40538689/4 ];
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6 := lst[5];
C := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); C;
p := [-3, 6757/2]; P4 := p; PP4 := [0,6367/2];

"\n The given Ew above and its point p_w := ", Rw_var!p, " correspond to the following Ed and p_ad \n";

ch := Ew2Ed_P4_eval(BF,lst,P4,p);

d := ch[1][1];
C := Curve(Rd_var, 1 + d*x^2*y^2 - x^2 - y^2); C;
p_ad := Rd_var![ch[2][1],ch[2][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in Ed? \n\t", tf;


"\n===============\n";

BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
d:= 4; lst := [d];
C := Curve(Rd_var, 1 + d*x^2*y^2 - x^2 - y^2); C;	
P4 := [-3, -4545/2];
p := [-1,0];

"\n The given Ed above and its point p_ad := ", Rd_var!p, " correspond to the following Ew and p_w \n";

ch := Ed2Ew_P4_eval(BF,lst,P4,p : a1 := 724, a3 := 6705);
a1 := ch[1][1]; a2 := ch[1][2]; a3 := ch[1][3]; a4 := ch[1][4]; a6 := ch[1][5]; 
C := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); C;
p_w := Rw_var![ch[2][1],ch[2][2]]; "\n p_w := ", p_w;
tf,pt := p_w in C; "\n Is p_w in Ew? \n\t", tf;


"\n===============\n";

p := 2^448 - 2^224 - 1;
BF := FiniteField(p);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
lst :=  [ 0, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018287279, 0,1527402724, 0 ];
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6 := lst[5];
C := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); C;
p := [ 161519716510134864566516401752889896523031413486070680062553377595691628481495727263643640821928121317797529880409697191345221781815562,717225809973557004946798282930174671419297352528066850822411768352295942696132428073016719856723894754274572935636178936371882792255446];
P4 := [39082, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018287275];

"\n The given Ew above and its point p_w := ", Rw_var!p, " correspond to the following Ed and p_ad \n";

ch := Ew2Ed_P4_eval(BF,lst,P4,p);

d := ch[1][1];
C := Curve(Rd_var, 1 + d*x^2*y^2 - x^2 - y^2); C;
p_ad := Rd_var![ch[2][1],ch[2][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in Ed? \n\t", tf;



"\n===============\n";


p := 2^448 - 2^224 - 1;
BF := FiniteField(p);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
lst :=  [BF| 746, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018148150, 1012, 1527025248, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018109403 ];
a1 := lst[1]; a2 := lst[2]; a3 := lst[3]; a4 := lst[4]; a6 := lst[5];
C := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); C;
p := [161519716510134864566516401752889896523031413486070680062553377595691628481495727263643640821928121317797529880409697191345221781815562, 71146943933017546680732675918615085327171698583784130572198257969529428770149519018442156569011409483084453413997390170990995681016312];
P4_l := [39082, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498003709183];
PP4_l := [0, 726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018364933];

"\n The given Ew above and its point p_w := ", Rw_var!p, " correspond to the following Ed and p_ad \n";

ch := Ew2Ed_P4_eval(BF,lst,P4_l, p: PP4 := PP4_l);

lst_d := ch[1]; d := lst_d[1];
C := Curve(Rd_var, 1 + d*x^2*y^2 - x^2 - y^2); C;
p_ad := Rd_var![ch[2][1],ch[2][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in Ed? \n\t", tf;

ch := Ed2Ew_P4_eval(BF,lst_d,P4_l,p_ad : PP4 := PP4_l, a1 := a1, a3 := a3);
a1 := ch[1][1]; a2 := ch[1][2]; a3 := ch[1][3]; a4 := ch[1][4]; a6 := ch[1][5]; 
C := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); C;
p_w := Rw_var![ch[2][1],ch[2][2]]; "\n p_w := ", p_w;
tf,pt := p_w in C; "\n Is p_w in Ew? \n\t", tf;

quit;
