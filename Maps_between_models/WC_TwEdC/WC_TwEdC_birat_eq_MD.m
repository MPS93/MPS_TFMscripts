//*** Implementing the birat. eq. between Twisted Edwards and Weierstrass by composing maps of Ew2Ead and Ead2Ed ***//
//-> Here we will use the MD birat. eq. <-//

load "../functions.m";
// load ProjCoord, AffCoord, Ew2Ead_eval, Ead2Ew_eval


printf "\n\n MD BIRATIONAL EQUIVALENCE BETWEEN Ead AND Ew \n\n";

//--> Data to begin Ew->Ead example <--//
BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
lst_w :=  [ 126, -3959, -6367, 401130, -40538689/4 ];
a1 := lst_w[1]; a2 := lst_w[2]; a3 := lst_w[3]; a4 := lst_w[4]; a6 := lst_w[5];
C := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); C;
p := [-3, 6757/2]; P4 := p; PP4 := [0,6367/2];

"\n The given Ew above and its point p_w := ", Rw_var!p, " correspond to the following Ead and p_ad: \n";

ch := Ew2Ead_eval(BF,lst_w,p : PP4 := PP4);
lst_ad := ch[2]; d := lst_ad[2];
C := Curve(Rd_var, 1 + d*x^2*y^2 - x^2 - y^2); C;
p_ad := Rd_var![ch[3][1],ch[3][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in Ead? \n\t", tf;

ch2 := Ead2Ew_eval(BF,ch[2],p_ad : PP4 := [0,6367/2], a1 := lst_w[1], a3 := lst_w[3]);
lst_w2 := ch2[2]; p_w2 := ch2[3]; 
"\n Do we get the original Ew if we apply the inverse map to the obtained Ead? \n\t", lst_w eq lst_w2;
"\n Is phi(p_ad) = p_w? \n\t", p eq p_w2;

"\n The above correspondence is made assuming a = 1. If we prefer to work with a twisted Edwards curve that it's not an Ead, we can always specify a different value for 'a'. \n";
"\n So for example, the above Ew is also birat. equivalent to the following Ead and p_ad: \n";

ch := Ew2Ead_eval(BF,lst_w,p : PP4 := PP4, a := 5);
BF := ch[1]; 
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
lst_ad := ch[2]; a:= lst_ad[1]; d := lst_ad[2];
C := Curve(Rd_var, 1 + d*x^2*y^2 - a*x^2 - y^2); C;
p_ad := Rd_var![ch[3][1],ch[3][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in Ead? \n\t", tf;


ch2 := Ead2Ew_eval(BF,ch[2],p_ad : PP4 := [0,6367/2], a1 := lst_w[1], a3 := lst_w[3]);
lst_w2 := ch2[2]; p_w2 := ch2[3]; 
"\n Do we get the original Ew if we apply the inverse map to the obtained Ead? \n\t", lst_w eq lst_w2;
"\n Is phi(p_ad) = p_w? \n\t", p eq p_w2;

"\n===============\n";

// --> Data to begin the FF Ew->Ead example <--//
BF := FiniteField(23);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
a1 := 0; a2 := 12; a3 := 0; a4 := 16; a6 := 0;
lst_w := [a1,a2,a3,a4,a6];
Cw := EllipticCurve([BF|a1,a2,a3,a4,a6]);
u := 1; v := 12; p := [u,v];

"\n The given Ew above and its point p_w := ", Rw_var!p, " correspond to the following Ead and p_ad: \n";

ch := Ew2Ead_eval(BF,lst_w,p);
lst_ad := ch[2]; d := lst_ad[2];
C := Curve(Rd_var, 1 + d*x^2*y^2 - x^2 - y^2); C;
p_ad := Rd_var![ch[3][1],ch[3][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in Ead? \n\t", tf;

ch2 := Ead2Ew_eval(BF,ch[2],p_ad : PP4 := [0,6367/2], a1 := lst_w[1], a3 := lst_w[3]);
lst_w2 := ch2[2]; p_w2 := ch2[3]; 
"\n Do we get the original Ew if we apply the inverse map to the obtained Ead? \n\t", lst_w eq lst_w2;
"\n Is phi(p_ad) = p_w? \n\t", p eq p_w2;

"\n===============\n";

//--> Data to begin Ead->Ew example <--//
BF := Rationals();
R1<a,d> := FunctionField(BF,2);
R2<u,v> := AffineSpace(R1,2);
R3<x,y> := AffineSpace(R1,2);
a := 82; d := 17;
lst_ad := [a,d];
Cad := Curve(R3, -1 - d*x^2*y^2 + a*x^2 + y^2); C;
x := 1; y := 9/4; p := [x,y]; 

"\n The given Ead above and its point p_ad := ", R2!p, " correspond to the following Ew and p_w: \n";

ch := Ead2Ew_eval(BF,lst_ad,p);
BF := ch[1];
"\nBF is ", BF;
R<a1,a2,a3,a4,a6,d> := FunctionField(BF,6);
Rw_var<u,v> := AffineSpace(R,2);
a1 := ch[2][1]; a2 := ch[2][2]; a3 := ch[2][3]; a4 := ch[2][4]; a6 := ch[2][5]; 
C := Curve(Rw_var, v^2 - u^3 - a2*u^2 - a4*u - a6); C;
p_w := Rw_var![BF!ch[3][1],BF!ch[3][2]]; "\n p_w := ", p_w;
tf,pt := p_w in C; "\n Is p_w in Ew? \n\t", tf;


ch2 := Ew2Ead_eval(BF,ch[2],p_w : a := 82);
lst_ad2 := ch2[2]; p_ad2 := ch2[3]; 
"\n Do we get the original Ead if we apply the inverse map to the obtained Ew? \n\t", lst_ad eq lst_ad2;
"\n Is phi_inv(p_w) = p_ad? \n\t", p eq p_ad2;

"\n===============\n";


// --> Data to begin the FF Ead->Ew example <--//
p := 2^382 - 105;
BF := FiniteField(p);
R1<a,d> := FunctionField(BF,2);
R2<u,v> := AffineSpace(R1,2);
R3<x,y> := AffineSpace(R1,2);
a := 1; d := -67254;
lst_ad := [a,d];
Cad := Curve(R3, -1 - d*x^2*y^2 + a*x^2 + y^2); Cad;
p := [BF|3914921414754292646847594472454013487047137431784830634731377862923477302047857640522480241298429278603678181725699,17];
// p := [67255, 9850501549098619803069760025035903451269934817616361666987073351061430442874302652853566563721228910201656997442089]; is of order 4

"\n The given Ead above and its point p_ad := ", R2!p, " correspond to the following Ew and p_w: \n";

ch := Ead2Ew_eval(BF,lst_ad,p);
BF := ch[1];
R<a1,a2,a3,a4,a6,d> := FunctionField(BF,6);
Rw_var<u,v> := AffineSpace(R,2);
a1 := ch[2][1]; a2 := ch[2][2]; a3 := ch[2][3]; a4 := ch[2][4]; a6 := ch[2][5]; 
C := Curve(Rw_var, v^2 - u^3 - a2*u^2 - a4*u - a6); C;
p_w := Rw_var![BF!ch[3][1],BF!ch[3][2]]; "\n p_w := ", p_w;
tf,pt := p_w in C; "\n Is p_w in Ew? \n\t", tf;

ch2 := Ew2Ead_eval(BF,ch[2],p_w : a := a);
lst_ad2 := ch2[2]; p_ad2 := ch2[3]; 
"\n Do we get the original Ead if we apply the inverse map to the obtained Ew? \n\t", lst_ad eq lst_ad2;
"\n Is phi_inv(p_w) = p_ad? \n\t", p eq p_ad2;

quit;
