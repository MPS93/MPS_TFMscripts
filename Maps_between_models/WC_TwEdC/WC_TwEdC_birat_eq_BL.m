//*** Implementing the birat. eq. between Twisted Edwards and Weierstrass by composing maps of Ew2EdC and Ead2Ed ***//
//-> Here we will use the BL birat. eq. <-//

load "../functions.m";
// load ProjCoord, AffCoord, Ew2Ead_P4_eval, Ead2Ew_P4_eval


printf "\n\n BL BIRATIONAL EQUIVALENCE BETWEEN Ead AND Ew \n\n";

//--> Data to begin the Ead-> Ew example <--//
BF := Rationals();
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
lst_w :=  [ 126, -3959, -6367, 401130, -40538689/4 ];
a1 := lst_w[1]; a2 := lst_w[2]; a3 := lst_w[3]; a4 := lst_w[4]; a6 := lst_w[5]; 
C := Curve(Rw_var, v^2 - u^3 - a2*u^2 - a4*u - a6); C;
P4 := [-3, 6757/2]; PP4 := [0,6367/2];
p := [-3, 6733/2];

"\n The given Ew above and its point p_w := ", Rw_var!p, " correspond to the following TwEC and p_ad: \n";

ch := Ew2Ead_P4_eval(BF,lst_w,P4,p : PP4 := [0,6367/2], a := 16);

a := ch[2][1]; d := ch[2][2];
C := Curve(Rd_var, 1 + d*x^2*y^2 - a*x^2 - y^2); C;
p_ad := Rd_var![ch[3][1],ch[3][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in the Ead? \n\t", tf;

ch2 := Ead2Ew_P4_eval(BF,ch[2],P4,p_ad : PP4 := [0,6367/2], a1 := lst_w[1], a3 := lst_w[3]);
lst_w2 := ch2[2]; p_w2 := ch2[3]; 
"\n Do we get the original Ew if we apply the inverse map to the obtained Ead? \n\t", lst_w eq lst_w2;
"\n Is phi(p_ad) = p_w? \n\t", p eq p_w2;

"\n The above correspondence is made with a = 16, which is a square over the rationals. But we also have the chance to use a non-square 'a': \n";

ch := Ew2Ead_P4_eval(BF,lst_w,P4,p : PP4 := [0,6367/2], a := 29);

BF := ch[1];
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
a := ch[2][1]; d := ch[2][2];
C := Curve(Rd_var, 1 + d*x^2*y^2 - a*x^2 - y^2); C;
p_ad := Rd_var![ch[3][1],ch[3][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in the TwEC? \n\t", tf;

ch2 := Ead2Ew_P4_eval(BF,ch[2],P4,p_ad : PP4 := [0,6367/2], a1 := lst_w[1], a3 := lst_w[3]);
lst_w2 := ch2[2]; p_w2 := ch2[3]; 
"\n Do we get the original Ew if we apply the inverse map to the obtained Ead? \n\t", lst_w eq lst_w2;
"\n Is phi(p_ad) = p_w? \n\t", p eq p_w2;

"\n===============\n";
//--> Data to begin the FF Ead-> Ew example <--//

BF := FiniteField(23);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rw_var<u,v> := AffineSpace(R,2);
Rd_var<x,y> := AffineSpace(R,2);
a1 := 0; a2 := 12; a3 := 0; a4 := 16; a6 := 0;
lst_w := [a1,a2,a3,a4,a6];
C := Curve(Rw_var, v^2 - u^3 - a2*u^2 - a4*u - a6); C;
P4 := [19,8]; PP4 := [0,0];  //[19,15]
p := [15,17]; //[1,12]

"\n The given Ew above and its point p_w := ", Rw_var!p, " correspond to the following TwEC and p_ad: \n";

ch := Ew2Ead_P4_eval(BF,lst_w,P4,p : a := 13);

a := ch[2][1]; d := ch[2][2];
C := Curve(Rd_var, 1 + d*x^2*y^2 - a*x^2 - y^2); C;
p_ad := Rd_var![ch[3][1],ch[3][2]]; "\n p_ad := ", p_ad;
tf,pt := p_ad in C; "\n Is p_ad in the TwEC? \n\t", tf;

ch2 := Ead2Ew_P4_eval(BF,ch[2],P4,p_ad : PP4 := [0,6367/2], a1 := lst_w[1], a3 := lst_w[3]);
lst_w2 := ch2[2]; p_w2 := ch2[3]; 
"\n Do we get the original Ew if we apply the inverse map to the obtained Ead? \n\t", lst_w eq lst_w2;
"\n Is phi(p_ad) = p_w? \n\t", p eq p_w2;

"\n===============\n";
//--> Data to begin the Ead-> Ew example <--//
BF := AlgebraicClosure(Rationals());
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
tf, sqrt_a := IsSquare(BF!37);
a := 37; d:= 148; lst_ad := [a,d];
C := Curve(Rd_var, 1 + d*x^2*y^2 - a*x^2 - y^2); C;	
p := [ -4/109187*sqrt_a,  -1]; 
P4 := [-3, -4545/2];

"\n The given TwEC above and its point p_ad := ", Rd_var!p, " correspond to the following Ew and p_w: \n";

ch := Ead2Ew_P4_eval(BF,lst_ad,P4,p: a1 := 2645, a3 := 705/17);
BF := ch[1];
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
a1 := ch[2][1]; a2 := ch[2][2]; a3 := ch[2][3]; a4 := ch[2][4]; a6 := ch[2][5];
C := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); C;
p_w := Rw_var![ch[3][1],ch[3][2]]; "\n p_w := ", p_w;
tf,pt := p_w in C; "\n Is p_w in the Ew? \n\t", tf;

"\n Inverse check can't be done here since the inverse is not defined for this p_w point.";

"\n===============\n";
//--> Data to begin the Ead-> Ew example <--//
p := 2^382 - 105;
BF := FiniteField(p);
R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
Rd_var<x,y> := AffineSpace(R,2);
Rw_var<u,v> := AffineSpace(R,2);
a := 18935; d := -1273454490;
lst_ad := [a,d];
C := Curve(Rd_var, 1 + d*x^2*y^2 - a*x^2 - y^2); C;	
P4 := [67255, 9850501549098619803069760025035903451269934817616361666987073351061430442874302652853566563721228910201656997442089];
p := [BF|2489759715079200041202855701086126438986119374073451597835957207123822721556226164156606228525545369658577863518126,17];

"\n The given TwEC above and its point p_ad := ", Rd_var!p, " correspond to the following Ew and p_w: \n";

ch := Ead2Ew_P4_eval(BF,lst_ad,P4,p);

a1 := ch[2][1]; a2 := ch[2][2]; a3 := ch[2][3]; a4 := ch[2][4]; a6 := ch[2][5]; 
C := Curve(Rw_var, v^2 + a1*u*v + a3*v - u^3 - a2*u^2 - a4*u - a6); C;
p_w := Rw_var![ch[3][1],ch[3][2]]; "\n p_w := ", p_w;
tf,pt := p_w in C; "\n Is p_w in the Ew? \n\t", tf;

ch2 := Ew2Ead_P4_eval(BF,ch[2],P4,p_w : a := 18935);
lst_ad2 := ch2[2]; p_ad2 := ch2[3]; 
"\n Do we get the original Ead if we apply the inverse map to the obtained Ew? \n\t", lst_ad eq lst_ad2;
"\n Is phi_inv(p_w) = p_ad? \n\t", p eq p_ad2;

quit;
