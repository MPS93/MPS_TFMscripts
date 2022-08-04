//** Script to implement the change of coordinates between EdC and Ead **//

printf "\n === CHANGE OF COORD. BETWEEN EDWARDS AND TWISTED EDWARDS === \n";

BF := Rationals();

// Twisted Edwards model -Ead-
AD<a,d,sa,sd,x1,y1,x2,y2> := FunctionField(BF,8);
Rxy<x,y> := PolynomialRing(AD,2);
a := sa^2; d := sd^2;
Ead := a*x^2 + y^2 - 1 - d*(x^2)*(y^2); //Ead;
Ead1 := a*x1^2 + y1^2 - 1 - d*(x1^2)*(y1^2); //Ead1;
Ead2 := a*x2^2 + y2^2 - 1 - d*(x2^2)*(y2^2); //Ead2;
Iad := ideal<Rxy|Ead1,Ead2>; QIad := Rxy/Iad;

// Edwards model -EdC-
B<b,sb,u1,v1,u2,v2> := FunctionField(BF,6);
Ruv<u,v> := PolynomialRing(B,2);
b := sb^2;
Eb := u^2 + v^2 - 1 - b*(u^2)*(v^2); //Eb;
Eb1 := u1^2 + v1^2 - 1 - b*(u1^2)*(v1^2); //Eb1;
Eb2 := u2^2 + v2^2 - 1 - b*(u2^2)*(v2^2); //Eb2;
//Iuv := ideal<Ruv|Eb1,Eb2>; QIuv := Ruv/Iuv;


// Map from EdC to Ead
Ed2Ead_b := hom<B -> AD | d/a,sd/sa,sa*x1,y1,sa*x2,y2>;
Ed2Ead := hom<Ruv -> Rxy | Ed2Ead_b,sa*x,y>;

// Map from Ead to EdC
//Ead2E_ad := hom<AD -> B | 1,b,1,sb,u1,v1,u2,v2>;
//Ead2E := hom<Rxy -> Ruv | Ead2E_ad,u,v>;

// Alternative way of defining Ead2E, in case one wants to illustrate the Sqrt(a) inverse.
Ead2E_ad := hom<AD -> B | Denominator(b),Numerator(b),Denominator(sb),Numerator(sb),u1/Denominator(sb),v1,u2/Denominator(sb),v2>;
Ead2E := hom<Rxy -> Ruv | Ead2E_ad, u/Denominator(sb),v>;


// Checking that addition laws are preserved through these maps:

// Addition laws: 
x3 := (x1*y2 + x2*y1)/(1 + d*x1*x2*y1*y2); y3 := (y1*y2 - a*x1*x2)/(1 - d*x1*x2*y1*y2); // Ead
u3 := (u1*v2 + u2*v1)/(1 + b*u1*u2*v1*v2); v3 := (v1*v2 - u1*u2)/(1 - b*u1*u2*v1*v2); // EdC

//  EdC to Ead
hx1 := Ed2Ead(u1);hx2 := Ed2Ead(u2);
hy1 := Ed2Ead(v1);hy2 := Ed2Ead(v2);
ha := 1; hd := Ed2Ead_b(b);
hsa := 1; hsd := Ed2Ead_b(sb);

tf_u3 := Evaluate(x3, [ha,hd,hsa,hsd,hx1,hy1,hx2,hy2]) eq Ed2Ead(u3); 
tf_v3 := Evaluate(y3, [ha,hd,hsa,hsd,hx1,hy1,hx2,hy2]) eq Ed2Ead(v3); 

// Ead to EdC
gu1 := Ead2E(x1); gu2 := Ead2E(x2);
gv1 := Ead2E(y1); gv2 := Ead2E(y2);
gb := Ead2E(d)/Ead2E(a); gsb := Ead2E(sd)/Ead2E(sa);

tf_x3 := Evaluate(u3, [gb,gsb,gu1,gv1,gu2,gv2]) eq Ead2E(x3); 
tf_y3 := Evaluate(v3, [gb,gsb,gu1,gv1,gu2,gv2]) eq Ead2E(y3); 

printf "\n Here we have the maps that relate the Edwards model and the twisted Edwards model. \n";
printf "\n We use b to denote the parameter of the Edwards curve since we can't use d in two places at the same time. The variables sa,sd,sb denote the square roots of a,d and b respectively.\n";
printf "\n The map 'h' to go from Edwards to twisted Edwards is the following: \n \t%o \n\tb   ->   a = %o, d = %o \n\t(u,v) -> (x,y) = (%o,%o)\n",Ed2Ead,Ead2E_ad(a),Ead2E_ad(d),Ead2E(x),Ead2E(y);
printf "\n This map preserves addition: \n\tGiven (x1,y1) + (x2,y2) = (x3,y3), does it hold that \n\t(h(x1),h(y1)) + (h(x2),h(y2)) = (h(x3), h(y3))? %o \n",tf_u3 and tf_v3;
printf "\n And the map 'g' to go from twisted Edwards to Edwards is the following: \n \t%o \n\ta,d -> b = %o \n\t(x,y) -> (u,v) = (%o,%o)\n",Ead2E,Ed2Ead_b(b),Ed2Ead(u),Ed2Ead(v);
printf "\n This map also preserves addition: \n\tGiven (u1,v1) + (u2,v2) = (u3,v3), does it hold that \n\t(g(u1),g(v1)) + (g(u2),g(v2)) = (g(u3), g(v3))? %o \n",tf_u3 and tf_v3;

quit;
