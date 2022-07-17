//*** About the exceptional points of the 001 addition law for WC ***//
//** complete form *//

BF := Rationals();
R1<a1,a2,a3,a4,a6> := FunctionField(BF,5); 
R2<X1,Y1,Z1, X2,Y2,Z2> := PolynomialRing(R1,6);

// Putting together the addition law for a=b=0, c=1    
mons := [ X1^2*X2^2,X1^2*X2*Y2,X1^2*X2*Z2, X1^2*Y2^2, X1^2*Y2*Z2,
          X1^2*Z2^2, X1*Y1*X2^2, X1*Y1*X2*Y2, X1*Y1*X2*Z2, X1*Y1*Y2^2,
          X1*Y1*Y2*Z2, X1*Y1*Z2^2, X1*Z1*X2^2, X1*Z1*X2*Y2, X1*Z1*X2*Z2,
          X1*Z1*Y2^2, X1*Z1*Y2*Z2, X1*Z1*Z2^2, Y1^2*X2^2, Y1^2*X2*Y2, 
          Y1^2*X2*Z2, Y1^2*Y2^2, Y1^2*Y2*Z2, Y1^2*Z2^2, Y1*Z1*X2^2, 
          Y1*Z1*X2*Y2, Y1*Z1*X2*Z2, Y1*Z1*Y2^2, Y1*Z1*Y2*Z2, Y1*Z1*Z2^2, 
          Z1^2*X2^2, Z1^2*X2*Y2, Z1^2*X2*Z2, Z1^2*Y2^2, Z1^2*Y2*Z2, Z1^2*Z2^2];
  
X3_001_ci := [ <-a2, 3>, <a1, 5>, <-a4, 6>, <2, 11>, <a3, 12>, <a2, 13>, 
               <1, 16>, <2*a3, 17>, <-3*a6, 18>, <-1, 21>, <-a1, 25>, 
               <-2, 26>, <-2*a3, 27>, <a4, 31>, <-a3, 32>, <3*a6, 33>];
Y3_001_ci := [ <-3, 2>, <a1*a2 - 3*a3, 3>, <-a1^2 - a2, 5>, 
               <a1*a4 - a2*a3, 6>, <3, 7>, <2*a2, 9>, <-2*a1, 11>, 
               <a4, 12>, <-a1*a2 + 3*a3, 13>, <-2*a2, 14>,
               <-2*a1*a3 - 2*a4, 17>, <3*a1*a6 - a3*a4, 18>, <-1, 23>, 
               <a1^2 + a2, 25>, <2*a1, 26>, <2*a1*a3 + 2*a4, 27>, <1, 28>, 
               <a3^2 + 3*a6, 30>, <-a1*a4 + a2*a3, 31>, <-a4, 32>, 
               <-3*a1*a6 + a3*a4, 33>, <-a3^2 - 3*a6, 35>];

Z3_001_ci := [ <3, 3>, <a2, 6>, <-a1, 12>, <-3, 13>, <a4, 18>, <-1, 24>, 
               <-a3, 30>, <-a2, 31>, <a1, 32>, <-a4, 33>, <1, 34>, <a3, 35>];
  
 X3_001 := R2!&+[R2!(X3_001_ci[i][1])*R2!(mons[Integers()!X3_001_ci[i][2]]): i in [1..#X3_001_ci]];
 Y3_001 := R2!&+[R2!(Y3_001_ci[i][1])*R2!(mons[Integers()!Y3_001_ci[i][2]]): i in [1..#Y3_001_ci]];
 Z3_001 := R2!&+[R2!(Z3_001_ci[i][1])*R2!(mons[Integers()!Z3_001_ci[i][2]]): i in [1..#Z3_001_ci]];

x3 := X3_001; y3 := Y3_001; z3 := Z3_001;

//--> set up to do checks <--//
WC1 := Y1^2*Z1+a1*X1*Y1*Z1+a3*Y1*Z1^2-X1^3-a2*X1^2*Z1-a4*X1*Z1^2-a6*Z1^3; 
WC2 := Y2^2*Z2+a1*X2*Y2*Z2+a3*Y2*Z2^2-X2^3-a2*X2^2*Z2-a4*X2*Z2^2-a6*Z2^3;
wc := { WC1, WC2 };

I0 := ideal< R2 | wc >;
// QQ<x1,y1,z1,x2,y2,z2> := quo<R2 | I0 >;

I3 := ideal<  R2 | wc join {x3, y3, z3} >; 
time GB3 := GroebnerBasis(I3); #GB3;

xy := X1*Y2-X2*Y1;
xz := X1*Z2-X2*Z1;
yz := Y1*Z2-Y2*Z1;
xyz := { xy, xz, yz };

I10 := ideal< R2 | wc join xyz >;
time GB10 := GroebnerBasis(I10); #GB10;

//--> Checks <--//
I3 subset I10;
//  I10 subset I3; 
for k := 1 to #GB10 do
erin := GB10[k] in I3;
k, erin;
// if not erin then break; end if;
end for;

xyz2 := {a*b : a, b in xyz};
I10b := ideal< R2 | wc join xyz2 >;
time GB10b := GroebnerBasis(I10b); #GB10b;
I3 subset I10b;
//  I10b subset I3; 
for k := 1 to #GB10b do
erin := GB10b[k] in I3;
k, erin;
// if not erin then break; end if;
end for;

/*
xyz3 := {a*b : a in xyz, b in xyz2};
I10c := ideal< R2 | wc join xyz3 >;
time GB10c := GroebnerBasis(I10c); #GB10c;
I3 subset I10c;
I10c subset I3; 
for k := 1 to #GB10c do
erin := GB10c[k] in I3;
k, erin;
// if not erin then break; end if;
end for;
*/

print "Gelijk door naar xyz4";

xyz4 := {a*b : a, b in xyz2};
I10d := ideal< R2 | wc join xyz4 >;
SetVerbose("Groebner", true);
time GB10d := GroebnerBasis(I10d); #GB10d;
I3 subset I10d;
//  I10d subset I3; 
for k := 1 to #GB10d do
erin := GB10d[k] in I3;
k, erin;
// if not erin then break; end if;
end for;

