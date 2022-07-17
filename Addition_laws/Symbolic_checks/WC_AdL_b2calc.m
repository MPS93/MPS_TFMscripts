//*** Function that outputs bidegree (2,2) add laws given any abc values ***//

load "../../functions.m";
//load Ew_AddL_b2calc


BF := Rationals();
AAA<a1,a2,a3,a4,a6,	A880880, A880110, A880101, A880011, A880088, A880002, A110880, A110110, A110101, A110011, A110088, A110002, A101880, A101110, A101101, A101011, A101088, A101002, A011880, A011110, A011101, A011011, A011088, A011002, A088880, A088110, A088101, A088011, A088088, A088002, A002880, A002110, A002101, A008811, A008888, A008802 > :=  FunctionField(BF, 41);
RRR<X1,Y1,Z1, X2,Y2,Z2> := PolynomialRing(AAA,6);	  

AdL_001 := Ew_AddL_b2calc(Rationals(),[a1,a2,a3,a4,a6],[0,0,1]);
AdL_010 := Ew_AddL_b2calc(Rationals(),[a1,a2,a3,a4,a6],[0,1,0]);

printf"\n X3_001 := %o; \n", AdL_001[1][1];
printf"\n Y3_001 := %o; \n", AdL_001[1][2];
printf"\n Z3_001 := %o; \n", AdL_001[1][3];

printf"\n X3_010 := %o; \n", AdL_010[1][1];
printf"\n Y3_010 := %o; \n", AdL_010[1][2];
printf"\n Z3_010 := %o; \n", AdL_010[1][3];

//Ew_AddL_b2calc(Rationals(),[a1,a2,a3,a4,a6],[3,7,-8]);
//&+[R2!(X3_001_ci[i][1])*R2!(mons[Integers()!X3_001_ci[i][2]]): i in [1..#X3_001_ci]]


quit;
