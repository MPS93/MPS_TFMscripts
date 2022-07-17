//*** Script to solve the system of equations defined by the Groebner basis of each addition law ***//

load "../../functions.m";
//load Ead_AddL, bfz

printf "\n === CHECKING THE COMPLETENESS OF THE SYSTEM OF ADDITION LAWS === \n";

//procedure to print a list nicely
procedure PrintMyList(lst)

	printf"\n[";
	for i in [1..#lst-1] do
		pair := lst[i];
		printf"[";
		for j in [1..#pair-1] do
			printf "%o,",pair[j];
		end for;
		printf "%o], \n",pair[#pair];
	end for;
	pair := lst[#lst];
	printf"[";
	for j in [1..#pair-1] do
		printf "%o,",pair[j];
	end for;
	printf "%o]] \n",pair[#pair];

end procedure;


//set up to get the formulae to work with them:
BF := Rationals();
R1<a,d> := FunctionField(BF,2);
R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8);
P := [[X1,Z1],[Y1,T1]]; Q := [[X2,Z2],[Y2,T2]];

Ead1 := Z1^2*T1^2 + d*X1^2*Y1^2 - a*X1^2*T1^2 - Y1^2*Z1^2; 
Ead1 := Z1^2*T1^2 + d*X1^2*Y1^2 - a*X1^2*T1^2 - Y1^2*Z1^2; 
Ead2 := Z2^2*T2^2 + d*X2^2*Y2^2 - a*X2^2*T2^2 - Y2^2*Z2^2;
I := ideal< R2 | Ead1, Ead2>;

X3_o := (R2 ! Ead_AddL(BF,[a,d])[1][1][1]); Y3_o := (R2 ! Ead_AddL(BF,[a,d])[1][2][1]);
Z3_o := (R2 ! Ead_AddL(BF,[a,d])[1][1][2]); T3_o := (R2 ! Ead_AddL(BF,[a,d])[1][2][2]);

X3_d := (R2 ! Ead_AddL(BF,[a,d])[2][1][1]); Y3_d := (R2 ! Ead_AddL(BF,[a,d])[2][2][1]);
Z3_d := (R2 ! Ead_AddL(BF,[a,d])[2][1][2]); T3_d := (R2 ! Ead_AddL(BF,[a,d])[2][2][2]);



//--> Solving the system for the original addition law
II_o := ideal<R2 | X3_o,Y3_o,Z3_o,T3_o>;
GB_o := GroebnerBasis(II_o); //GB doesn't work on ideals over QQ; embedding the GB1 into QQ gives the same polys though
//gets all f in GB equal to 0, solve the system and you'll have all the zeroes
printf "\n---> Groebner basis GB_o of the original addition law:\n"; GB_o;

lst_o := bfz(BF,GB_o); 

printf"These are the solutions for GB_o: \n";

PrintMyList(lst_o);

printf"\n\n ==> Note that none of these sets of solutions can happen so there are no pairs of points P1, P2 such that the four polynomials of the original addition law evaluate to zero simultaneously. \n";

//--> Solving the system for the dual addition law
II_d := ideal<R2 | X3_d,Y3_d,Z3_d,T3_d>;
GB_d := GroebnerBasis(II_d); 
printf "\n---> Groebner basis GB_d of the dual addition law:\n"; GB_d;

lst_d := bfz(BF,GB_d); 

printf"These are the solutions for GB_d: \n";

PrintMyList(lst_d);

printf"\n\n ==> Note that none of these sets of solutions can happen so there are no pairs of points P1, P2 such that the four polynomials of the dual addition law evaluate to zero simultaneously. \n";

//--> Solving the system for the complete addition law system
II := ideal<R2 | X3_o,Y3_o,Z3_o,T3_o,X3_d,Y3_d,Z3_d,T3_d>;
GB := GroebnerBasis(II); 
printf "\n---> Groebner basis GB of the complete addition law system:\n"; GB;


lst := bfz(BF,GB); 

printf"These are the solutions for GB: \n";

PrintMyList(lst);

printf"\n\n ==> Note that all these sets of solutions are not valid projective coordinates so there are no pairs of points P1, P2 such that the eight polynomials of the whole addition law system evaluate to zero simultaneously.\n";

printf"\n CONCLUSION: the system of addition laws is complete. \n\n";

quit;
