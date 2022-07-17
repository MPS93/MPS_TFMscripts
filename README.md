

Hello and welcome to the repository of the scripts from my master thesis: 'Elliptic curves: various models and their addition laws'.  Have some fun playing around with the scripts and hopefully you can find something helpful for you here.

All '.m' files are ready to be run in \MAGMA and, in most cases, their computational time also allows to run them in MAGMA's online calculator. For the script where the computational time is too long the the output is available in this repository.

--> Link to MAGMA's online calculator: https://magma.maths.usyd.edu.au/calc

For those interested in running some of the scripts in the online calculator, be aware that you will most likely need to fetch some functions from the 'functions.m' file. So at the beginning of each file substitute these lines with the functions they mention:

	   load "../../functions.m";
	   //load Ead_AddL, Ead_AddL_evals, ProjCoord, AffCoord

	   
//~ Notation ~//	   

--> Ew refers to Weierstrass model.
--> Ed refers to Edwards model.
--> Ead refers to twisted Edwards model.


//~ Collection of all the functions ~//


//--> To change points from projective to affine and viceversa <--//


//--> ProjCoord(PR, p) <--//
To change points from affine to projective P1xP1, with input a polynomial ring 'PR' and the point 'p'.


//--> AffCoord(BF,P) <--//
To change points from projective P1xP1 to affine, with input a base field 'BF' and the point 'p'.


//--> bfz(BF,GB) <--//
Function that takes a base field 'BF' and a Groebner basis 'GB' and outputs the list of solutions of the form U = V = ... = 0, computed by brute force.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//--> Ew_AddL_geom(BF,lst,P,Q) <--//
Addition on the Weierstrass model, using the chord and tangent method. Takes as input a base field 'BF', a list 'lst' of the curve parameters and two points 'P','Q' in P2 projective form.


//--> Ead_AddL_geom(BF,lst,P,Q) <--//
Addition on the twisted Edwards model, using the conic method. Takes as input a base field 'BF', a list 'lst' of the curve parameters and two points 'P','Q' in P2 projective form.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//--> Ew_AddL_b2calc(BF,lst,abc_lst) <--//
Regarding the addition on the Weierstrass model, calculates any of the additions with bidegree (2,2). As input takes a base field 'BF', the list 'lst' of curve parameters and the list 'abc_lst' for the a, b, c values to determine which addition law is calculated.


//--> Ew_AddL(BF,lst) <--//
Complete system of addition laws on the Weierstrass model. Takes as input the base field 'BF' and the list 'lst' of curve parameters and outputs the complete system of addition laws of bidegree (2,2).


//--> Ew_AddL_eval(BF,lst,P,Q) <--//
Addition on the Weierstrass model, using the complete system of addition laws of bidegree (2,2). With the input of the base field 'BF', the list 'lst' of curve parameters and the points 'P', 'Q' in P2 projective form, it outputs P+Q

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//--> Ead_AddL(BF,lst : except := false) <--//
Complete system of addition laws on the twisted Edwards model. Takes as input the base field 'BF', the list 'lst' of curve parameters and as optional parameter the boolean label 'except'. It outputs the complete system of addition laws. If 'except' is true, then the output is in terms of the square roots sa, sd of the a, d parameters (useful for some symbolic checks).


//--> Ead_AddL_eval(BF,lst,P,Q) <--//
Addition on the twisted Edwards model, using the complete system of addition laws. With the input of the base field 'BF', the list 'lst' of curve parameters and the points 'P', 'Q' in P1xP1 projective form, it outputs P+Q


//--> Ead_AddL_evals(BF,lst,P,Q) <--//
Same as 'Ead_AddL_eval' but it outputs the result of both addition laws.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//--> Ew2Ed_P4_symb(BF,lst,P4 : PP4 := [0,0]) <--//
Ew to Ed, using the BL birat. eq., outputs the parameter 'd' of the corresponding Edwards model and the symbolic map. As input, the base field 'BF', the list 'lst' of the ai parameters of the Weierstrass model and the point of order 4 'P4'. As optional parameter also the double of P4 'PP4', in case the input is numerical and in complete Weierstrass form.


//--> Ed2Ew_P4_symb(BF,lst,P4 : PP4 := [0,0], aa1 := 0, aa3 := 0) <--//
EdC to Ew, using the BL birat. eq., outputs the corresponding Weierstrass parameters and the map symbolically. As input, the base field 'BF', the list 'lst' of the parameter d of the Edwards model and the point of order 4 'P4'. As optional parameter also the double of P4 'PP4', 'aa1' and 'aa3' in case that we want the output in complete Weierstrass form.


//--> Ew2Ed_symb(BF,lst) <--//
Ew to Ed, using the MD birat. eq., outputs the parameter 'd' of the corresponding Edwards model and the symbolic map. As input, the base field 'BF' and the list 'lst' of the ai parameters of the Weierstrass model.


//--> Ed2Ew_symb(BF,lst) <--//
Ed to Ew, using the BL birat. eq., outputs the corresponding Weierstrass parameters and the map symbolically. As input, the base field 'BF' and the list 'lst' of the parameter d of the Edwards model.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//--> Ew2Ed_P4_eval(BF,lst,P4,p : PP4 := [0,0]) <--//
Ew to Ed, using the BL birat. eq., outputs the Edwards model parameter 'd' and  the coordinates of the corresponding point. As input: the base field 'BF', the list 'lst' of the ai parameters of the Weierstrass model, the point of order 4 'P4' and a point 'p'. As optional parameter also the double of P4 'PP4', in case the input is in complete Weierstrass form.


//--> Ed2Ew_P4_eval(BF,lst,P4,p : PP4 := [0,0], a1 := 0, a3 := 0) <--//
Ed to Ew, using the BL birat. eq., outputs the Weierstrass parameters and the coordinates of the corresponding point. As input, the base field 'BF', the list 'lst' of the parameter d of the Edwards model, the point 'P4' of order 4 and the point 'p'. As optional parameter also the double of P4 'PP4', 'a1' and 'a3' in case that we want the output in complete Weierstrass form.


//--> Ew2Ed_eval(BF,lst,p: PP4 := [0,0]) <--//
Ew to Ed, using the MD birat. eq., outputs the Edwards model parameter 'd' and  the coordinates of the corresponding point. As input: the base field 'BF', the list 'lst' of the ai parameters of the Weierstrass model and a point 'p'. As optional parameter also the double of P4 'PP4', in case that the input is in complete Weierstrass form.


//--> Ed2Ew_eval(BF,lst,p: PP4 := [0,0], a1 := 0, a3 := 0) <--//
Ed to Ew, using the MD birat. eq., outputs the Weierstrass parameters and the coordinates of the corresponding point. As input, the base field 'BF', the list 'lst' of the parameter d of the Edwards model and the point 'p'. As optional parameter also the double of P4 'PP4', 'a1' and 'a3' in case that we want the output in complete Weierstrass form.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//--> Ew2Ead_P4_eval(BF,lst,P4,p : PP4 := [0,0], a := 1) <--//
Ew to Ead, using the BL birat. eq., outputs the twisted Edwards model parameters 'a', 'd' and  the coordinates of the corresponding point. As input: the base field 'BF', the list 'lst' of the ai parameters of the Weierstrass model, the point of order 4 'P4' and a point 'p'. As optional parameter also the double of P4 'PP4', in case the input is in complete Weierstrass form.

//--> Ead2Ew_P4_eval(BF,lst,P4,p : PP4 := [0,0], a1 := 0, a3 := 0) <--//
Ead to Ew, using the BL birat. eq., outputs the Weierstrass parameters and the coordinates of the corresponding point. As input, the base field 'BF', the list 'lst' of the parameters of the twisted Edwards model, the point 'P4' of order 4 and the point 'p'. As optional parameter also the double of P4 'PP4', 'a1' and 'a3' in case that we want the output in complete Weierstrass form.

//--> Ew2Ead_eval(BF,lst,p: PP4 := [0,0], a := 1) <--//
Ew to Ead, using the MD birat. eq., outputs the twisted Edwards model parameters 'a', 'd' and  the coordinates of the corresponding point. As input: the base field 'BF', the list 'lst' of the ai parameters of the Weierstrass model and a point 'p'. As optional parameter also the double of P4 'PP4', in case that the input is in complete Weierstrass form.

//--> Ead2Ew_eval(BF,lst,p: PP4 := [0,0], a1 := 0, a3 := 0) <--//
Ead to Ew, using the MD birat. eq., outputs the Weierstrass parameters and the coordinates of the corresponding point. As input, the base field 'BF', the list 'lst' of the parameters of the twisted Edwards model and the point 'p'. As optional parameter also the double of P4 'PP4', 'a1' and 'a3' in case that we want the output in complete Weierstrass form.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//--> Ew2Ead_iso_eval(BF,lst,p : PP4 := [0,0], u_lst := [0,0,0]) <--//
Ew to Ead, using the 2-isogeny, outputs the twisted Edwards model parameters 'a', 'd' and the coordinates of the corresponding point. As input: the base field 'BF', the list 'lst' of the ai parameters of the Weierstrass model and a point 'p'. As optional parameter also the double of P4 'PP4', in case that the input is in complete Weierstrass form, and the list 'u_lst' of ui values.


//--> Ead2Ew_iso_eval(BF,lst,p : PP4 := [0,0], a1 := 0, a3 := 0, u0 := 0) <--//
Ead to Ew, using the 2-isogeny, outputs the Weierstrass parameters and the coordinates of the corresponding point. As input, the base field 'BF', the list 'lst' of the parameters of the twisted Edwards model and the point 'p'. As optional parameter also the double of P4 'PP4', 'a1' and 'a3' in case that we want the output in complete Weierstrass form and the 'u0'.