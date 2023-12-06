//*** Collection of all the functions ***//


//--> To change points from projective to affine and viceversa <--//
ProjCoord := function(PR, p)
	return [[PR ! p[1], PR ! 1],[PR ! p[2],PR ! 1]];
end function;

AffCoord := function(BF,P)
	if (P[1][2] eq 0) or (P[2][2] eq 0) then
		return "Error: point at infinity"; 
	else
		return [BF ! (P[1][1] div P[1][2]), BF ! (P[2][1]div P[2][2])];
	end if;
end function;


//--> function that outputs the list of solutions of the form U = V = ... = 0, computed by brute force: <--//
bfz := function(BF,GB)
	
	R1<a,d> := FunctionField(BF,2);
    R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8);

	GB := [R2!GB[i]:i in [1..#GB]];
	nGB := #GB;
	S8 := Subsets({1,2,3,4,5,6,7,8});
	S8 := Setseq(S8);
	sols := [];
	for idx_list in S8 do
		listvar := [X1, Y1, Z1, T1, X2, Y2, Z2, T2];
		for idx in idx_list do
		    listvar[idx] := 0;
		end for;
		zrs := [];
		for f in GB do
		    if Evaluate(f,listvar) eq 0 then
		        Append(~zrs,0);
		    end if;
		end for;
		if #zrs eq nGB then Append(~sols,Setseq(idx_list));end if;
		listvar := [X1, Y1, Z1, T1, X2, Y2, Z2, T2];
	end for;
	sols1 := [];sols2 := [];sols3 := [];sols4 := [];sols5 := [];sols6 := [];sols7 := [];sols8 := [];
	for s in sols do
		if #s eq 1 then Append(~sols1,s); end if;
		if #s eq 2 then Append(~sols2,s); end if;
		if #s eq 3 then Append(~sols3,s); end if;
		if #s eq 4 then Append(~sols4,s); end if;
		if #s eq 5 then Append(~sols5,s); end if;
		if #s eq 6 then Append(~sols6,s); end if;
		if #s eq 7 then Append(~sols7,s); end if;
		if #s eq 8 then Append(~sols8,s); end if;
	end for;
	for s2 in sols2 do
		for s3 in sols3 do
		    if s2[1] in s3 and s2[2] in s3 then Exclude(~sols3,s3);end if;
		end for;
		for s4 in sols4 do
		    if s2[1] in s4 and s2[2] in s4 then Exclude(~sols4,s4);end if;
		end for;
		for s5 in sols5 do
		    if s2[1] in s5 and s2[2] in s5 then Exclude(~sols5,s5);end if;
		end for;
		for s6 in sols6 do
		    if s2[1] in s6 and s2[2] in s6 then Exclude(~sols6,s6);end if;
		end for;
		for s7 in sols7 do
		    if s2[1] in s7 and s2[2] in s7 then Exclude(~sols7,s7);end if;
		end for;
		for s8 in sols8 do
		    if s2[1] in s8 and s2[2] in s8 then Exclude(~sols8,s8);end if;
		end for;
	end for;
	for s3 in sols3 do
		for s4 in sols4 do
		    if s3[1] in s4 and s3[2] in s4 and s3[3] in s4 then Exclude(~sols4,s4);end if;
		end for;
	end for;
	sols := sols2 cat sols3; Sort(~sols);

	var_sols := [];
	for s in sols do
		var_s := [];
		for idx in s do
			Append(~var_s,listvar[idx]);
		end for;
		Append(~var_sols,var_s);
	end for;
	return Reverse(var_sols);
end function;


//--> WC addition, using chord and tangent method <--//

WC_AddL_geom := function(BF,lst,P,Q)
	
	R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
	R2<X,Y,Z> := ProjectiveSpace(R1,2);
	O := R2![0,1,0]; 
	
	//Get the coeffs
	if #lst eq 2 then
        a1 := R1 ! 0; a2 := R1 ! 0; a3 := R1 ! 0; 
        a4 := R1 ! lst[1]; a6 := R1 ! lst[2];
    elif #lst eq 5 then 
        a1 := R1 ! lst[1]; a2 := R1 ! lst[2]; a3 := R1 ! lst[3]; 
        a4 := R1 ! lst[4]; a6 := R1 ! lst[5];
    else
        return "Error: wrong number of coeffs";
    end if;
	
    // Get the coordinates
	Xp := R1 ! P[1]; Yp := R1 ! P[2]; Zp := R1 ! P[3]; P := R2![Xp,Yp,Zp];
    Xq := R1 ! Q[1]; Yq := R1 ! Q[2]; Zq := R1 ! Q[3]; Q := R2![Xq,Yq,Zq];

	// Define the WC 
	C := Curve(R2, (Y^2)*Z + a1*X*Y*Z + a3*Y*Z^2 - X^3-a2*(X^2)*Z-a4*X*Z^2 - a6*Z^3);
	
	// Compute the coeffs of the line and define it
	// We define the line as a curve because it makes the intersection easier
	if (P eq Q) and (Type(Rationals()) eq Type(BF)) and (a1 eq 0) and (a3 eq 0) then 
	
		L := Curve(R2,X*Zp - Xp*Z);
	
	elif P eq Q then
		l := (3*Xp^2+2*a2*Xp*Zp+a4*Zp^2-a1*Yp*Zp)/(2*Yp*Zp+a1*Xp*Zp+a3*Zp^2); n := (-Xp^3+a4*Xp*Zp^2+2*a6*Zp^3-a3*Yp*Zp^2)/(2*Yp*Zp^2+a1*Xp*Zp^2+a3*Zp^3);

		L := Curve(R2,Y - l*X - n*Z);

	elif (Xp eq Xq) and (Zp eq Zq) and ((Yp eq - Yq-a1*Xq-a3*Zq) or (Yq eq - Yp-a1*Xp-a3*Zp)) then

		L := Curve(R2,X*Zp - Xp*Z);

	elif P eq O then
	
		L := Curve(R2,X*Zq - Xq*Z);
	
	elif Q eq O then

		L := Curve(R2,X*Zp - Xp*Z);
	
	else
		l := (Zp*Yq - Yp*Zq)/(Xq*Zp-Xp*Zq); n := (Yp*Xq - Yq*Xp)/(Xq*Zp-Xp*Zq);

		L := Curve(R2,Y - l*X - n*Z);
	end if;

	// Get the points in the intersection
	pts := Points(Intersection(C,L)); 
	
	// This should only contain R
	if P eq Q then
		R_lst := pts diff Seqset([P,Q,O]);
	else
		R_lst := pts diff Seqset([P,Q]);
	end if;
	
	// In case we have a special situation, we need this.
	if #R_lst ne 1 then
		nmrs := [IntersectionNumber(L,C,pts[i]) : i in [1..#pts]];
		i := Index(nmrs,2); 
		if i ne 0 then
			R := pts[i];
		else
			return "Error: something went wrong in the intersection.";
		end if;
	else
		R := R_lst[1];
	end if;
	
	if R eq O then
		return R;
	else 
		Xr := R[1]; Zr := R[3];
		Lp := Curve(R2,X*Zr - Xr*Z);
		pts := Points(Intersection(C,Lp)); 
		Rp_lst := pts diff Seqset([R,O]);
		if #Rp_lst ne 1 then
			nmrs := [IntersectionNumber(Lp,C,pts[i]) : i in [1..#pts]];
			i := Index(nmrs,2); 
			if i ne 0 then
				Rp := pts[i];
			else
				return "Error: something went wrong in the intersection.";
			end if;
		else
			Rp := Rp_lst[1];
		end if;		
	end if;

	return Rp;

end function;



//--> TwEdC addition using the conic method <--//

TwEd_AddL_geom := function(BF,lst,P,Q)
	
	R1<a,d> := FunctionField(BF,2);
	R2<X,Y,Z> := ProjectiveSpace(R1,2);
	Op := R2 ! [0,-1,1]; O1 := R2![1,0,0]; O2 := R2![0,1,0]; 
	
	// Get the coeffs
    if #lst eq 2 then
        a := R1!lst[1]; d := R1!lst[2];
    else
        return "Error: wrong number of coeffs";
    end if;

    // Get the coordinates
	Xp := R1 ! P[1]; Yp := R1 ! P[2]; Zp := R1 ! P[3]; P := R2![Xp,Yp,Zp];
    Xq := R1 ! Q[1]; Yq := R1 ! Q[2]; Zq := R1 ! Q[3]; Q := R2![Xq,Yq,Zq];

	// Define the TwEdCurve
	Crv := Curve(R2, (a*(X^2) + (Y^2))*Z^2 - Z^4 - d*(X^2)*(Y^2));
	
	// Compute the coeffs of the conic and define it
	// We define the conic as a curve because in some cases its polynomial has a degree too low for the magma command Conic to work
	if P eq Op then
		cz := -Xq; cxy := Zq; cxz := Zq;

		Cnc := Curve(R2,cz*(Z^2 + Y*Z)+cxy*X*Y+cxz*X*Z);

	elif Q eq Op then
		cz := -Xp; cxy := Zp; cxz := Zp;

		Cnc := Curve(R2,cz*(Z^2 + Y*Z)+cxy*X*Y+cxz*X*Z);		

	elif P eq Q then
		cz := Xp*Zp*(Zp-Yp);
		cxy := d*(Xp^2)*Yp-Zp^3;
		cxz := Zp*(Zp*Yp - a*Xp^2);

		Cnc := Curve(R2,cz*(Z^2 + Y*Z)+cxy*X*Y+cxz*X*Z);
	
	else
		cz := Xp*Xq*(Yp*Zq-Yq*Zp);
		cxy := Zp*Zq*(Xp*Zq-Xq*Zp + Xp*Yq - Xq*Yp);
		cxz := Xq*Yq*(Zp^2)-Xp*Yp*(Zq^2)+Yp*Yq*(Xq*Zp-Xp*Zq);

		Cnc := Curve(R2,cz*(Z^2 + Y*Z)+cxy*X*Y+cxz*X*Z);
	end if;

	// Get the points in the intersection
	pts := Points(Intersection(Crv,Cnc)); 
	
	// This should only contain -P3
	P3_lst := pts diff Seqset([P,Q,Op,O1,O2]); 
	
	// In case we have a special situation, we need this.
	if #P3_lst ne 1 then
		og_pts := [P,Q,Op,O1,O2]; nmrs := [1,1,1,2,2];
		if #pts eq #Seqset(og_pts) then
			if P eq Op then nmrs[3] +:= 1; end if;
			if Q eq Op then nmrs[3] +:= 1; end if;
			for i in [1..2] do
				if IntersectionNumber(Cnc,Crv,og_pts[i]) gt nmrs[i] then
					return [-og_pts[i][1],og_pts[i][2],og_pts[i][3]];
				end if;
			end for;
			for i in [3..5] do
				if IntersectionNumber(Cnc,Crv,og_pts[i]) gt nmrs[i] then
					return og_pts[i];
				end if;
			end for;
		return "Error: something went wrong in the intersection.";
		end if;
	end if;
	
	// Note that the output is the opposite of the 8th intersection point.
	P3 := [-P3_lst[1][1], P3_lst[1][2], P3_lst[1][3]];
	
	return P3;

end function;


//--> WC addition, calculates any of the additions with bidegree (2,2) <--//

WC_AddL_b2calc := function(BF,lst,abc_lst)
	
	//abc_lst --> case 1: [0,0,1]; case 2: [0,1,0]; case 3: [1,0,0];

	AAA<a1,a2,a3,a4,a6,	A200200, A200110, A200101, A200011, A200020, A200002,
	A110200, A110110, A110101, A110011, A110020, A110002, A101200, A101110, A101101, A101011, A101020, A101002,
	A011200, A011110, A011101, A011011, A011020, A011002, A020200, A020110, A020101, A020011, A020020, A020002,
	A002200, A002110, A002101, A002011, A002020, A002002 > :=   FunctionField(BF, 41);
	  
    if #lst eq 2 then
        a1 := AAA ! 0; a2 := AAA ! 0; a3 := AAA ! 0; 
        a4 := AAA ! lst[1]; a6 := AAA ! lst[2];
    elif #lst eq 5 then
        a1 := AAA ! lst[1]; a2 := AAA ! lst[2]; a3 := AAA ! lst[3]; 
        a4 := AAA ! lst[4]; a6 := AAA ! lst[5];
    else
        return "Error: wrong number of coeffs";
    end if;

	RRR<X1,Y1,Z1, X2,Y2,Z2> := PolynomialRing(AAA,6);
	I := ideal< RRR |
		   Y1^2*Z1+a1*X1*Y1*Z1+a3*Y1*Z1^2-X1^3-a2*X1^2*Z1-a4*X1*Z1^2-a6*Z1^3,
		   Y2^2*Z2+a1*X2*Y2*Z2+a3*Y2*Z2^2-X2^3-a2*X2^2*Z2-a4*X2*Z2^2-a6*Z2^3 >;
	Q<x1,y1,z1,x2,y2,z2> := quo<RRR | I>;

	F := A200200*X1^2*X2^2+ A200110*X1^2*X2*Y2+A200101*X1^2*X2*Z2+A200011*X1^2*Y2*Z2+A200020*X1^2*Y2^2+A200002*X1^2*Z2^2+A110200*X1*Y1*X2^2+
		A110110*X1*Y1*X2*Y2+A110101*X1*Y1*X2*Z2+A110011*X1*Y1*Y2*Z2+A110020*X1*Y1*Y2^2+A110002*X1*Y1*Z2^2+A101200*X1*Z1*X2^2+A101110*X1*Z1*X2*Y2+
		A101101*X1*Z1*X2*Z2+A101011*X1*Z1*Y2*Z2+A101020*X1*Z1*Y2^2+A101002*X1*Z1*Z2^2+A011200*Y1*Z1*X2^2+A011110*Y1*Z1*X2*Y2+A011101*Y1*Z1*X2*Z2+
		A011011*Y1*Z1*Y2*Z2+A011020*Y1*Z1*Y2^2+A011002*Y1*Z1*Z2^2+A020200*Y1^2*X2^2+A020110*Y1^2*X2*Y2+A020101*Y1^2*X2*Z2+A020011*Y1^2*Y2*Z2+A020020*Y1^2*Y2^2+
		A020002*Y1^2*Z2^2+A002200*Z1^2*X2^2+A002110*Z1^2*X2*Y2+A002101*Z1^2*X2*Z2+A002011*Z1^2*Y2*Z2+A002020*Z1^2*Y2^2+A002002*Z1^2*Z2^2;

        monF := Monomials(F);
		
	B<b1,b2,b3,b4,b6> := FunctionField(BF,5);
	HHH<
	H200200, H200110, H200101, H200011, H200020, H200002,
	H110200, H110110, H110101, H110011, H110020, H110002,
	H101200, H101110, H101101, H101011, H101020, H101002,
	H011200, H011110, H011101, H011011, H011020, H011002,
	H020200, H020110, H020101, H020011, H020020, H020002,
	H002200, H002110, H002101, H002011, H002020, H002002> := PolynomialRing(B, 36);

	h := hom< AAA->HHH | b1,b2,b3,b4,b6,
	H200200, H200110, H200101, H200011, H200020, H200002,
	H110200, H110110, H110101, H110011, H110020, H110002,
	H101200, H101110, H101101, H101011, H101020, H101002,
	H011200, H011110, H011101, H011011, H011020, H011002,
	H020200, H020110, H020101, H020011, H020020, H020002,
	H002200, H002110, H002101, H002011, H002020, H002002 >;

	binv := hom<B->AAA |a1,a2,a3,a4,a6>;
	hinv := hom< HHH->RRR | binv,
	X1^2*X2^2, X1^2*X2*Y2, X1^2*X2*Z2, X1^2*Y2*Z2, X1^2*Y2^2, X1^2*Z2^2,
	X1*Y1*X2^2, X1*Y1*X2*Y2, X1*Y1*X2*Z2, X1*Y1*Y2*Z2, X1*Y1*Y2^2, X1*Y1*Z2^2,
	X1*Z1*X2^2, X1*Z1*X2*Y2, X1*Z1*X2*Z2, X1*Z1*Y2*Z2, X1*Z1*Y2^2, X1*Z1*Z2^2,
	Y1*Z1*X2^2, Y1*Z1*X2*Y2, Y1*Z1*X2*Z2, Y1*Z1*Y2*Z2, Y1*Z1*Y2^2, Y1*Z1*Z2^2,
	Y1^2*X2^2, Y1^2*X2*Y2, Y1^2*X2*Z2, Y1^2*Y2*Z2, Y1^2*Y2^2, Y1^2*Z2^2,
	Z1^2*X2^2, Z1^2*X2*Y2, Z1^2*X2*Z2, Z1^2*Y2*Z2, Z1^2*Y2^2, Z1^2*Z2^2>;

	V := VectorSpace(B, 36);

	kappa := (Y1*Z2+Y2*Z1+a1*X2*Z1+a3*Z1*Z2)/(X1*Z2-X2*Z1);
	mu := -(Y1*X2+Y2*X1+a1*X1*X2+a3*X1*Z2)/(X1*Z2-X2*Z1);  
	sXZ := kappa^2+a1*kappa-(X1*Z2+X2*Z1)/(Z1*Z2)-a2;
	sYZ := -(kappa+a1)*sXZ-mu-a3;
	lambda := (Y1*Z2-Y2*Z1)/(X1*Z2-X2*Z1);
	nu := -(Y1*X2-Y2*X1)/(X1*Z2-X2*Z1);
	f := lambda^2+a1*lambda-(X1*Z2+X2*Z1)/(Z1*Z2)-a2;
	g := -(lambda+a1)*f-nu-a3;
	Z0 := (X1*Z2-Z1*X2)^3/(Z1*Z2);
	a := RRR!abc_lst[1]; b := RRR!abc_lst[2]; c := RRR!abc_lst[3];

	ABC := ( RRR ! Q ! (RRR ! (a*Numerator(sXZ*Z0) + b*Numerator(sYZ*Z0) + c*Z1*Z2*Numerator(Z0))) div (Z1^2*Z2^2));
	Z3 := -ABC;

	coefZ3 := Coefficients(Z3);
	monZ3 := Monomials(Z3);
	Z3_ci := [<coefZ3[i],Index(monF,monZ3[i])>:i in [1..#coefZ3]];

	Nx3 := (RRR ! Q ! (Numerator(f)*ABC)) div (Z1*Z2);
	Fx3 := F*(X1*Z2-Z1*X2)^2 + Nx3;
	Qx3 := Q ! Fx3;
	Cx3 := Coefficients(Qx3);
	Ex3 := [ [B | Coefficient(h(Cx3[i]),j,1) :j in [1..36]] : i in [1..88] ];
	Mx3 := Matrix(B,88,36, Ex3);

	Wx3 := [ h(Cx3[i])-&+[Ex3[i][j]*HHH.j : j in [1..36]] : i in [1..88]];
	Vx3 := Vector(B, 88, Wx3);

	Sx3, Kx3 := Solution(Transpose(Mx3), Vx3);
	X3 := -hinv(&+[ Sx3[i]*HHH.i : i in [1..36]]);

	coefX3 := Coefficients(X3);
	monX3 := Monomials(X3);
	X3_ci := [<coefX3[i],Index(monF,monX3[i])>:i in [1..#coefX3]];

	Ny3 := RRR ! Q ! (Numerator(g)*ABC) div (Z1*Z2);
	Fy3 := F*(X1*Z2-Z1*X2)^3 + Ny3;
	Qy3 := Q ! Fy3;
	Cy3 := Coefficients(Qy3);

	Ey3 := [ [B | Coefficient(h(Cy3[i]),j,1) :j in [1..36]] : i in [1..88] ];
	My3 := Matrix(B,88,36, Ey3);

	Wy3 := [ h(Cy3[i])-&+[Ey3[i][j]*HHH.j : j in [1..36]] : i in [1..88]];
	Vy3 := Vector(B, 88, Wy3);

	Sy3, Ky3 := Solution(Transpose(My3), Vy3);
	Y3 := -hinv(&+[ Sy3[i]*HHH.i : i in [1..36]]);

	coefY3 := Coefficients(Y3);
	monY3 := Monomials(Y3);
	Y3_ci := [<coefY3[i],Index(monF,monY3[i])>:i in [1..#coefY3]];

	return [*[X3,Y3,Z3],[X3_ci,Y3_ci,Z3_ci], monF*];
end function;

//--> WC complete system of addition laws, outputs the system of add. laws <--//

WC_AddL := function(BF,lst)
    R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
    R2<X1,Y1,Z1, X2,Y2,Z2> := PolynomialRing(R1,6);
	
    if #lst eq 2 then
        a1 := R1 ! 0; a2 := R1 ! 0; a3 := R1 ! 0; 
        a4 := R1 ! lst[1]; a6 := R1 ! lst[2];
    elif #lst eq 5 then 
        a1 := R1 ! lst[1]; a2 := R1 ! lst[2]; a3 := R1 ! lst[3]; 
        a4 := R1 ! lst[4]; a6 := R1 ! lst[5];
    else
        return "Error: wrong number of coeffs";
    end if;

	X3_001 := R2!(-a2*X1^2*X2*Z2 + a1*X1^2*Y2*Z2 - a4*X1^2*Z2^2 + 2*X1*Y1*Y2*Z2 + a3*X1*Y1*Z2^2 + a2*X1*Z1*X2^2 + X1*Z1*Y2^2 + 2*a3*X1*Z1*Y2*Z2 -
            3*a6*X1*Z1*Z2^2 - Y1^2*X2*Z2 - a1*Y1*Z1*X2^2 - 2*Y1*Z1*X2*Y2 - 2*a3*Y1*Z1*X2*Z2 + a4*Z1^2*X2^2 - a3*Z1^2*X2*Y2 + 3*a6*Z1^2*X2*Z2);
        
	Y3_001 := R2!(-3*X1^2*X2*Y2 + (a1*a2 - 3*a3)*X1^2*X2*Z2 + (-a1^2 - a2)*X1^2*Y2*Z2 + (a1*a4 - a2*a3)*X1^2*Z2^2 + 3*X1*Y1*X2^2 + 2*a2*X1*Y1*X2*Z2 -
            2*a1*X1*Y1*Y2*Z2 + a4*X1*Y1*Z2^2 + (-a1*a2 + 3*a3)*X1*Z1*X2^2 - 2*a2*X1*Z1*X2*Y2 + (-2*a1*a3 - 2*a4)*X1*Z1*Y2*Z2 + (3*a1*a6 -
            a3*a4)*X1*Z1*Z2^2 - Y1^2*Y2*Z2 + (a1^2 + a2)*Y1*Z1*X2^2 + 2*a1*Y1*Z1*X2*Y2 + (2*a1*a3 + 2*a4)*Y1*Z1*X2*Z2 + Y1*Z1*Y2^2 + (a3^2
            + 3*a6)*Y1*Z1*Z2^2 + (-a1*a4 + a2*a3)*Z1^2*X2^2 - a4*Z1^2*X2*Y2 + (-3*a1*a6 + a3*a4)*Z1^2*X2*Z2 + (-a3^2 - 3*a6)*Z1^2*Y2*Z2);
    
	Z3_001 := R2!(3*X1^2*X2*Z2 + a2*X1^2*Z2^2 - a1*X1*Y1*Z2^2 - 3*X1*Z1*X2^2 + a4*X1*Z1*Z2^2 - Y1^2*Z2^2 - a3*Y1*Z1*Z2^2 - a2*Z1^2*X2^2 +
            a1*Z1^2*X2*Y2 - a4*Z1^2*X2*Z2 + Z1^2*Y2^2 + a3*Z1^2*Y2*Z2);
        
        
    X3_010 := R2!(-a1*a2*X1^2*X2^2 - a2*X1^2*X2*Y2 + (-a1^2*a3 - 2*a1*a4)*X1^2*X2*Z2 + (-a1*a3 - a4)*X1^2*Y2*Z2 + (-a1*a3^2 - 3*a1*a6)*X1^2*Z2^2 + (a1^2 -
            a2)*X1*Y1*X2^2 + 2*a1*X1*Y1*X2*Y2 - 2*a4*X1*Y1*X2*Z2 + X1*Y1*Y2^2 + (-a3^2 - 3*a6)*X1*Y1*Z2^2 + (-a1*a4 - a2*a3)*X1*Z1*X2^2 -
            2*a4*X1*Z1*X2*Y2 + (-2*a1*a3^2 - 6*a1*a6 - 2*a3*a4)*X1*Z1*X2*Z2 + (-2*a3^2 - 6*a6)*X1*Z1*Y2*Z2 + (-a1^3*a6 + a1^2*a3*a4 - a1*a2*a3^2 -
            4*a1*a2*a6 + a1*a4^2 - a3^3 - 3*a3*a6)*X1*Z1*Z2^2 + a1*Y1^2*X2^2 + Y1^2*X2*Y2 + a3*Y1^2*X2*Z2 + (a1*a3 - a4)*Y1*Z1*X2^2 +
            2*a3*Y1*Z1*X2*Y2 - 6*a6*Y1*Z1*X2*Z2 + (-a1^2*a6 + a1*a3*a4 - a2*a3^2 - 4*a2*a6 + a4^2)*Y1*Z1*Z2^2 - a3*a4*Z1^2*X2^2 - 3*a6*Z1^2*X2*Y2 +
            (-a3^3 - 6*a3*a6)*Z1^2*X2*Z2 + (-a1^2*a6 + a1*a3*a4 - a2*a3^2 - 4*a2*a6 + a4^2)*Z1^2*Y2*Z2 + (-a1^2*a3*a6 + a1*a3^2*a4 - a2*a3^3 -
            4*a2*a3*a6 + a3*a4^2)*Z1^2*Z2^2);
			
	Y3_010 := R2!((-a2^2 + 3*a4)*X1^2*X2^2 + (a1^2*a4 - 2*a1*a2*a3 - a2*a4 + 3*a3^2 + 9*a6)*X1^2*X2*Z2 + (3*a1^2*a6 - 2*a1*a3*a4 + a2*a3^2 + 3*a2*a6 -
            a4^2)*X1^2*Z2^2 + (a1*a2 - 3*a3)*X1*Y1*X2^2 + (2*a1*a4 - 2*a2*a3)*X1*Y1*X2*Z2 + (3*a1*a6 - a3*a4)*X1*Y1*Z2^2 + (-a2*a4 +
            9*a6)*X1*Z1*X2^2 + (6*a1^2*a6 - 4*a1*a3*a4 + 2*a2*a3^2 + 12*a2*a6 - 4*a4^2)*X1*Z1*X2*Z2 + (a1^4*a6 - a1^3*a3*a4 + a1^2*a2*a3^2 +
            5*a1^2*a2*a6 - a1^2*a4^2 - a1*a2*a3*a4 - a1*a3^3 - 3*a1*a3*a6 + a2^2*a3^2 + 4*a2^2*a6 - a2*a4^2 - a3^2*a4 - 3*a4*a6)*X1*Z1*Z2^2 +
            a1*Y1^2*X2*Y2 + Y1^2*Y2^2 + a3*Y1^2*Y2*Z2 + (a1*a4 - a2*a3)*Y1*Z1*X2^2 + (6*a1*a6 - 2*a3*a4)*Y1*Z1*X2*Z2 + (a1^3*a6 -
            a1^2*a3*a4 + a1*a2*a3^2 + 4*a1*a2*a6 - a1*a4^2 - a3^3 - 3*a3*a6)*Y1*Z1*Z2^2 + (3*a2*a6 - a4^2)*Z1^2*X2^2 + (a1^2*a2*a6 -
            a1*a2*a3*a4 + 3*a1*a3*a6 + a2^2*a3^2 + 4*a2^2*a6 - a2*a4^2 - 2*a3^2*a4 - 3*a4*a6)*Z1^2*X2*Z2 + (a1^3*a3*a6 - a1^2*a3^2*a4 +
            a1^2*a4*a6 + a1*a2*a3^3 + 4*a1*a2*a3*a6 - 2*a1*a3*a4^2 + a2*a3^2*a4 + 4*a2*a4*a6 - a3^4 - 6*a3^2*a6 - a4^3 - 9*a6^2)*Z1^2*Z2^2);
    Z3_010 := R2!(3*a1*X1^2*X2^2 + 3*X1^2*X2*Y2 + (a1^3 + 2*a1*a2)*X1^2*X2*Z2 + (a1^2 + a2)*X1^2*Y2*Z2 + (a1^2*a3 + a1*a4)*X1^2*Z2^2 + 3*X1*Y1*X2^2 +
            (2*a1^2 + 2*a2)*X1*Y1*X2*Z2 + 2*a1*X1*Y1*Y2*Z2 + (2*a1*a3 + a4)*X1*Y1*Z2^2 + (a1*a2 + 3*a3)*X1*Z1*X2^2 + 2*a2*X1*Z1*X2*Y2 +
            (2*a1^2*a3 + 2*a1*a4 + 2*a2*a3)*X1*Z1*X2*Z2 + (2*a1*a3 + 2*a4)*X1*Z1*Y2*Z2 + (2*a1*a3^2 + 3*a1*a6 + a3*a4)*X1*Z1*Z2^2 +
            a1*Y1^2*X2*Z2 + Y1^2*Y2*Z2 + a3*Y1^2*Z2^2 + a2*Y1*Z1*X2^2 + (2*a1*a3 + 2*a4)*Y1*Z1*X2*Z2 + Y1*Z1*Y2^2 + 2*a3*Y1*Z1*Y2*Z2 + (2*a3^2 +
            3*a6)*Y1*Z1*Z2^2 + a2*a3*Z1^2*X2^2 + a4*Z1^2*X2*Y2 + (a1*a3^2 + 2*a3*a4)*Z1^2*X2*Z2 + (a3^2 + 3*a6)*Z1^2*Y2*Z2 + (a3^3 + 3*a3*a6)*Z1^2*Z2^2);

    return [[X3_001,Y3_001,Z3_001],[X3_010,Y3_010,Z3_010]];
end function;


//--> WC complete system of addition laws, outputs P+Q <--//

WC_AddL_eval := function(BF,lst,P,Q)

	// Things are a bit different for decimal numbers, so we need a flag for it
	if (Type(BF) eq Type(RealField())) or (Type(BF) eq Type(ComplexField())) then
		fl := true;
	else
		fl := false;
	end if;

	if fl then
		R1<a1,a2,a3,a4,a6> := PolynomialRing(BF,5);
	else
		R1<a1,a2,a3,a4,a6> := FunctionField(BF,5);
	end if;
    R2<X1,Y1,Z1, X2,Y2,Z2> := PolynomialRing(R1,6);

	// Get the coeffs
    if #lst eq 2 then
        a1 := R1 ! 0; a2 := R1 ! 0; a3 := R1 ! 0; 
        a4 := R1 ! lst[1]; a6 := R1 ! lst[2];
    elif #lst eq 5 then 
        a1 := R1 ! lst[1]; a2 := R1 ! lst[2]; a3 := R1 ! lst[3]; 
        a4 := R1 ! lst[4]; a6 := R1 ! lst[5];
    else
        return "Error: wrong number of coeffs";
    end if;

	// Get the points 
    Xp := R2!P[1]; Yp := R2!P[2]; Zp := R2!P[3]; 
    Xq := R2!Q[1]; Yq := R2!Q[2]; Zq := R2!Q[3];

	
	// Choose which to use
	if fl then 
		if (Abs(BF!(Xp - Xq)) gt 1E-10) and (Abs(BF!(Yp - Yq)) gt 1E-10) and (Abs(BF!(Zp - Zq)) gt 1E-10) then
			//here we use (a,b,c) = (0,1,0)
			
			X3 := R2!(-a1*a2*Xp^2*Xq^2 - a2*Xp^2*Xq*Yq + (-a1^2*a3 - 2*a1*a4)*Xp^2*Xq*Zq + (-a1*a3 - a4)*Xp^2*Yq*Zq + (-a1*a3^2 - 3*a1*a6)*Xp^2*Zq^2 + (a1^2 -
            a2)*Xp*Yp*Xq^2 + 2*a1*Xp*Yp*Xq*Yq - 2*a4*Xp*Yp*Xq*Zq + Xp*Yp*Yq^2 + (-a3^2 - 3*a6)*Xp*Yp*Zq^2 + (-a1*a4 - a2*a3)*Xp*Zp*Xq^2 -
            2*a4*Xp*Zp*Xq*Yq + (-2*a1*a3^2 - 6*a1*a6 - 2*a3*a4)*Xp*Zp*Xq*Zq + (-2*a3^2 - 6*a6)*Xp*Zp*Yq*Zq + (-a1^3*a6 + a1^2*a3*a4 - a1*a2*a3^2 -
            4*a1*a2*a6 + a1*a4^2 - a3^3 - 3*a3*a6)*Xp*Zp*Zq^2 + a1*Yp^2*Xq^2 + Yp^2*Xq*Yq + a3*Yp^2*Xq*Zq + (a1*a3 - a4)*Yp*Zp*Xq^2 +
            2*a3*Yp*Zp*Xq*Yq - 6*a6*Yp*Zp*Xq*Zq + (-a1^2*a6 + a1*a3*a4 - a2*a3^2 - 4*a2*a6 + a4^2)*Yp*Zp*Zq^2 - a3*a4*Zp^2*Xq^2 - 3*a6*Zp^2*Xq*Yq +
            (-a3^3 - 6*a3*a6)*Zp^2*Xq*Zq + (-a1^2*a6 + a1*a3*a4 - a2*a3^2 - 4*a2*a6 + a4^2)*Zp^2*Yq*Zq + (-a1^2*a3*a6 + a1*a3^2*a4 - a2*a3^3 -
            4*a2*a3*a6 + a3*a4^2)*Zp^2*Zq^2);
			
			Y3 := R2!((-a2^2 + 3*a4)*Xp^2*Xq^2 + (a1^2*a4 - 2*a1*a2*a3 - a2*a4 + 3*a3^2 + 9*a6)*Xp^2*Xq*Zq + (3*a1^2*a6 - 2*a1*a3*a4 + a2*a3^2 + 3*a2*a6 -
            a4^2)*Xp^2*Zq^2 + (a1*a2 - 3*a3)*Xp*Yp*Xq^2 + (2*a1*a4 - 2*a2*a3)*Xp*Yp*Xq*Zq + (3*a1*a6 - a3*a4)*Xp*Yp*Zq^2 + (-a2*a4 +
            9*a6)*Xp*Zp*Xq^2 + (6*a1^2*a6 - 4*a1*a3*a4 + 2*a2*a3^2 + 12*a2*a6 - 4*a4^2)*Xp*Zp*Xq*Zq + (a1^4*a6 - a1^3*a3*a4 + a1^2*a2*a3^2 +
            5*a1^2*a2*a6 - a1^2*a4^2 - a1*a2*a3*a4 - a1*a3^3 - 3*a1*a3*a6 + a2^2*a3^2 + 4*a2^2*a6 - a2*a4^2 - a3^2*a4 - 3*a4*a6)*Xp*Zp*Zq^2 +
            a1*Yp^2*Xq*Yq + Yp^2*Yq^2 + a3*Yp^2*Yq*Zq + (a1*a4 - a2*a3)*Yp*Zp*Xq^2 + (6*a1*a6 - 2*a3*a4)*Yp*Zp*Xq*Zq + (a1^3*a6 -
            a1^2*a3*a4 + a1*a2*a3^2 + 4*a1*a2*a6 - a1*a4^2 - a3^3 - 3*a3*a6)*Yp*Zp*Zq^2 + (3*a2*a6 - a4^2)*Zp^2*Xq^2 + (a1^2*a2*a6 -
            a1*a2*a3*a4 + 3*a1*a3*a6 + a2^2*a3^2 + 4*a2^2*a6 - a2*a4^2 - 2*a3^2*a4 - 3*a4*a6)*Zp^2*Xq*Zq + (a1^3*a3*a6 - a1^2*a3^2*a4 +
            a1^2*a4*a6 + a1*a2*a3^3 + 4*a1*a2*a3*a6 - 2*a1*a3*a4^2 + a2*a3^2*a4 + 4*a2*a4*a6 - a3^4 - 6*a3^2*a6 - a4^3 - 9*a6^2)*Zp^2*Zq^2);
			
			Z3 := R2!(3*a1*Xp^2*Xq^2 + 3*Xp^2*Xq*Yq + (a1^3 + 2*a1*a2)*Xp^2*Xq*Zq + (a1^2 + a2)*Xp^2*Yq*Zq + (a1^2*a3 + a1*a4)*Xp^2*Zq^2 + 3*Xp*Yp*Xq^2 +
            (2*a1^2 + 2*a2)*Xp*Yp*Xq*Zq + 2*a1*Xp*Yp*Yq*Zq + (2*a1*a3 + a4)*Xp*Yp*Zq^2 + (a1*a2 + 3*a3)*Xp*Zp*Xq^2 + 2*a2*Xp*Zp*Xq*Yq +
            (2*a1^2*a3 + 2*a1*a4 + 2*a2*a3)*Xp*Zp*Xq*Zq + (2*a1*a3 + 2*a4)*Xp*Zp*Yq*Zq + (2*a1*a3^2 + 3*a1*a6 + a3*a4)*Xp*Zp*Zq^2 +
            a1*Yp^2*Xq*Zq + Yp^2*Yq*Zq + a3*Yp^2*Zq^2 + a2*Yp*Zp*Xq^2 + (2*a1*a3 + 2*a4)*Yp*Zp*Xq*Zq + Yp*Zp*Yq^2 + 2*a3*Yp*Zp*Yq*Zq + (2*a3^2 +
            3*a6)*Yp*Zp*Zq^2 + a2*a3*Zp^2*Xq^2 + a4*Zp^2*Xq*Yq + (a1*a3^2 + 2*a3*a4)*Zp^2*Xq*Zq + (a3^2 + 3*a6)*Zp^2*Yq*Zq + (a3^3 + 3*a3*a6)*Zp^2*Zq^2);
		
		else
			//here we use (a,b,c) = (0,0,1)
			X3 := R2!(-a2*Xp^2*Xq*Zq + a1*Xp^2*Yq*Zq - a4*Xp^2*Zq^2 + 2*Xp*Yp*Yq*Zq + a3*Xp*Yp*Zq^2 + a2*Xp*Zp*Xq^2 + Xp*Zp*Yq^2 + 2*a3*Xp*Zp*Yq*Zq -
            3*a6*Xp*Zp*Zq^2 - Yp^2*Xq*Zq - a1*Yp*Zp*Xq^2 - 2*Yp*Zp*Xq*Yq - 2*a3*Yp*Zp*Xq*Zq + a4*Zp^2*Xq^2 - a3*Zp^2*Xq*Yq + 3*a6*Zp^2*Xq*Zq);
        
			Y3 := R2!(-3*Xp^2*Xq*Yq + (a1*a2 - 3*a3)*Xp^2*Xq*Zq + (-a1^2 - a2)*Xp^2*Yq*Zq + (a1*a4 - a2*a3)*Xp^2*Zq^2 + 3*Xp*Yp*Xq^2 + 2*a2*Xp*Yp*Xq*Zq -
            2*a1*Xp*Yp*Yq*Zq + a4*Xp*Yp*Zq^2 + (-a1*a2 + 3*a3)*Xp*Zp*Xq^2 - 2*a2*Xp*Zp*Xq*Yq + (-2*a1*a3 - 2*a4)*Xp*Zp*Yq*Zq + (3*a1*a6 -
            a3*a4)*Xp*Zp*Zq^2 - Yp^2*Yq*Zq + (a1^2 + a2)*Yp*Zp*Xq^2 + 2*a1*Yp*Zp*Xq*Yq + (2*a1*a3 + 2*a4)*Yp*Zp*Xq*Zq + Yp*Zp*Yq^2 + (a3^2
            + 3*a6)*Yp*Zp*Zq^2 + (-a1*a4 + a2*a3)*Zp^2*Xq^2 - a4*Zp^2*Xq*Yq + (-3*a1*a6 + a3*a4)*Zp^2*Xq*Zq + (-a3^2 - 3*a6)*Zp^2*Yq*Zq);
    
			Z3 := R2!(3*Xp^2*Xq*Zq + a2*Xp^2*Zq^2 - a1*Xp*Yp*Zq^2 - 3*Xp*Zp*Xq^2 + a4*Xp*Zp*Zq^2 - Yp^2*Zq^2 - a3*Yp*Zp*Zq^2 - a2*Zp^2*Xq^2 +
            a1*Zp^2*Xq*Yq - a4*Zp^2*Xq*Zq + Zp^2*Yq^2 + a3*Zp^2*Yq*Zq);

		end if;
	else 
		if (Xp eq Xq) and (Yp eq Yq) and (Zp eq Zq) then
			//here we use (a,b,c) = (0,1,0)
			
			X3 := R2!(-a1*a2*Xp^2*Xq^2 - a2*Xp^2*Xq*Yq + (-a1^2*a3 - 2*a1*a4)*Xp^2*Xq*Zq + (-a1*a3 - a4)*Xp^2*Yq*Zq + (-a1*a3^2 - 3*a1*a6)*Xp^2*Zq^2 + (a1^2 -
            a2)*Xp*Yp*Xq^2 + 2*a1*Xp*Yp*Xq*Yq - 2*a4*Xp*Yp*Xq*Zq + Xp*Yp*Yq^2 + (-a3^2 - 3*a6)*Xp*Yp*Zq^2 + (-a1*a4 - a2*a3)*Xp*Zp*Xq^2 -
            2*a4*Xp*Zp*Xq*Yq + (-2*a1*a3^2 - 6*a1*a6 - 2*a3*a4)*Xp*Zp*Xq*Zq + (-2*a3^2 - 6*a6)*Xp*Zp*Yq*Zq + (-a1^3*a6 + a1^2*a3*a4 - a1*a2*a3^2 -
            4*a1*a2*a6 + a1*a4^2 - a3^3 - 3*a3*a6)*Xp*Zp*Zq^2 + a1*Yp^2*Xq^2 + Yp^2*Xq*Yq + a3*Yp^2*Xq*Zq + (a1*a3 - a4)*Yp*Zp*Xq^2 +
            2*a3*Yp*Zp*Xq*Yq - 6*a6*Yp*Zp*Xq*Zq + (-a1^2*a6 + a1*a3*a4 - a2*a3^2 - 4*a2*a6 + a4^2)*Yp*Zp*Zq^2 - a3*a4*Zp^2*Xq^2 - 3*a6*Zp^2*Xq*Yq +
            (-a3^3 - 6*a3*a6)*Zp^2*Xq*Zq + (-a1^2*a6 + a1*a3*a4 - a2*a3^2 - 4*a2*a6 + a4^2)*Zp^2*Yq*Zq + (-a1^2*a3*a6 + a1*a3^2*a4 - a2*a3^3 -
            4*a2*a3*a6 + a3*a4^2)*Zp^2*Zq^2);
			
			Y3 := R2!((-a2^2 + 3*a4)*Xp^2*Xq^2 + (a1^2*a4 - 2*a1*a2*a3 - a2*a4 + 3*a3^2 + 9*a6)*Xp^2*Xq*Zq + (3*a1^2*a6 - 2*a1*a3*a4 + a2*a3^2 + 3*a2*a6 -
            a4^2)*Xp^2*Zq^2 + (a1*a2 - 3*a3)*Xp*Yp*Xq^2 + (2*a1*a4 - 2*a2*a3)*Xp*Yp*Xq*Zq + (3*a1*a6 - a3*a4)*Xp*Yp*Zq^2 + (-a2*a4 +
            9*a6)*Xp*Zp*Xq^2 + (6*a1^2*a6 - 4*a1*a3*a4 + 2*a2*a3^2 + 12*a2*a6 - 4*a4^2)*Xp*Zp*Xq*Zq + (a1^4*a6 - a1^3*a3*a4 + a1^2*a2*a3^2 +
            5*a1^2*a2*a6 - a1^2*a4^2 - a1*a2*a3*a4 - a1*a3^3 - 3*a1*a3*a6 + a2^2*a3^2 + 4*a2^2*a6 - a2*a4^2 - a3^2*a4 - 3*a4*a6)*Xp*Zp*Zq^2 +
            a1*Yp^2*Xq*Yq + Yp^2*Yq^2 + a3*Yp^2*Yq*Zq + (a1*a4 - a2*a3)*Yp*Zp*Xq^2 + (6*a1*a6 - 2*a3*a4)*Yp*Zp*Xq*Zq + (a1^3*a6 -
            a1^2*a3*a4 + a1*a2*a3^2 + 4*a1*a2*a6 - a1*a4^2 - a3^3 - 3*a3*a6)*Yp*Zp*Zq^2 + (3*a2*a6 - a4^2)*Zp^2*Xq^2 + (a1^2*a2*a6 -
            a1*a2*a3*a4 + 3*a1*a3*a6 + a2^2*a3^2 + 4*a2^2*a6 - a2*a4^2 - 2*a3^2*a4 - 3*a4*a6)*Zp^2*Xq*Zq + (a1^3*a3*a6 - a1^2*a3^2*a4 +
            a1^2*a4*a6 + a1*a2*a3^3 + 4*a1*a2*a3*a6 - 2*a1*a3*a4^2 + a2*a3^2*a4 + 4*a2*a4*a6 - a3^4 - 6*a3^2*a6 - a4^3 - 9*a6^2)*Zp^2*Zq^2);
			
			Z3 := R2!(3*a1*Xp^2*Xq^2 + 3*Xp^2*Xq*Yq + (a1^3 + 2*a1*a2)*Xp^2*Xq*Zq + (a1^2 + a2)*Xp^2*Yq*Zq + (a1^2*a3 + a1*a4)*Xp^2*Zq^2 + 3*Xp*Yp*Xq^2 +
            (2*a1^2 + 2*a2)*Xp*Yp*Xq*Zq + 2*a1*Xp*Yp*Yq*Zq + (2*a1*a3 + a4)*Xp*Yp*Zq^2 + (a1*a2 + 3*a3)*Xp*Zp*Xq^2 + 2*a2*Xp*Zp*Xq*Yq +
            (2*a1^2*a3 + 2*a1*a4 + 2*a2*a3)*Xp*Zp*Xq*Zq + (2*a1*a3 + 2*a4)*Xp*Zp*Yq*Zq + (2*a1*a3^2 + 3*a1*a6 + a3*a4)*Xp*Zp*Zq^2 +
            a1*Yp^2*Xq*Zq + Yp^2*Yq*Zq + a3*Yp^2*Zq^2 + a2*Yp*Zp*Xq^2 + (2*a1*a3 + 2*a4)*Yp*Zp*Xq*Zq + Yp*Zp*Yq^2 + 2*a3*Yp*Zp*Yq*Zq + (2*a3^2 +
            3*a6)*Yp*Zp*Zq^2 + a2*a3*Zp^2*Xq^2 + a4*Zp^2*Xq*Yq + (a1*a3^2 + 2*a3*a4)*Zp^2*Xq*Zq + (a3^2 + 3*a6)*Zp^2*Yq*Zq + (a3^3 + 3*a3*a6)*Zp^2*Zq^2);
		else
			//here we use (a,b,c) = (0,0,1)
			
			X3 := R2!(-a2*Xp^2*Xq*Zq + a1*Xp^2*Yq*Zq - a4*Xp^2*Zq^2 + 2*Xp*Yp*Yq*Zq + a3*Xp*Yp*Zq^2 + a2*Xp*Zp*Xq^2 + Xp*Zp*Yq^2 + 2*a3*Xp*Zp*Yq*Zq -
            3*a6*Xp*Zp*Zq^2 - Yp^2*Xq*Zq - a1*Yp*Zp*Xq^2 - 2*Yp*Zp*Xq*Yq - 2*a3*Yp*Zp*Xq*Zq + a4*Zp^2*Xq^2 - a3*Zp^2*Xq*Yq + 3*a6*Zp^2*Xq*Zq);
        
			Y3 := R2!(-3*Xp^2*Xq*Yq + (a1*a2 - 3*a3)*Xp^2*Xq*Zq + (-a1^2 - a2)*Xp^2*Yq*Zq + (a1*a4 - a2*a3)*Xp^2*Zq^2 + 3*Xp*Yp*Xq^2 + 2*a2*Xp*Yp*Xq*Zq -
            2*a1*Xp*Yp*Yq*Zq + a4*Xp*Yp*Zq^2 + (-a1*a2 + 3*a3)*Xp*Zp*Xq^2 - 2*a2*Xp*Zp*Xq*Yq + (-2*a1*a3 - 2*a4)*Xp*Zp*Yq*Zq + (3*a1*a6 -
            a3*a4)*Xp*Zp*Zq^2 - Yp^2*Yq*Zq + (a1^2 + a2)*Yp*Zp*Xq^2 + 2*a1*Yp*Zp*Xq*Yq + (2*a1*a3 + 2*a4)*Yp*Zp*Xq*Zq + Yp*Zp*Yq^2 + (a3^2
            + 3*a6)*Yp*Zp*Zq^2 + (-a1*a4 + a2*a3)*Zp^2*Xq^2 - a4*Zp^2*Xq*Yq + (-3*a1*a6 + a3*a4)*Zp^2*Xq*Zq + (-a3^2 - 3*a6)*Zp^2*Yq*Zq);
    
			Z3 := R2!(3*Xp^2*Xq*Zq + a2*Xp^2*Zq^2 - a1*Xp*Yp*Zq^2 - 3*Xp*Zp*Xq^2 + a4*Xp*Zp*Zq^2 - Yp^2*Zq^2 - a3*Yp*Zp*Zq^2 - a2*Zp^2*Xq^2 +
            a1*Zp^2*Xq*Yq - a4*Zp^2*Xq*Zq + Zp^2*Yq^2 + a3*Zp^2*Yq*Zq);
			
		end if;
	end if;
		
    return [X3,Y3,Z3];
end function;


//--> TwEdC complete system of addition laws, outputs the system of add. laws <--//

TwEd_AddL := function(BF,lst : except := false)

	if except then 
	    R1<a,d,sa,sd> := FunctionField(BF,4);
	else
	    R1<a,d> := FunctionField(BF,2);
	end if;
	
    R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8);
    if #lst eq 2 then
        a := R1 ! lst[1]; d := R1 ! lst[2];
    else
        return "Error: wrong number of coeffs";
    end if;

    X3_o := R2 ! X1*Y2*Z2*T1 + X2*Y1*Z1*T2;
    Z3_o := R2 ! Z1*Z2*T1*T2 + d*X1*X2*Y1*Y2;
    Y3_o := R2 ! Y1*Y2*Z1*Z2 - a*X1*X2*T1*T2;
    T3_o := R2 ! Z1*Z2*T1*T2 - d*X1*X2*Y1*Y2;

    X3_d := R2 ! X1*Y1*Z2*T2 + X2*Y2*Z1*T1;
    Z3_d := R2 ! a*X1*X2*T1*T2 + Y1*Y2*Z1*Z2;
    Y3_d := R2 ! X1*Y1*Z2*T2 - X2*Y2*Z1*T1;
    T3_d := R2 ! X1*Y2*Z2*T1 - X2*Y1*Z1*T2;

    return [[[X3_o, Z3_o], [Y3_o, T3_o]],[[X3_d, Z3_d], [Y3_d, T3_d]]];
end function;


//--> TwEdC complete system of addition laws, outputs P+Q <--//

TwEd_AddL_eval := function(BF,lst,P,Q)

	// Things are a bit different for decimal numbers, so we need a flag for it
	if (Type(BF) eq Type(RealField())) or (Type(BF) eq Type(ComplexField())) then
		fl := true;
	else
		fl := false;
	end if;

	if fl then
		R1<a,d> := PolynomialRing(BF,2);
	else
		R1<a,d> := FunctionField(BF,2);
	end if;
    R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8);

	// Get the coeffs
    if #lst eq 2 then
        a := R2!lst[1]; d := R2!lst[2];
    else
        return "Error: wrong number of coeffs";
    end if;

    // Get the points
    Xp := R2 ! P[1][1]; Zp := R2 ! P[1][2]; Yp := R2 ! P[2][1]; Tp := R2 ! P[2][2]; 
    Xq := R2 ! Q[1][1]; Zq := R2 ! Q[1][2]; Yq := R2 ! Q[2][1]; Tq := R2 ! Q[2][2]; 

	// Compute the sums
    X_o := Evaluate(R2 ! X1*Y2*Z2*T1 + X2*Y1*Z1*T2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    Z_o := Evaluate(R2 ! Z1*Z2*T1*T2 + d*X1*X2*Y1*Y2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    Y_o := Evaluate(R2 ! Y1*Y2*Z1*Z2 - a*X1*X2*T1*T2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    T_o := Evaluate(R2 ! Z1*Z2*T1*T2 - d*X1*X2*Y1*Y2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);

    X_d := Evaluate(R2 ! X1*Y1*Z2*T2 + X2*Y2*Z1*T1,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    Z_d := Evaluate(R2 ! a*X1*X2*T1*T2 + Y1*Y2*Z1*Z2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    Y_d := Evaluate(R2 ! X1*Y1*Z2*T2 - X2*Y2*Z1*T1,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    T_d := Evaluate(R2 ! X1*Y2*Z2*T1 - X2*Y1*Z1*T2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);

	
	// Choose which to use
	if fl then 
		if (Abs(BF!X_o) gt 1E-10 or Abs(BF!Z_o) gt 1E-10) and (Abs(BF!Y_o) gt 1E-10 or Abs(BF!T_o) gt 1E-10) then
			X := X_o; Y := Y_o; Z := Z_o; T := T_o;
		elif (Abs(BF!X_d) gt 1E-10 or Abs(BF!Z_d) gt 1E-10) and (Abs(BF!Y_d) gt 1E-10 or Abs(BF!T_d) gt 1E-10) then
			X := X_d; Y := Y_d; Z := Z_d; T := T_d;
		else
			return "Error: all are zero";
		end if;
	else 
		if (X_o ne 0 or Z_o ne 0) and (Y_o ne 0 or T_o ne 0) then
			X := X_o; Y := Y_o; Z := Z_o; T := T_o;
		elif (X_d ne 0 or Z_d ne 0) and (Y_d ne 0 or T_d ne 0) then
			X := X_d; Y := Y_d; Z := Z_d; T := T_d;
		else
			return "Error: all are zero";
		end if;
	end if;
		
    return [[X,Z],[Y,T]];
end function;


//--> TwEdC complete system of addition laws, outputs the result of both addition laws <--//

TwEd_AddL_evals := function(BF,lst,P,Q)

	// Things are a bit different for decimal numbers, so we need a flag for it
	if (Type(BF) eq Type(RealField())) or (Type(BF) eq Type(ComplexField())) then
		fl := true;
	else
		fl := false;
	end if;

	if fl then
		R1<a,d> := PolynomialRing(BF,2);
	else
		R1<a,d> := FunctionField(BF,2);
	end if;
    R2<X1,Y1,Z1,T1, X2,Y2,Z2,T2> := PolynomialRing(R1,8);

	// Get the coeffs
    if #lst eq 2 then
        a := R2!lst[1]; d := R2!lst[2];
    else
        return "Error: wrong number of coeffs";
    end if;

    // Get the points with z-coordinate 1
    Xp := R2 ! P[1][1]; Zp := R2 ! P[1][2]; Yp := R2 ! P[2][1]; Tp := R2 ! P[2][2]; 
    Xq := R2 ! Q[1][1]; Zq := R2 ! Q[1][2]; Yq := R2 ! Q[2][1]; Tq := R2 ! Q[2][2]; 

	// Compute the sums
    X_o := Evaluate(R2 ! X1*Y2*Z2*T1 + X2*Y1*Z1*T2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    Z_o := Evaluate(R2 ! Z1*Z2*T1*T2 + d*X1*X2*Y1*Y2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    Y_o := Evaluate(R2 ! Y1*Y2*Z1*Z2 - a*X1*X2*T1*T2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    T_o := Evaluate(R2 ! Z1*Z2*T1*T2 - d*X1*X2*Y1*Y2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);

    X_d := Evaluate(R2 ! X1*Y1*Z2*T2 + X2*Y2*Z1*T1,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    Z_d := Evaluate(R2 ! a*X1*X2*T1*T2 + Y1*Y2*Z1*Z2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    Y_d := Evaluate(R2 ! X1*Y1*Z2*T2 - X2*Y2*Z1*T1,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);
    T_d := Evaluate(R2 ! X1*Y2*Z2*T1 - X2*Y1*Z1*T2,[Xp,Yp,Zp,Tp, Xq,Yq,Zq,Tq]);

	
	// Choose which to use
	if fl then 
		if (Abs(BF!X_o) gt 1E-10 or Abs(BF!Z_o) gt 1E-10) and (Abs(BF!Y_o) gt 1E-10 or Abs(BF!T_o) gt 1E-10) then
			X := X_o; Y := Y_o; Z := Z_o; T := T_o;
		elif (Abs(BF!X_d) gt 1E-10 or Abs(BF!Z_d) gt 1E-10) and (Abs(BF!Y_d) gt 1E-10 or Abs(BF!T_d) gt 1E-10) then
			X := X_d; Y := Y_d; Z := Z_d; T := T_d;
		else
			return "Error: all are zero";
		end if;
	else 
		if (X_o ne 0 or Z_o ne 0) and (Y_o ne 0 or T_o ne 0) then
			X := X_o; Y := Y_o; Z := Z_o; T := T_o;
		elif (X_d ne 0 or Z_d ne 0) and (Y_d ne 0 or T_d ne 0) then
			X := X_d; Y := Y_d; Z := Z_d; T := T_d;
		else
			return "Error: all are zero";
		end if;
	end if;
		
    return [[[X_o,Z_o],[Y_o,T_o]],[[X_d,Z_d],[Y_d,T_d]]];
end function;


//--> WC to EdC, using BL birat. eq. , outputs the map symbolically <--//

W2E_P4_symb := function(BF,lst,P4 : PP4 := [0,0])
	// lst := [a1,a2,a3,a4,a6] can be numerical coeffs or letters

	R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
	Rw_var<u,v> := FunctionField(R,2);

	a1 := R ! lst[1]; a2 := (R ! lst[2]); a3 := (R ! lst[3]); a4 := (R ! lst[4]); a6 := (R ! lst[5]);
	up := (R ! P4[1]); vp := (R ! P4[2]);

	u2p := R ! PP4[1]; v2p := R ! PP4[2];
	
	up +:= u2p; vp +:= (a1*up + a3)/2;
	u +:= u2p; v +:= (a1*u + a3)/2;
	
	a4 +:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - ((a1^2)*u2p)/2;
	a2 +:= -3*u2p + (a1^2)/4;
	
	d := 1-4*(up^3)/(vp^2);
	c_x := (vp*u)/(up*v); c_y := (u-up)/(u+up);

	return [*[d],[c_x,c_y]*];

end function;


//--> EdC to WC, using BL birat. eq. , outputs the map symbolically <--//

E2W_P4_symb := function(BF,lst,P4 : PP4 := [0,0], aa1 := 0, aa3 := 0)
	//lst := [d] can be numerical coeff or letter
	
	R<a1,a2,a3,a4,a6,up,vp,d> := FunctionField(BF,8);
	Rd_var<x,y> := FunctionField(R,2);
	
	d := R ! lst[1];
	up := R ! P4[1]; vp := R ! P4[2];
	u2p := R ! PP4[1]; v2p := R ! PP4[2];

	a1 := R ! aa1; a2 := 2*(1+d)/(1-d)*up; a3 := R ! aa3; a4 := up^2; a6 := 0;
	c_u := up*(1+y)/(1-y); c_v := vp*(1+y)/(x*(1-y));

	c_u -:= u2p; c_v -:= (a1*(c_u-u2p)+a3)/2;
	a2 -:= -3*u2p + (a1^2)/4;
	a4 -:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - ((a1^2)*u2p)/2;	
	a6 := - u2p^3 - a2*u2p^2 + a4*u2p - (a1*u2p - a3)^2/4;
	
	return [*[a1,a2,a3,a4,a6],[c_u,c_v]*];

end function;



//--> WC to EdC, using MD birat. eq. , outputs the map symbolically <--//

W2E_symb := function(BF,lst)
	// lst := [a1,a2,a3,a4,a6] can be numerical coeffs or letters

	R<a1,a2,a3,a4,a6,d> := FunctionField(BF,6);
	Rw_var<u,v> := FunctionField(R,2);

	a1 := R ! lst[1]; a2 := (R ! lst[2]);a3 := (R ! lst[3]);a4 := (R ! lst[4]);a6 := (R ! lst[5]);

	d := (a2/2)-1;
	c_x := -2*u/v; c_y := (v^2-(2+2*d)*u^2 - 2*u^3)/(4*d*u^2-v^2);

	return [*[d],[c_x,c_y]*];

end function;


//--> EdC to WC, using MD birat. eq. , outputs the map symbolically <--//

E2W_symb := function(BF,lst)
	//lst := [d] can be numerical coeff or letter
	
	R<a1,a2,a3,a4,a6,d> := FunctionField(BF,6);
	Rd_var<x,y> := FunctionField(R,2);
	
	d := R ! lst[1];

	a1 := 0; a2 := 2*(d+1); a3 := 0; a4 := (d-1)^2; a6 := 0;
	
	A := 2*y - (2*d*y + d +1)*x^2 + 2;
	c_u := A/(x^2); c_v := -2*A/(x^3);

	return [*[a1,a2,a3,a4,a6],[c_u,c_v]*];

end function;



//--> WC to EdC, using BL birat. eq. <--//

W2E_P4_eval := function(BF,lst,P4,p : PP4 := [0,0])

    if #lst eq 5 then
        a1 := BF!lst[1]; a2 := BF!lst[2]; a3 := BF!lst[3]; a4 := BF!lst[4]; a6 := BF!lst[5];
    elif #lst eq 2 then 
        a1 := BF!0; a2 := BF!0; a3 := BF!0; a4 := BF!lst[1]; a6 := BF!lst[2];
    else
        return "Error: wrong number of coeffs";
    end if;
	u2p := BF ! PP4[1];
	up := (BF ! P4[1]) + u2p; vp := (BF ! P4[2]) + (a1*up + a3)/2;
	u := BF!p[1] + u2p; v := BF!p[2] + (a1*u + a3)/2;
	a4 +:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - (a1^2)*u2p/2; a2 +:= - 3*u2p + (a1^2)/4; 
	
	denom_x := up*v; denom_y := (u+up); denom_d := vp^2;
	if denom_x eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_x eq 0).";
	elif denom_y eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_y eq 0).";
	elif denom_y eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_d eq 0).";
	else
		d := 1-4*(up^3)/(vp^2);
		c_x := vp*u/denom_x; c_y := (u-up)/denom_y;

		return [*[BF|d],[BF|c_x,c_y]*];
	end if;

end function;


//--> EdC to WC, using BL birat. eq. <--//

E2W_P4_eval := function(BF,lst,P4,p : PP4 := [0,0], a1 := 0, a3 := 0)

	up := (BF ! P4[1]) + BF!PP4[1]; vp := (BF ! P4[2]) + (a1*up + a3)/2;
	u2p := BF ! PP4[1]; v2p := BF ! PP4[2];
	d := BF!lst[1]; x := BF!p[1]; y := BF!p[2];

	denom_x := (1-y); denom_y := x*(1-y);
	if denom_x eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_x eq 0).";
	elif denom_y eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_y eq 0).";
	else
		a2 := 2*(1+d)/(1-d)*up; a4 := up^2; a6 := 0;
		c_u := up*(1+y)/denom_x; c_v:= vp*(1+y)/denom_y;

		c_u -:= u2p; c_v -:= (a1*c_u+a3)/2; 

		a2 -:= -3*u2p + (a1^2)/4;
		a4 -:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - ((a1^2)*u2p)/2;	
		a6 := - u2p^3 - a2*u2p^2 + a4*u2p - (a1*u2p - a3)^2/4;
			
		return [*[a1,a2,a3,a4,a6],[c_u,c_v]*];
	end if;

end function;


//--> WC to EdC, using MD birat. eq. <--//

W2E_eval := function(BF,lst,p: PP4 := [0,0])
	
    if #lst eq 5 then
        a1 := BF!lst[1]; a2 := BF!lst[2]; a3 := BF!lst[3]; a4 := BF!lst[4]; a6 := BF!lst[5];
    elif #lst eq 2 then 
        a1 := BF!0; a2 := BF!0; a3 := BF!0; a4 := BF!lst[1]; a6 := BF!lst[2];
    else
        return "Error: wrong number of coeffs";
    end if;
	u := BF!p[1]; v := BF!p[2];
	
	u2p := BF ! PP4[1];
	u := BF!p[1] + u2p; v := BF!p[2] + (a1*u + a3)/2;
	a4 +:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - (a1^2)*u2p/2; a2 +:= - 3*u2p + (a1^2)/4; 
	
	if a2 ne 0 then
		d := (a2/2)-1;
	elif (a2 eq 0) and (IsSquare(BF!(a4))) then
		d := Sqrt(a4) + 1;
	else
	   return "Error: d is a not a square over the basis field"; 
	end if;

	denom_x := v; denom_y := 2*(a2-2)*u^2-v^2;
	if denom_x eq 0 then
		return "Error: map not defined for this point since (denom_x eq 0).";
	elif denom_y eq 0 then
		return "Error: map not defined for this point since (denom_y eq 0).";
	else
		c_x := -2*u/denom_x; c_y := (v^2-a2*u^2 - 2*u^3)/denom_y;

		return [*[d],[c_x,c_y]*];
	end if;

end function;


//--> EdC to WC, using MD birat. eq. <--//

E2W_eval := function(BF,lst,p: PP4 := [0,0], a1 := 0, a3 := 0)

	d := BF!lst[1]; x := BF!p[1]; y := BF!p[2]; u2p := BF!PP4[1];

	denom_x := x^2; denom_y := x^3;
	if denom_x eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_x eq 0).";
	elif denom_y eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_y eq 0).";
	else
		a2 := 2*(d+1); a4 := (d-1)^2; a6 := 0;
		A := 2*y - (2*d*y + d +1)*x^2 + 2;
		c_u := A/(x^2); c_v:= -2*A/(x^3);

		c_u -:= u2p; c_v -:= (a1*c_u+a3)/2; 

		a2 -:= -3*u2p + (a1^2)/4;
		a4 -:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - ((a1^2)*u2p)/2;	
		a6 := - u2p^3 - a2*u2p^2 + a4*u2p - (a1*u2p - a3)^2/4;
		
		
		return [*[a1,a2,a3,a4,a6],[c_u,c_v]*];
	end if;

end function;


//--> WC to TwEC, using BL birat. eq. <--//

W2TwE_P4_eval := function(BF,lst,P4,p : PP4 := [0,0], a := 1)

	// If 'a' is not a square on BF, then we need to work on an extension field
	tf, sqrt_a := IsSquare(BF!a);
	if not tf then 
		// QF_int := Integers()!(Numerator(a)*Denominator(a));
		// BF<sqrt_a> := QuadraticField(QF_int);
		BF := AlgebraicClosure(BF);
		tf, sqrt_a := IsSquare(BF!a);
		// set a flag that we changed base field?;
	end if;
	
	a := BF!a;
	if #lst eq 5 then
        a1 := BF!lst[1]; a2 := BF!lst[2]; a3 := BF!lst[3]; a4 := BF!lst[4]; a6 := BF!lst[5];
    elif #lst eq 2 then 
        a1 := BF!0; a2 := BF!0; a3 := BF!0; a4 := BF!lst[1]; a6 := BF!lst[2];
    else
        return "Error: wrong number of coeffs";
    end if;
	u2p := BF ! PP4[1];
	up := (BF ! P4[1]) + u2p; vp := (BF ! P4[2]) + (a1*up + a3)/2;
	u := BF!p[1] + u2p; v := BF!p[2] + (a1*u + a3)/2;
	a4 +:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - (a1^2)*u2p/2; a2 +:= - 3*u2p + (a1^2)/4; 
	
	denom_x := up*v; denom_y := (u+up); denom_d := vp^2;
	if denom_x eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_x eq 0).";
	elif denom_y eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_y eq 0).";
	elif denom_y eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_d eq 0).";
	else
		d := (1-4*(up^3)/(vp^2))*a; //<---- is this right?
		c_x := vp*u/(denom_x*sqrt_a); c_y := (u-up)/denom_y; //<---- is this right?

		return [*BF,[BF|a,d],[BF|c_x,c_y]*];
	end if;

end function;


//--> TwEC to WC, using BL birat. eq. <--//

TwE2W_P4_eval := function(BF,lst,P4,p : PP4 := [0,0], a1 := 0, a3 := 0)

	a := BF!lst[1];
	
	// If 'a' is not a square on BF, then we need to work on an extension field
	tf, sqrt_a := IsSquare(BF!a);
	if not tf then 
		// QF_int := Integers()!(Numerator(a)*Denominator(a));
		// BF<sqrt_a> := QuadraticField(QF_int);
		BF := AlgebraicClosure(BF);
		tf, sqrt_a := IsSquare(BF!a);
		// set a flag that we changed base field?;
	end if;

	up := (BF ! P4[1]) + BF!PP4[1]; vp := (BF ! P4[2]) + (a1*up + a3)/2;
	u2p := BF ! PP4[1]; v2p := BF ! PP4[2];
	a := BF!a; d := (BF!lst[2])/a; x := (BF!p[1])*sqrt_a; y := BF!p[2]; //<---- is this right?
	
	denom_x := (1-y); denom_y := x*(1-y);
	if denom_x eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_x eq 0).";
	elif denom_y eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_y eq 0).";
	else
		a2 := 2*(1+d)/(1-d)*up; a4 := up^2;
		c_u := up*(1+y)/denom_x; c_v:= vp*(1+y)/denom_y;

		c_u -:= u2p; c_v -:= (a1*c_u+a3)/2; 

		a2 -:= -3*u2p + (a1^2)/4;
		a4 -:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - ((a1^2)*u2p)/2;	
		a6 := - u2p^3 - a2*u2p^2 + a4*u2p - (a1*u2p - a3)^2/4;
			
		return [*BF,[a1,a2,a3,a4,a6],[c_u,c_v]*];
	end if;

end function;



//--> WC to TwEdC, using the MD birat. eq. <--//

W2TwE_eval := function(BF,lst,p: PP4 := [0,0], a := 1)
	
	// If 'a' is not a square on BF, then we need to work on the algebraic closure of BF
    tf, sqrt_a := IsSquare(BF!a);
	if not tf then 
		//QF_int := Integers()!(Numerator(a)*Denominator(a));
		//BF<sqrt_a> := QuadraticField(QF_int); 
		BF := AlgebraicClosure(BF);
		tf, sqrt_a := IsSquare(BF!a);
		// set a flag that we changed base field?;
	end if;	

	a := BF!a;
	if #lst eq 5 then
        a1 := BF!lst[1]; a2 := BF!lst[2]; a3 := BF!lst[3]; a4 := BF!lst[4]; a6 := BF!lst[5];
    elif #lst eq 2 then 
        a1 := BF!0; a2 := BF!0; a3 := BF!0; a4 := BF!lst[1]; a6 := BF!lst[2];
    else
        return "Error: wrong number of coeffs";
    end if;
	u := BF!p[1]; v := BF!p[2];
	
	u2p := BF ! PP4[1];
	u := BF!p[1] + u2p; v := BF!p[2] + (a1*u + a3)/2;
	a4 +:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - (a1^2)*u2p/2; a2 +:= - 3*u2p + (a1^2)/4; 
	
	if a2 ne 0 then
		d := (a2/2)-1; 
	elif (a2 eq 0) and (IsSquare(BF!(a4))) then
		d := Sqrt(a4) + 1;
	else
	   return "Error: d is a not a square over the basis field"; 
	end if;
	d *:= a; //<---- is this right? 


	denom_x := v; denom_y := 2*(a2-2)*u^2-v^2;
	if denom_x eq 0 then
		return "Error: map not defined for this point since (denom_x eq 0).";
	elif denom_y eq 0 then
		return "Error: map not defined for this point since (denom_y eq 0).";
	else
		c_x := -2*u/(denom_x*sqrt_a); c_y := (v^2-a2*u^2 - 2*u^3)/denom_y; //<---- is this right?

		return [*BF,[a,d],[c_x,c_y]*];
	end if;

end function;

//--> TwEdC to WC, using the MD birat. eq. <--//

TwE2W_eval := function(BF,lst,p: PP4 := [0,0], a1 := 0, a3 := 0)

        a := BF!lst[1];
	// If 'a' is not a square on BF, then we need to work on the algebraic closure of BF
        tf, sqrt_a := IsSquare(BF!a);
	if not tf then 
		BF := AlgebraicClosure(BF);
		tf, sqrt_a := IsSquare(BF!a);
	end if;

	a := BF!a; d := (BF!lst[2])/a; x := (BF!p[1])*sqrt_a; y := BF!p[2]; u2p := BF!PP4[1];

	denom_x := x^2; denom_y := x^3;
	if denom_x eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_x eq 0).";
	elif denom_y eq 0 then
		return "Error: birat. eq. not defined for this point since (denom_y eq 0).";
	else
		a2 := 2*(d+1); a4 := (d-1)^2;
		A := 2*y - (2*d*y + d +1)*x^2 + 2;
		c_u := A/(x^2); c_v:= -2*A/(x^3);

		c_u -:= u2p; c_v -:= (a1*c_u+a3)/2; 

		a2 -:= -3*u2p + (a1^2)/4;
		a4 -:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - ((a1^2)*u2p)/2;	
		a6 := - u2p^3 - a2*u2p^2 + a4*u2p - (a1*u2p - a3)^2/4;
		
		
		return [*BF,[a1,a2,a3,a4,a6],[c_u,c_v]*];
	end if;

end function;


//--> WC to TwEdC, using 2-isogeny <--//

W2TwE_iso_eval := function(BF,lst,p : PP4 := [0,0], u_lst := [0,0,0])

    if #lst eq 5 then
        a1 := BF!lst[1]; a2 := BF!lst[2]; a3 := BF!lst[3]; a4 := BF!lst[4]; a6 := BF!lst[5];
    elif #lst eq 2 then 
        a1 := BF!0; a2 := BF!0; a3 := BF!0; a4 := BF!lst[1]; a6 := BF!lst[2];
    else
        return "Error: wrong number of coeffs";
    end if;
	u2p := BF ! PP4[1];
	a4 +:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - (a1^2)*u2p/2; a2 +:= - 3*u2p + (a1^2)/4; 	
	
	if u_lst eq [0,0,0] then
		PR<s> := PolynomialRing(BF);
		if a1 ne 0 or a3 ne 0 then
			a6 := 0;
		end if;
		rts := Roots(s^3 + a2*s^2 + a4*s + a6);
		if #rts lt 3 then
			BF := AlgebraicClosure(BF);
			PR<s> := PolynomialRing(BF);
			rts := Roots(s^3 + BF!a2*s^2 + BF!a4*s + BF!a6);
		end if;
		rts_abs := [Abs(rts[i][1]): i in [1..3]];
		ParallelSort(~rts_abs,~rts);
		u0 := rts[1][1]; u1 := rts[2][1] - u0; u2 := rts[3][1] - u0;
	else
		u0 := BF!u_lst[1]; u1 := BF!u_lst[2]; u2 := BF!u_lst[3];
	end if;
	u := BF!p[1] + u2p - u0; v := BF!p[2] + (a1*u + a3)/2;
			
	denom_x := u1*u2 - u^2; denom_y := v^2 + (u1-u2)*u^2;
	if denom_x eq 0 then
		return "Error: isogeny not defined for this point since (denom_x eq 0).";
	elif denom_y eq 0 then
		return "Error: isogeny not defined for this point since (denom_y eq 0).";
	else
		a := 4*u1; d := 4*u2;
		c_x := v/denom_x; c_y := (v^2 - (u1-u2)*u^2)/denom_y;

		return [* BF,[a,d],[c_x,c_y], u0 *];
	end if;

end function;


//--> TwEdC to WC, using 2-isogeny <--//

TwE2W_iso_eval := function(BF,lst,p : PP4 := [0,0], a1 := 0, a3 := 0, u0 := 0)

	u2p := PP4[1];
	if (u2p ne 0 or a1 ne 0 or a3 ne 0) and u0 ne 0 then
		return "Error: either specify u2p,a1,a3 != 0 or u0 != 0";
	end if;
	
	a := BF!lst[1]; d := BF!lst[2]; x := BF!p[1]; y := BF!p[2];

	denom_u := 4*x^2; denom_v := 32*x*(1+y)*(1-y);
	if denom_u eq 0 then
		return "Error: isogeny not defined for this point since (denom_u eq 0).";
	elif denom_v eq 0 then
		return "Error: isogeny not defined for this point since (denom_v eq 0).";
	else
		u1 := a/4 + u0; u2 := d/4 + u0;
		a2 := -(u0+u1+u2); a4 := u0*u1 + u0*u2 + u1*u2; a6 := -u0*u1*u2;
		c_u := 1/denom_u +u0; c_v := (a-d)*((1-y)^2-(1+y)^2)/denom_v;

		c_u -:= u2p; c_v -:= (a1*c_u+a3)/2; 		
		a2 -:= -3*u2p + (a1^2)/4;
		a4 -:= 3*u2p^2 - 2*a2*u2p + (a1*a3)/2 - ((a1^2)*u2p)/2;	
		a6 := - u2p^3 - a2*u2p^2 + a4*u2p - (a1*u2p - a3)^2/4;

		return [*[a1,a2,a3,a4,a6],[u0,u1,u2],[c_u,c_v]*];
	end if;

end function;
