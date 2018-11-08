
with(LinearAlgebra):
with(linalg):
with(SolveTools):
with(plots):
with(ArrayTools):

Digits := 20:
interface(rtablesize = 100):




AGenerator := proc(r,t)

	local L,Line1,Line2,Line3;

	Line1 := convert(Concatenate(2, Vector[row](<1,0>) ,Vector[row](2*r-2,0) ),list);
	Line2 := convert(Concatenate(2, Vector[row](<0,1>) ,Vector[row](2*r-2,0) ),list);
	
	L := [Line1, Line2];
	
	for i from 1 to r-1 do

		Line1 := convert(Concatenate(2, Vector[row](2*(i-1),0), Vector[row](<-t,0,1>), Vector[row](2*(r-i)-1,0) ),list);
		Line2 := convert(Concatenate(2, Vector[row](2*(i-1),0), Vector[row](<0,-t,1>), Vector[row](2*(r-i)-2,0) ),list);
		Line3 := convert(Concatenate(2, Vector[row](2*(i-1),0), Vector[row](<-(t^(1-(1/(2^i)))),-(t^(1-(1/(2^i)))),0,1>), Vector[row](2*(r-i-2),0) ),list);
			
		L := [op(L),Line1]:
		L := [op(L),Line2]:
		L := [op(L),Line3]:

	end do;
	

	Line1 := convert(Concatenate(2, Vector[row](2*r-2,0), Vector[row](<-1,0>)),list);
	Line2 := convert(Concatenate(2, Vector[row](2*r-2,0), Vector[row](<0,-1>)),list);
	
	L := [op(L),Line1]:
	L := [op(L),Line2]:

	
	return Matrix(3*r+1,2*r,L);
	
	
end proc:


BGenerator := proc(r,t)

	return Concatenate(1,Vector[column](<t^2,t>),Vector[column](3*r-1,0));
	
end proc:


Tropical_x := proc(lambda,r)
	

	local j,x, x1_old,x2_old,x1_new,x2_new ;

	x1_old := min(2,lambda);
	x2_old := 1;

	x := [x1_old,x2_old];
	
	for j from 1 to r-1 do:

		x1_new := 1 + min(x1_old,x2_old);
		x2_new := (1-1/(2^j)) + max(x1_old,x2_old);

		x := [op(x),x1_new,x2_new];

		x1_old := x1_new;
		x2_old := x2_new;
	end do:


	return convert(x, Vector[column]);


end proc:


Tropical_y:= proc(lambda,r)
	
	local j,x,y ;

	x := Tropical_x(lambda,r);

	y := [lambda - 2, lambda - 1];
	
	for j from 1 to r-1 do:
		
		y := [ op(y),lambda - 1 - x[2*j-1] ];
		y := [ op(y),lambda - 1 - x[ 2*j ] ];
		y := [ op(y),lambda - x[2*j+2] ];

	end do:
	
	y:= [ op(y),lambda - x[2*r-1] ];
	y:= [ op(y),lambda - x[2*r] ];

	return convert(y, Vector[column]);

end proc:











EquatlityConstraintsAGenerator := proc(r,t)
	local A;
	
	A := AGenerator(r,t);

	return Matrix(blockmatrix(1,2,[A,IdentityMatrix(3*r+1)]));

end proc;


EquatlityConstraintsBGenerator := proc(r,t)
	
local B;
	B := BGenerator(r,t);
	return B;
	
end proc;

#################################################################################################
#################################################################################################
#################################################################################################

MuBar := proc(z,n,m)
	local x,s;
	x  :=  z(1..n);
	s  :=  z(n+m+1..2*n+m);
	return DotProduct(Transpose(x),s) / n
end proc;


IsAdmissible := proc(z, A,b,c, n,m)
	local x,y,s, b1,b2,b3;
	x  :=  z(1..n);
	y  := z(n+1 .. n+m);
	s  :=  z(n+m+1..2*n+m);

	b1:= Equal(Multiply(A,x) , b);
	b2:= Equal( Multiply(Transpose(A) , y) + s , c);
	b3:= true;
	for i from 1 to n do 
		if evalf(x(i))<=0 then b3:= false end if;
		if evalf(s(i))<=0 then b3:= false end if;
		end do;
	return (b1 and b2 and b3);	
end proc;

IsInV := proc(theta,z, A,b,c, n,m)
	
	local x,y,s, b1,muBar;
	x  :=  z(1..n);
	y  := z(n+1 .. n+m);
	s  :=  z(n+m+1..2*n+m);

	b1:= IsAdmissible(z,A,b,c,n,m);
	muBar := MuBar(z,n,m);
	return (b1 and evalb( Norm( Multiply( Matrix(DiagonalMatrix(x)), s) - muBar * Vector(n,1) ,2) <= theta * muBar));
end proc;




NewtonMatrix := proc(n,m,A,z)
	local x,s, X,S;
	x  :=  z(1..n):
	s  :=  z(n+m+1..2*n+m):

	X := Matrix(DiagonalMatrix(x)):
	S := Matrix(DiagonalMatrix(s)):

	return Matrix(blockmatrix(3,3,[ZeroMatrix(n,n),  Transpose(A), IdentityMatrix(n), A, ZeroMatrix(m,m), ZeroMatrix(m,n), S, ZeroMatrix(n,m), X])):

end proc: 



GetNewtonDirection := proc(n,m,A,z,mu)
	
	local x,s, M,X,e,Vect;
	x  :=  z(1..n):
	s  :=  z(n+m+1..2*n+m):

	M  := NewtonMatrix(n,m,A,z):

	X      :=  Matrix(DiagonalMatrix(x)):
	e      :=  Vector(n,1):
	Vect   :=  convert(<Vector(n,0),Vector(m,0),mu.e - X.s>,Vector[column]):

	return LinearSolve(M,Vect):

end proc:



GetMaxAlpha := proc(n,m,z,dz,theta_prime)

	local dx,dX,ds,dS,x,X,s,S,e, v1,v2,v3,a0,a1,a2,a3,a4, Sol;
	dx :=  dz(1..n):
	dX :=  Matrix(DiagonalMatrix(dx)):

	ds :=  dz(n+m+1..2*n+m):
	dS :=  Matrix(DiagonalMatrix(ds)):

	x  :=  z(1..n):
	X  :=  Matrix(DiagonalMatrix(x)):

	s  :=  z(n+m+1..2*n+m):
	S  :=  Matrix(DiagonalMatrix(s)):

	e  :=  Vector(n,1):

	v1 := (X.s - MuBar(z,n,m).e):
	v2 := (X.ds+dX.s -  ( (DotProduct(x,ds) + DotProduct(s,dx)) / n ).e ):
	v3 := dX.ds-MuBar(dz,n,m).e:


	a0 := DotProduct(v1,v1) - (theta_prime.MuBar(z,n,m))^2:

	a1 := 2*( DotProduct(v1,v2) - (theta_prime^2) * MuBar(z,n,m)*((DotProduct(x,ds) + DotProduct(s,dx))/n) ):
	
	a2 := DotProduct(v2,v2)+ 2*DotProduct(v1,v3) - (theta_prime^2)* (2*MuBar(z,n,m)*MuBar(dz,n,m) +  ( (DotProduct(x,ds) + DotProduct(s,dx))/n  )^2 ) :

	a3 := 2*DotProduct(v2,v3)   - 2 * theta_prime * MuBar(dz,n,m) * ((DotProduct(x,ds) + DotProduct(s,dx))/n)  :
	
	a4 := DotProduct(v3,v3) - (theta_prime*MuBar(dz,n,m))^2:

	Sol := fsolve(a0 + a1*rho + a2*rho^2 + a3*rho^3 + a4*rho^4, rho, fulldigits, 0 ..infinity):
	
	if nops([Sol]) = 0 then print("alpha: Makaynch SOLOTION :( :(") end if;
	if nops([Sol]) = 1 then return Sol end if;
	return Sol[1];

end	proc:




GetNextPoint := proc(n,m,A,z,theta,theta_prime)

	########## Prediction ###########
	local d,rho,z_prime, d_prime,z_plus;
	d := GetNewtonDirection(n,m,A,z,0):

	rho := GetMaxAlpha(n,m,z,d,theta_prime);
	
	z_prime := z + rho.d:
	
	print("Prediction "):
	
	print("rho"):
	print(rho):
	



	print("z_prime"):
	print(z_prime(n-1..n)):
	

	########## Correction ############


	d_prime := GetNewtonDirection(n,m,A,z_prime,MuBar(z_prime,n,m)):
	z_plus := z_prime + d_prime:

	print("Correction "):
	print(z_plus(n-1..n)):

	return evalf(z_plus,1000000);

end proc:
	



SolvePL := proc(n,m,A,z,theta,theta_prime,N)
	local zz, L, z1;
	L := [z(1..n)]:
	zz:= z;

	for i from 1 to N do
			z1 := GetNextPoint(n,m,A,zz,theta,theta_prime):
			zz := z1:
			L := [op(L),zz(1..n)]:
	end do:
end proc:










t:=2;
r:=3:

n := 5*r+1:
m := 3*r+1:
lambda:= 3/2:
theta := 0.25:
theta_prime := 0.5:



A := AGenerator(r, t):
At := Transpose(A):
b := BGenerator(r,t):
c := convert(<1,Vector(2*r-1,0)>,Vector):




A_eq := EquatlityConstraintsAGenerator(r,t):
A_eqt:= Transpose(A_eq):
b_eq := b:
c_eq := Vector[column](<c,Vector(3*r+1,0)>):


alpha := <2^6,2^6,2^5,2^5,2^4,2^4,2^3,2^3,2^2,2^2,2^1,2^1,2^0,2^0>/(2^8): 
l := 2.01:
beta := <l^6, l^6, l^5, l^5, l^5,l^4, l^4, l^4, l^3, l^3, l^3, l^2, l^2, l^2, l, l, l, l, l>/(l^7):




Y:= Tropical_y(lambda,r):
X:= Tropical_x(lambda,r):

f := x -> t^x:

v_y:= convert(map(f, Y),Vector):


v_x:= convert(map(f, X),Vector):


lifted_x := (alpha *~ v_x):
lifted_w := (b - A.lifted_x):
lifted_y := (beta *~ v_y):
lifted_s := (Multiply(At,lifted_y) + c):




x_eq := convert(<lifted_x,lifted_w>,Vector[column]):
y_eq := convert(<-lifted_y>,Vector[column]):
s_eq := convert(<lifted_s,lifted_y>,Vector[column]):


c_eq := convert(<c,Vector[column](3*r+1,0)>,Vector[column]):
z_eq := convert(<x_eq,y_eq,s_eq>,Vector[column]):








IsAdmissible := proc(z, A,b,c, n,m)
	local x,y,s, b1,b2,b3;
	x  :=  z(1..n);
	y  := z(n+1 .. n+m);
	s  :=  z(n+m+1..2*n+m);

	b1:= Equal(Multiply(A,x) , b);
	b2:= Equal( Multiply(Transpose(A) , y) + s , c);
	b3:= true;
	for i from 1 to n do 
		if evalf(x(i))<=0 then b3:= false end if;
		if evalf(s(i))<=0 then b3:= false end if;
		end do;
	return (b1 and b2 and b3);	
end proc;


IsInV := proc(theta,z, A,b,c, n,m)
	
	local x,y,s, b1,muBar;
	x  :=  z(1..n);
	y  := z(n+1 .. n+m);
	s  :=  z(n+m+1..2*n+m);

	b1:= IsAdmissible(z,A,b,c,n,m);
	muBar := MuBar(z,n,m);
	return (b1 and evalb( Norm( Multiply( Matrix(DiagonalMatrix(x)), s) - muBar * Vector(n,1) ,2) <= theta * muBar));
end proc;









































































##################################################
#################
mu  := evalf(MuBar(z_eq,n,m));
evalf(Norm( DiagonalMatrix(x_eq).s_eq - mu.Vector(n,1)  ,2 ) / mu);


d_p := GetNewtonDirection(n,m,A_eq,z_eq,mu):

z   := evalf(z_eq + d_p):

x  :=  z(1..n):
y  := z(n+1 .. n+m):
s  :=  z(n+m+1..2*n+m):


z_eq :=z: 
x_eq :=x:
y_eq :=y:
s_eq :=s:





##################################################
#################
mu  := evalf(MuBar(z_eq,n,m));
evalf(Norm( DiagonalMatrix(x_eq).s_eq - mu.Vector(n,1)  ,2 ) / mu);


d_p := GetNewtonDirection(n,m,A_eq,z_eq,mu):

z   := evalf(z_eq + d_p):

x  :=  z(1..n):
y  := z(n+1 .. n+m):
s  :=  z(n+m+1..2*n+m):


z_eq :=z: 
x_eq :=x:
y_eq :=y:
s_eq :=s:




##################################################
#################
mu  := evalf(MuBar(z_eq,n,m));
evalf(Norm( DiagonalMatrix(x_eq).s_eq - mu.Vector(n,1)  ,2 ) / mu);


d_p := GetNewtonDirection(n,m,A_eq,z_eq,mu):

z   := evalf(z_eq + d_p):

x  :=  z(1..n):
y  := z(n+1 .. n+m):
s  :=  z(n+m+1..2*n+m):


z_eq :=z: 
x_eq :=x:
y_eq :=y:
s_eq :=s:




##################################################
#################
mu  := evalf(MuBar(z_eq,n,m));
evalf(Norm( DiagonalMatrix(x_eq).s_eq - mu.Vector(n,1)  ,2 ) / mu);


d_p := GetNewtonDirection(n,m,A_eq,z_eq,mu):

z   := evalf(z_eq + d_p):

x  :=  z(1..n):
y  := z(n+1 .. n+m):
s  :=  z(n+m+1..2*n+m):


z_eq :=z:

x_eq :=x;
y_eq :=y;
s_eq :=s;



###########################################################
###########################################################
###########################################################
##########################################################
###########################################################
###########################################################
##########################################################
###########################################################
###########################################################
##########################################################
###########################################################
###########################################################


z := z_eq:

z(n-1..n);


L:=SolvePL(n,m,A_eq,z,theta,theta_prime,20):









