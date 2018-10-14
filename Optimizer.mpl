with(LinearAlgebra):
with(linalg):
with(SolveTools):
with(plots):



MuBar := proc(z,n,m)
	local x,s;
	x  :=  z(1..n);
	s  :=  z(n+m+1..2*n+m);
	return DotProduct(Transpose(x),s) / n
end proc;


IsAdmissible := proc(z, A,b,c, n,m)
	local x,y,s;
	x  :=  z(1..n);
	y  := z(n+1 .. n+m);
	s  :=  z(n+m+1..2*n+m);

	local b1,b2,b3;
	b1:= Equal(Multiply(A,x) , b);
	b2:= Equal( Multiply(Transpose(A) , y) + s , c);
	b3:= true;
	for i from 1 to n do 
		if x(i)<=0 then b3:= false end if;
		if s(i)<=0 then b3:= false end if;
		end do;
	return (b1 and b2 and b3);	
end proc;


IsInV := proc(theta,z, A,b,c, n,m)
	
	local x,y,s;
	x  :=  z(1..n);
	y  := z(n+1 .. n+m);
	s  :=  z(n+m+1..2*n+m);

	local b1,muBar;
	b1:= IsAdmissible(z,A,b,c,n,m);
	muBar := MuBar(z,n,m);
	return (b1 and evalb( Norm( Multiply( Matrix(DiagonalMatrix(x)), s) - muBar * Vector(n,1) ,2) <= theta * muBar);)
end proc;


NewtonMatrix := proc(n,m,A,z)
	local x,s;
	x  :=  z(1..n):
	s  :=  z(n+m+1..2*n+m):

	local X,S;
	X := Matrix(DiagonalMatrix(x)):
	S := Matrix(DiagonalMatrix(s)):

	return Matrix(blockmatrix(3,3,[ZeroMatrix(n,n),  Transpose(A), IdentityMatrix(n), A, ZeroMatrix(m,m), ZeroMatrix(m,n), S, ZeroMatrix(n,m), X])):

end proc: 



GetNewtonDirection := proc(n,m,A,z,mu)
	local x,s;
	x  :=  z(1..n):
	s  :=  z(n+m+1..2*n+m):

	local M,X,e,Vect;
	M  := NewtonMatrix(n,m,A,z):

	X      :=  Matrix(DiagonalMatrix(x)):
	e      :=  Vector(n,1):
	Vect   :=  convert(blockmatrix(3,1,[Vector(n,0),Vector(m,0),mu.e - X.s]),Vector):

	return LinearSolve(M,Vect):


end proc:



GetMaxAlpha := proc(n,m,z,dz,theta_prime)

	local dx,dX,ds,dS,x,X,s,S,e,v1,v2,v3,a0,a1,a2,a3,a4;
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

	local S;
	S := fsolve(a0 + a1*alpha + a2*alpha^2 + a3*alpha^3 + a4*alpha^4, alpha, fulldigits, -infinity ..0):
	
	if nops([S]) = 0 then print("alpha: Makaynch SOLOTION :( :(") end if;

	return S[1];

end	proc:



GetNextPoint := proc(n,m,A,z,theta,theta_prime)

	########## Prediction ###########
	local d,alpha,z_prime;
	d := GetNewtonDirection(n,m,A,z,0):
	alpha := GetMaxAlpha(n,m,z,d,theta_prime);

	z_prime := z + alpha.d:

	print("Prediction "):
	
	print("alpha"):
	print(alpha):
	
	print("d"):
	print(d):


	print("z_prime"):
	print(z_prime):
	

	########## Correction ############

	local d_prime,z_plus;
	d_prime := GetNewtonDirection(n,m,A,z_prime,MuBar(z_prime,n,m)):
	z_plus := z_prime + d_prime:

	print("Correction "):
	print(z_plus):

	return evalf(z_plus,1000000);

end proc:
	


SolvePL := proc(n,m,A,z,theta,theta_prime,N);
	
	L := [z(1..n)]:


	for i from 1 to N do

			z1 := GetNextPoint(n,m,A,z,theta,theta_prime):
			z := z1:
			L := [op(L),z(1..n)]:

	end do:
end proc:







A := Matrix(2,4,[[1,0,0,0],[0,1,0,0]]):

n := 4:
m := 2:	
b := Vector(2,1): 
z := <1,1,100,102,  -2,-2,  2,2,1,1>;
theta := 1/4;
theta_prime := 1/2;


SolvePL(n,m,A,z,theta,theta_prime,3);




(*





t := 2;
n := 11;
m := 7;
theta := 1/4;
theta_prime := 1/2;
A_prime := Matrix(7,4,[[1       ,0       ,0  , 0],
				      [0       ,1        ,0  , 0],
				 	  [-t      ,0        ,1  , 0],
				      [0       ,-t       ,1  , 0],
				      [-t^(1/2),-t^(1/2) ,0  , 1],
				      [0       ,0        ,-1 , 0],
				      [0       ,0        ,0  ,-1]]);

A := Matrix(blockmatrix(1,2,[A_prime,IdentityMatrix(m)]));
b := convert(blockmatrix(2,1,[<t^2,t>,Vector(m-2,0)]), Vector );
c := <1,0,0,0,0,0,0,0,0,0,0>;



x_prime := <1,1,1,1>;
s_prime := b - A_prime.x_prime;
x := convert(blockmatrix(2,1,[x_prime,s_prime]),Vector);


y := -<7,7,2,2,2,1,1>;


s := c - Transpose(A).y;


z := convert(blockmatrix(3,1,[x,y,s]),Vector);





dx :=  dz(1..n);
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

S_exact := Polynomial(a0 + a1*alpha + a2*alpha^2 + a3*alpha^3 + a4*alpha^4,alpha, explicit = true):
S_floats := evalf(S_exact):

print(S_floats);

k:= nops(S_exact);

alpha_max :=0:

for i from 1 to k do 

	if type( S_floats[i],'float' ) then

				if S_floats[i]>=alpha_max then 
							alpha_max:= S_exact[i]:
				end if:
	end if:
end do:

alpha_max;


*)

