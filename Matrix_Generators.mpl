Lifter := proc(r,t,lambda,l_alpha,l_beta)
	
	
	alpha := AlphaGenerator(r,l_alpha):
	beta  := BetaGenerator(r,l_beta):

	A_ := AGenerator(r,t):
	At_:= Transpose(A_):
	b_ := BGenerator(r,t):
	c_ := convert(<1,Vector(2*r-1,0)>,Vector):

	f := x -> t^x:

	x_trop:= Tropical_x(lambda,r):
	y_trop:= Tropical_y(lambda,r):

	v_x:= convert(map(f, x_trop),Vector):
	v_y:= convert(map(f, y_trop),Vector):


	lifted_x := (alpha *~ v_x):
	lifted_w := (b_ - A_.lifted_x):
	lifted_y := (beta *~ v_y):
	lifted_s := (Multiply(At_,lifted_y) + c_):


	x_eq := convert(<lifted_x,lifted_w>,Vector[column]):
	y_eq := convert(<-lifted_y>,Vector[column]):
	s_eq := convert(<lifted_s,lifted_y>,Vector[column]):

	z_eq := convert(<x_eq,y_eq,s_eq>,Vector[column]):


end proc:







PointInitial:= proc(n,m,z_eq,A,theta)
	
	local x_eq,s_eq,y_eq ,z,z_prime,dz,mu,Cr,i:


	z := copy(z_eq):

	x  :=  convert(z(1..n),Vector[column]):
	y  :=  convert(z(n+1 .. n+m),Vector[column]):
	s  :=  convert(z(n+m+1..2*n+m),Vector[column]):

	
	mu  := evalf(MuBar(z,n,m)):
	
	Cr := evalf(Norm( DiagonalMatrix(x).s - mu.Vector[column](n,1) ,2 ) / mu):



	while(Cr >= 0.025 ) do

		dz := GetNewtonDirection(n,m,A,z,mu):

		z_prime    := z + dz:
		
		z := z_prime:

		x  :=  convert(z(1..n),Vector[column]):
		y  :=  convert(z(n+1 .. n+m),Vector[column]):
		s  :=  convert(z(n+m+1..2*n+m),Vector[column]):

		mu  := evalf(MuBar(z,n,m)):
		Cr := evalf(Norm( DiagonalMatrix(x).s - mu.Vector[column](n,1) ,2 ) / mu):

	end do:

	return z_eq:

end proc:
