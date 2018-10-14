with(LinearAlgebra):
with(linalg):
with(SolveTools):
with(plots):
with(ArrayTools):




AGenerator := proc(r,t);
	local L;
	 L := [Concatenate(2, Vector[row](<1,0>) ,Vector[row](2*r-2,0) ),
	       Concatenate(2, Vector[row](<0,1>) ,Vector[row](2*r-2,0) )];



	



	
	for i from 1 to N do

			z1 := GetNextPoint(n,m,A,z,theta,theta_prime):
			z := z1:
			L := [op(L),z(1..n)]:

	end do:
end proc:





