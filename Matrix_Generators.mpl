with(LinearAlgebra):
with(linalg):
with(SolveTools):
with(plots):
with(ArrayTools):



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

