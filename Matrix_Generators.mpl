with(LinearAlgebra):
with(linalg):
with(SolveTools):
with(plots):
with(ArrayTools):


AGenerator := proc(r,t)

	local L,Line1,Line2;

	Line1 := convert(Concatenate(2, Vector[row](<1,0>) ,Vector[row](2*r-2,0) ),list);
	Line2 := convert(Concatenate(2, Vector[row](<0,1>) ,Vector[row](2*r-2,0) ),list);
	
	L := [Line1, Line2];
	
	for i from 1 to r-1 do

		Line1 := convert(Concatenate(2, Vector[row](2*(i-1),0), Vector[row](<-t,0,1>), Vector[row](2*(r-i)-1,0) ),list);
		Line2 := convert(Concatenate(2, Vector[row](2*(i-1),0), Vector[row](<0,-t,1>), Vector[row](2*(r-i)-2,0) ),list);
			
			
		L := [op(L),Line1]:
		L := [op(L),Line2]:

	end do;
	

	for i from 1 to r-1 do
		Line1 := convert(Concatenate(2, Vector[row](2*(i-1),0), Vector[row](<-(t^(1-(1/(2^i)))),-(t^(1-(1/(2^i)))),0,1>), Vector[row](2*(r-i-2),0) ),list);
		
		L := [op(L),Line1]:
		
	end do;

	
	Line1 := convert(Concatenate(2, Vector[row](2*r-2,0), Vector[row](<-1,0>)),list);
	Line2 := convert(Concatenate(2, Vector[row](2*r-2,0), Vector[row](<0,-1>)),list);
	
	L := [op(L),Line1]:
	L := [op(L),Line2]:

	
	return Matrix(3*r+1,2*r,L);
	
	
end proc:


bGenerator := proc(r,t)

	return Concatenate(1,Vector[column](<t^2,t>),Vector[column](3*r-1,0));
	
end proc:





A := AGenerator(6,t):
A;

b := bGenerator(6,t):
b;

