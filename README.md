# maths-maxima-code-for-topics

/*1. Create and transform vectors and matrices (the transpose vector (matrix) conjugate transpose of a vector (matrix)).*/
load("eigen")$;
row:read("enter your matrix row:")$;
col:read("enter your matrix column :")$;
Matrix:entermatrix(row,col);
Vector:covect(read("Enter your vector in square bracket []:"));
Transpose_Matrix:transpose(Matrix);
Transpose_Vector:transpose(Vector);


/*2. Generate the matrix into echelon form and find its rank.*/
row:read("enter your matrix row:")$;
col:read("enter your matrix column :")$;
Matrix:entermatrix(row,col);
Echelon:echelon(Matrix);
Rank:rank(Echelon);


/*3. Find cofactors, determinant, adjoint and inverse of a matrix.*/ 
row:read("enter your matrix row:")$;
col:read("enter your matrix column :")$;
Matrix:entermatrix(row,col);
Cofactors:transpose(adjoint(Matrix));
Determinant:determinant(Matrix);
Adjoint:adjoint(Matrix);
Inverse:invert(Matrix);


/*4. Solve a system of Homogeneous and non-homogeneous equations using Gauss elimination method. 
Non-homogeneous equation
-5x-2y+2z=14
3x+1y-1z=-8
2x+2y-1z=-3*/
equation1:read("enter first equation:");
equation2:read("enter second equation:");
equation3:read("enter third equation:");
Matrix:augcoefmatrix([equation1,equation2,equation3],[x,y,z]);
Echelon:echelon(Matrix)$;
Solution:linsolve([Echelon[1][1]*x+Echelon[1][2]*y+Echelon[1][3]*z=Echelon[1][4],Echelon[2][2]*y+Echelon[2][3]*z=Echelon[2][4], Echelon[3][3]*z=Echelon[1][4]], [x,y,z]);

/*Homogeneous equation
x+y+7z=0
2x+3y+17z=0
x+2y+z=0*/
equation1:read("enter first equation:");
equation2:read("enter second equation:");
equation3:read("enter third equation:");
Matrix:augcoefmatrix([equation1,equation2,equation3],[x,y,z]);
Echelon:echelon(Matrix)$;
Solution:linsolve([Echelon[1][1]*x+Echelon[1][2]*y+Echelon[1][3]*z=0,Echelon[2][2]*y+Echelon[2][3]*z=0, Echelon[3][3]*z=0], [x,y,z]);


/*5. Solve a system of Homogeneous equations using the Gauss Jordan method.
2x+2y+2z=0
x+2y+z=0
3x+y-z=0*/
Matrix:matrix([2,2,2,0],[1,2,1,0],[3,1,-1,0]);
Echelon:echelon(Matrix);
Oper:row(Echelon,1)-(row(Echelon,2)+row(Echelon,3));
Rref:addrow(Oper,row(Echelon,2),row(Echelon,3));
Solution:linsolve([x=0, y=0, z=0], [x,y,z]);


/*6. Generate basis of column space, null space, row space and left null space of a matrix space. */
row:read("enter your matrix row:")$;
col:read("enter your matrix column :")$;
mat:entermatrix(row,col);
ColumnSpace:columnspace(mat);
NullSpace:nullspace(mat);
RowSpace:columnspace(transpose(mat));
leftNullSpace:nullspace(transpose(mat));


/*7. Check the linear dependence of vectors. Generate a linear combination of given vectors of Rn/ matrices of the same size and find the transition matrix of given matrix space. */
load("eigen")$;
vector1:covect(read("Enter your first vector in square bracket []:"))$;
vector2:covect(read("Enter your second vector in square bracket []:"))$;
vector3:covect(read("Enter your third vector in square bracket []:"))$;
mat:addcol(vector1,vector2,vector3);
Echelon:echelon(mat)$;
Rank:rank(mat);
if Rank=3 then("vector are linear independence")else("vector are linear dependence");
for i: 1 thru 3 do
        poly:print(mat[i][1]*"x"+mat[i][2]*"y"+mat[i][3]*"z"=0);


/*8. Find the orthonormal basis of a given vector space using the Gram-Schmidt orthogonalization process. */
load("eigen")$;
x: entermatrix (3,4);
y: gramschmidt (x);
a:y[1][4];
f1:(y[1][1]^2+y[1][2]^2+y[1][3]^2+y[1][4]^2)^0.5$;
row1:y[1]/%$;
f2:(y[2][1]^2+y[2][2]^2+y[2][3]^2+y[2][4]^2)^0.5$;
row2:y[2]/%$;
f3:(y[3][1]^2+y[3][2]^2+y[3][3]^2+y[3][4]^2)^0.5$;
row3:y[3]/%$;
Orthonormal_Basis:addrow(matrix(row1),matrix(row2),matrix(row3));


/*9. Check the diagonalizable property of matrices and find the corresponding eigenvalue and verify the Cayley- Hamilton theorem. */
load("eigen")$;
mat:matrix([-4,7,1,4],[6,-16,-3,-9],[12,-27,-4,-15],[-18,43,7,24]);
poly:expand(charpoly(mat,x));
Eigenvalues:eigenvalues(mat);
Eigenvectors:eivects(mat)[2];
P:determinant(addrow(matrix(Eigenvectors[1][1]),matrix(Eigenvectors[2][1]),matrix(Eigenvectors[2][2]),matrix(Eigenvectors[3][1])));
if P = 0 then("given matrix is not diagonalizable")
else ("given matrix is diagonalizable");
sol:-mat^^3+19*mat^^2-57*mat-9*ident(3);


/*10. Application of Linear algebra: Coding and decoding of messages using nonsingular matrices.   then decode it. */
"Enter your message in matrix form (numbers should be equal to alphabetical positions)";
"{A:1,B:2,C:3,D:4,E:5,F:6,G:7,H:8,I:9,J:10,K:11,L:12,M:13,N:14,O:15,P:16,Q:17,R:18,S:19,T:20,U:21,V:22,W:23,X:24,Y:25,Z:26, :27}";
row:read("enter your matrix row:")$;
col:read("enter your matrix column :")$;
h [i, j] :=2/ (i +6* j - 1)$;
mat_de:genmatrix(h,row,col)$;
mat_input:entermatrix(row,col);                           
decode_message:mat_de*mat_input;
coded_message:decodeMatrix/mat_de;


/*11. Compute Gradient of a scalar field. */
load("vect")$;
F(x,y,z):=read("Enter your Function (like x*y*exp(z^2)):");
Grad:grad(F(x,y,z));
Express_grad:express(Grad);                                     
Solution:ev(Express_grad,diff);


/*12. Compute Divergence of a vector field. */
load("vect")$;
F(x,y,z):=read("Enter your Function (like [x,2*y,exp(z*y)]):");
Div:div(F(x,y,z));                           
Express_div:express(Div);
Solution_div:ev(Express_div,diff);


/*13. Compute Curl of a vector field.*/
load("vect")$;
F(x,y,z):=read("Enter your Function (like [x,2*y,exp(z*y)]):");
Curl:curl(F(x,y,z));                      
Express_curl:express(Curl);
Solution_curl:ev(Express_curl,diff);
