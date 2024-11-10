function L=liouville_operator(H)
% L=liouville_operator(H)

N=size(H,1);

% implementation using kronecker products and sparse matrices, incompatible
% with multi-dimensional manipulations
L=spalloc(N^2,N^2,2*N*nnz(H)+N^2);

H_diag=diag(H);
H_off=H-diag(H_diag); 
% contribution from the diagonal terms of the Hamiltonian
[i1,i2]=meshgrid(1:N);

% L(1:N^2+1:N^4)=H_diag(i2)-H_diag(i1);
L=L+spdiags(H_diag(i2(:))-H_diag(i1(:)),0,N^2,N^2);

% contribution from the off-diagonal terms of the Hamiltonian
L=L+kron(speye(N),H_off)-kron(H_off',speye(N));
