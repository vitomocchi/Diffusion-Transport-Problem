function Uh=FDforward(alfa,h,mu,n,u0,uL)

    diag = (2*mu/(h^2)-alfa/h)*ones(1,n-2);
    offdiaginf = (-mu/(h^2))*ones(1,n-2); % spdiags ignora l'ultimo valore.
    offdiagsup = (-mu/(h^2)+alfa/h)*ones(1,n-2); % spdiags ignora il primo valore.

    % Create the sparse matrix using spdiags().
    A = spdiags([offdiaginf', diag', offdiagsup'],-1:1, n-2, n-2);

    b=ones(1,n-2); b=b'; % known term (column).
    b(1)=b(1)+u0/(h^2)+alfa/(2*h)*u0;  % modifies with Dirichlet BC in 0.
    b(end)=b(end)+uL*mu/(h^2)-alfa/(2*h)*uL; % modifies with Dirichlet BC in 1.

    Uh=A\b; % dim Uh: ((n-2)x1)
    Uh=[u0;Uh;uL];

end