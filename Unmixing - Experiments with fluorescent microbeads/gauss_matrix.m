function W02 = gauss_matrix(n,m)

if nargin<2
    m=n;
end

    I=[];
    for i=1:n
    J = zeros(1,m);    
    J=gaussmf(1:m,[m/4 i*(m/(n+1))]);
    I = [I;J];
    end    
    W02 = I';
    %W02 = power(2,-( (I-1)*r+1 -J)/r);
    %W02(W02>1) = min(W02(W02>0));
return


