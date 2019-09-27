function W02 = exponential_matrix(n,m)

if nargin<2
    m=n;
end

if n==m
    I = repmat([1:n]',1,n);
    J = I';

    W02 = power(2,-(I-J));
    W02(I<J) = power(2,-n);
    W02(I~=J) = W02(I~=J)/20; 
else
    I = repmat([1:m]',1,n);
    J = repmat([1:n],m,1);
    
    r = (n)/(m);
    
    W02 = power(2,-(I*r -(J/(n+1))*n)*1.2);
    W02(W02>1) =  2.^(1-(W02(W02>1)-1))-1;
    W02(W02<=0) = min(W02(W02>0));
    
end


return

W02(1,1)=n;
for i=2:n
    W02(i,:)=W02(i-1,:)/2;
    W02(i,i)=n;
end
W02=W02+triu(ones(n),1)*min(W02(W02>0))/2;
