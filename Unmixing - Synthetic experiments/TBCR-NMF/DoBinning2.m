function X = DoBinning2(I)
[m,n] = size(I);
X = zeros(m/2, n/2);

for i = 1:m/2
    for j = 1:n/2
        X(i,j) = (I(2*i-1, 2*j-1)+I(2*i, 2*j-1)+I(2*i-1, 2*j)+I(2*i, 2*j))/4;
    end
end

end