function SAD = SAD_distance(u,v)
for i =1:size(u,2)
    SAD(i)=abs(acos(u(:,i)'*v(:,i)/(norm(u(:,i))*norm(v(:,i)))));
end
SAD = sum(SAD);
