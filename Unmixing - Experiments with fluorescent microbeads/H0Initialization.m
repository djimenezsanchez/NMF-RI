function [H0] = H0Initialization(Y,SepChannel,nFl)

% H0 Initialization depending on number of fluorochromes 
if size(Y,1) == nFl
    H0 = Y;
else
    
    COUNT = 1;
    for i=SepChannel        
        if COUNT==1
        H0(COUNT,:) = sum(Y(1:i,:),1);
        else if COUNT>1
        H0(COUNT,:) = sum(Y(SepChannel(COUNT-1)+1:SepChannel(COUNT),:),1);
            end  
        end
        COUNT = COUNT+1; 
    end
    H0(COUNT,:) = sum(Y(SepChannel(COUNT-1)+1:end,:));
end 
    