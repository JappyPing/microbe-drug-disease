function [ result ] = cosSim( data )
%COSSIM Summary of this function goes here
%   Detailed explanation goes here
   
rows=size(data,1);
result=zeros(rows,rows);

    for i=1:rows
        
        for j=1:i
            
            if (norm(data(i,:))*norm(data(j,:))==0)
                
                result(i,j)=0;
                
            else
                result(i,j)=dot(data(i,:),data(j,:))/(norm(data(i,:))*norm(data(j,:)));
                
            end
            
            result(j,i)=result(i,j);
        end
        
        
    end

end

