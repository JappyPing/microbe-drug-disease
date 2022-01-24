function [ rank_result,known_assoccition_rank ] = Rank( result,interaction,microbes,diseases )

    [rows,cols]=size(interaction);
    num_ones=zeros(cols,1);
    for i=1:cols
       num_ones(i,1)=nnz(interaction(:,i)); 
    end
    num=rows-min(num_ones);
    rank_result=cell(num+1,cols);  
    known_assoccition_rank=zeros(size(interaction)); 
    for i=1:cols
       idx=find(interaction(:,i));  
       [~,idx_sort]=sort(result(:,i),'descend');
       
       known_assoccition_rank(:,i)=ismember(idx_sort,idx);  
       
       for j=1:length(idx)  
          del_idx= (idx(j,1)==idx_sort);
          idx_sort(del_idx,:)=[];
       end
       rank_result(1,i)=diseases(i,1);  
       for k=1:length(idx_sort)
          rank_result(k+1,i)=microbes(idx_sort(k,1)); 
       end
        
    end

end

