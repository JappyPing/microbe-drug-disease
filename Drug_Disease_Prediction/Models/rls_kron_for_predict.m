function candiRankList = rls_kron_for_predict(y,ka,kb,DiseaseIDs,DrugIDs,sigma)
	% Kronecker product Regularized Least Squares for association prediction.
	
	[va,la] = eig(ka);
	[vb,lb] = eig(kb);
	
	l = kron(diag(lb)',diag(la));
	l = l ./ (l + sigma);
	
	m1 = va' * y * vb;
	m2 = m1 .* l;
	y2 = va * m2 * vb';
    %candiRankList=cell(0,0);
	candiRankList = zeros(length(DrugIDs),length(DiseaseIDs));
	%count = 1;
	for i=1:length(DrugIDs)
		for j=1:length(DiseaseIDs)
			%candiRankList(count,1) = DrugIDs(i);
			%candiRankList(count,2) = DiseaseIDs(j);
			%candiRankList(count,3) = num2cell(y2(i,j));
			%count = count + 1;
            candiRankList(i,j)=y2(i,j);
		end
	end
	%[~,IndexRow]=sort(cell2mat(candiRankList(:,3)),'descend');
	%candiRankList = candiRankList(IndexRow,:);
end
