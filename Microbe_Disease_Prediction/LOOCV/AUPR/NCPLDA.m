function [ncp] = NCPLDA(MiFunSim,DiPheSim,mi2diNetwork)
	
	[MnRow,~] = size(MiFunSim);
	[~,DnCol] = size(DiPheSim);
	[RnRow,RnCol] = size(mi2diNetwork);
	%ncp_m = zeros(RnRow,RnCol)+10^-30;
	%ncp_d = zeros(RnRow,RnCol)+10^-30;
    ncp_m = zeros(RnRow,RnCol);
	ncp_d = zeros(RnRow,RnCol);
	ncp =zeros(RnRow,RnCol);
    
    for i = 1:RnRow
		for j = 1:RnCol
			
            if (mi2diNetwork(i,j)==0)
                mi2diNetwork(i,j)=mi2diNetwork(i,j)+10^-30;
            end
            
        end
   end
    
    	for j = 1:RnCol
        	ncp_m(:,j) = (MiFunSim*mi2diNetwork(:,j))/norm(mi2diNetwork(:,j));
    	end

	for i = 1:RnRow
		ncp_d(i,:) = (mi2diNetwork(i,:)*DiPheSim)/norm(mi2diNetwork(i,:));
	end
			
	for i = 1:MnRow
		for j = 1:DnCol
			ncp(i,j) = (ncp_m(i,j) + ncp_d(i,j))/(norm(MiFunSim(i,:))+norm(DiPheSim(:,j)));
		end
	end

	
end

