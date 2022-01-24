function Predcv()
A=textread('knowndiseasemicrobeinteraction.txt');
% nd:the number of diseases
% nm:the number of microbe
% pp:the number of known diseae-microbe associations
nd=max(A(:,1)); 
nm=max(A(:,2));
[pp,]=size(A);
%interaction: adjacency matrix for the disease-microbe association network
%interaction(i,j)=1 means microbe j is related to disease i
for i=1:pp
        interaction(A(i,1),A(i,2))=1;
end
save interaction interaction;

%implement leave-one-out cross validation
for cv=1:pp 
    % obtain training sample
    load interaction;
    interaction(A(cv,1),A(cv,2))=0;

    [nd,nc] = size(interaction); % nd:diseases number, nc:microbes number

circSim01  = GSM( interaction'); % microbe gaussian
disSim01  = GSD( interaction' ); % disease gaussian

disSim02  = cosSim( interaction); % disease cosine
circSim02  = cosSim( interaction');   % microbe cosine

w1 = 0.6;
w2 = 0.3;
circRNA_sim  = LKF(w1,circSim01,circSim02); % LKF microbe gaussian and cosine
dis_sim  = LKF(w2,disSim01,disSim02); % LKF disease gaussian and cosine 

F = NCPLDA(circRNA_sim,dis_sim,interaction');

F = F';
    
% obtain the score of tested  disease-microbe interaction
finalscore=F(A(cv,1),A(cv,2));
% make the score of seed  disease-microbe interactions as zero
for i=1:nd
    for j=1:nm
        if interaction(i,j)==1
           F(i,j)=-10000;
        end
    end
end

% obtain the position of tested disease-microbe interaction as variable globalposition(1,cv),
[ll1,mm1]=size(find(F>=finalscore));
[ll2,mm2]=size(find(F>finalscore));
positionPred(1,cv)=ll2+1+(ll1-ll2-1)/2;

end
 save('positionPred.mat','positionPred');   
end


        
        
        
    
   




