function Predcv(gamadd,gamall,p,k1,k2,b)
%cv(1,1,0.2,5,5,0.75)
A=textread('knowndiseasemicrobeinteraction.txt');
% nd:the number of diseases
% nm:the number of microbe
% pp:the number of known diseae-microbe associations
nd=max(A(:,1)); 
nm=max(A(:,2));
[pp,qq]=size(A);
%interaction: adjacency matrix for the disease-microbe association network
%interaction(i,j)=1 means microbe j is related to disease i
for i=1:pp
        interaction(A(i,1),A(i,2))=1;
end

save interaction interaction;

%implement 10-fold cross validation
x=randperm(pp)';
T=1;

for cv=1:10
load interaction interaction;
    if cv<10
        B=A(x((cv-1)*floor(pp/10)+1:floor(pp/10)*cv),:);
% obtain training sample
    for i=1:floor(pp/10)
        interaction(B(i,1),B(i,2))=0;
    end
    else B=A(x((cv-1)*floor(pp/10)+1:pp),:);
        % obtain training sample
    for i=1:pp-floor(pp/10)*9
        interaction(B(i,1),B(i,2))=0;
    end
    end
    
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

[size1B,size2B]=size(B);
% obtain the score of tested  disease-microbe interaction
for i=1:size1B
finalscore(i,1)=F(B(i,1),B(i,2));
end
% make the score of seed  disease-microbe interactions as zero
for i=1:nd
    for j=1:nm
        if interaction(i,j)==1
           F(i,j)=-10000;
        end
    end
end


for qq=1:size1B
% obtain the position of tested disease-microbe interaction as variable position(1,cv), 
[ll1,mm1]=size(find(F>=finalscore(qq)));
[ll2,mm2]=size(find(F>finalscore(qq)));
position(1,T)=ll2+1+(ll1-ll2-1)/2;
T=T+1;
end

end
save('position.mat','position');  
