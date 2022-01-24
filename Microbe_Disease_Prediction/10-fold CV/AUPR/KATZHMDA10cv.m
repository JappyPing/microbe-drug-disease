function KATZHMDA10cv(gamadd,gamall,beta,k)
%KATZHMDA10cv(1,1,0.01,2)
%predict disease-related microbe based on KATZHMDA  in the term of 5-fold cross validation
%A: Binary relations between disease and microbe, 1st column:disease, 2nd column:microbe
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
    
   
%calculate gamad for Gaussian kernel calculation
    for i=1:nd
        sd(i)=norm(interaction(i,:))^2;
    end
    gamad=nd/sum(sd')*gamadd;
    
 %calculate gamal for Gaussian kernel calculation
        for i=1:nm
        sl(i)=norm(interaction(:,i))^2;
    end
    gamal=nm/sum(sl')*gamall;
    
    %calculate Gaussian kernel for the similarity between disease: kd
    for i=1:nd
        for j=1:nd
    pkd(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
        end
    end
    
    %calculate Gaussian kernel for the similarity between microbe: km
        for i=1:nm
        for j=1:nm
    km(i,j)=exp(-gamal*(norm(interaction(:,i)-interaction(:,j)))^2);
        end
        end 

        for i=1:nd
for j=1:nd
    kd(i,j)=1/(1+exp(-15*pkd(i,j)+log(9999)));
end
        end
    
if k==2
F=beta*interaction'+(beta^2)*(km*interaction'+interaction'*kd);
else if k==3
F=beta*interaction'+(beta^2)*(km*interaction'+interaction'*kd)+(beta^3)*(interaction'*interaction*interaction'+km*km*interaction'+km*interaction'*kd+interaction'*kd*kd);
else if k==4
F=beta*interaction'+(beta^2)*(km*interaction'+interaction'*kd)+(beta^3)*(interaction'*interaction*interaction'+km*km*interaction'+km*interaction'*kd+interaction'*kd*kd)+(beta^4)*(km*km*km*interaction'+interaction'*interaction*km*interaction'+km*interaction'*interaction*interaction'+interaction'*kd*interaction*interaction')+(beta^4)*(interaction'*interaction*interaction'*kd+km*km*interaction'*kd+km*interaction'*kd*kd+interaction'*kd*kd*kd);
end
end
end

     F=F';

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
positionKATZ(1,T)=ll2+1+(ll1-ll2-1)/2;
T=T+1;
end

end
save('positionKATZ.mat','positionKATZ');  

end



        
        
        
    
   



