function NBLPIHMDA10cv(gamadd,gamall,p,k1,k2,b)
%NBLPIHMDA10cv(1,1,0.2,5,5,0.75)
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
   list_d=interaction;
    list_m=interaction;
    interaction1=interaction; 
   
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
    kd(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
        end
    end
    
    %calculate Gaussian kernel for the similarity between microbe: km
        for i=1:nm
        for j=1:nm
    km(i,j)=exp(-gamal*(norm(interaction(:,i)-interaction(:,j)))^2);
        end
        end 

    
 km1=km-diag(diag(km));
    for i=1:nm
        kmn=km1(i,:);
        [value,]=sort(kmn,'ascend');
        a=value(k2);
        index=find(kmn<=a);
        for j=1:length(index)
             km1(i,index(j))=0;
        end
        km1(i,:)=km1(i,:)./sum(km1(i,:));
    end

kd1=kd-diag(diag(kd));
    for i=1:nd
        kdn=kd1(i,:);
        [value,xi]=sort(kdn,'ascend');
        a=value(k1);
        index=find(kdn<=a);
        for j=1:length(index)
             kd1(i,index(j))=0;
        end
         kd1(i,:)=kd1(i,:)./sum(kd1(i,:));
    end

    count1=0;
    while true
        interaction_d=p*kd1*interaction1+interaction1*(1-p);
        list_d=list_d+interaction_d;
        C1=interaction_d-interaction1;
        C1=sum(sum(abs(C1)));
        count1=count1+1;
        interaction1=interaction_d;
        if C1<1e-12
            break
        end
    end
    save list_d list_d;
    
    count2=0;
    interaction1=interaction;
    while true
        interaction_m=(p*km1*interaction1')'+interaction1*(1-p);
        list_m=list_m+interaction_m;
        C2=interaction_m-interaction1;
        C2=sum(sum(abs(C2)));
        count2=count2+1;
        interaction1=interaction_m;
        if C2<1e-12
            break
        end
    end
    save list_m list_m;
  
    F=(1-b)*list_d+b*list_m;
    save F F;

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
positionNBLPI(1,T)=ll2+1+(ll1-ll2-1)/2;
T=T+1;
end

end
save('positionNBLPI.mat','positionNBLPI');  

end



        
        
        
    
   



