function overallauc=positiontooverallauc()
load positionKATZ.mat;
load interaction;
[n,m]=size(interaction);
sID=textread('knowndiseasemicrobeinteraction.txt');
[pp,qq]=size(sID);

for i=1:pp
    if i<floor(pp/10)*9+1

knowndt(i)=pp-floor(pp/10);
        unknown(i)=n*m-knowndt(i);
    else knowndt(i)=floor(pp/10)*9;
         unknown(i)=n*m-knowndt(i);
    end
end

for i=1:pp
    for j=1:m*n
        if j<unknown(i)+1
            rank(i,j)=j;
        else rank(i,j)=unknown(i)+1;
        end
    end
end

for k=1:m*n-floor(pp/10)*9
    tp=0;
    for t=1:pp
        if positionKATZ(1,t)<=k
            tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp;
    if k<m*n-pp+floor(pp/10)+1
    fp=k*pp-tp;
    else fp=floor(pp/10)*9*(m*n-pp+floor(pp/10))+(pp-floor(pp/10)*9)*k-tp;
    end

     fpr(1,k)=fp/(floor(pp/10)*9*(m*n-pp+floor(pp/10)-1)+(pp-floor(pp/10)*9)*(m*n-floor(pp/10)*9-1));
end
plot(fpr,tpr);
clear area;
area(1,1)=tpr(1,1)*fpr(1,1)/2;
for k=2:m*n-floor(pp/10)*9
    area(1,k)=[tpr(1,k-1)+tpr(1,k)]*[fpr(1,k)-fpr(1,k-1)]/2;
end
overallauc=sum(area);
end
          


