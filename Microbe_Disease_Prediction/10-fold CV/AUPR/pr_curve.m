load prKATZ.mat
load reKATZ.mat
load prNBLPI.mat
load reNBLPI.mat
load prNTS.mat
load reNTS.mat
load prBi.mat
load reBi.mat
load prBPNN.mat
load reBPNN.mat
load prPred.mat
load rePred.mat
load prLRLS.mat
load reLRLS.mat
load interaction
[n,m] = size(interaction);
%%
plot(reNTS,prNTS,':','LineWidth',1.4);
hold on;
plot(reKATZ,prKATZ,':y','LineWidth',1.4);
hold on;
plot(reNBLPI,prNBLPI,':r','LineWidth',1.4);
hold on;
plot(reBi,prBi,':g','LineWidth',1.4);
hold on;
plot(reBPNN,prBPNN,':m','LineWidth',1.4);
hold on;
plot(rePred,prPred,':c','LineWidth',1.4);
hold on;
plot(reLRLS,prLRLS,':k','LineWidth',1.4);
%%
p1=flip(prKATZ);
r1=flip(reKATZ);
area1(1,1)=p1(1,1)*r1(1,1)/2;
for k=2:10939
    area1(1,k)=[p1(1,k-1)+p1(1,k)]*[r1(1,k)-r1(1,k-1)]/2;
end
auprKATZ=sum(area1);
save auprKATZ auprKATZ;
%%
p2=flip(prNBLPI);
r2=flip(reNBLPI);
area2(1,1)=p2(1,1)*r2(1,1)/2;
for k=2:10939
    area2(1,k)=[p2(1,k-1)+p2(1,k)]*[r2(1,k)-r2(1,k-1)]/2;
end
auprNBLPI=sum(area2);
save auprNBLPI auprNBLPI;
%%
p3=flip(prNTS);
r3=flip(reNTS);
area3(1,1)=p3(1,1)*r3(1,1)/2;
for k=2:10939
    area3(1,k)=[p3(1,k-1)+p3(1,k)]*[r3(1,k)-r3(1,k-1)]/2;
end
auprNTS=sum(area3);
save auprNTS auprNTS;
%%
p4=flip(prBi);
r4=flip(reBi);
area4(1,1)=p4(1,1)*r4(1,1)/2;
for k=2:10939
    area4(1,k)=[p4(1,k-1)+p4(1,k)]*[r4(1,k)-r4(1,k-1)]/2;
end
auprBi=sum(area4);
save auprBi auprBi;

%%
p5=flip(prBPNN);
r5=flip(reBPNN);
area5(1,1)=p5(1,1)*r5(1,1)/2;
for k=2:10939
    area5(1,k)=[p5(1,k-1)+p5(1,k)]*[r5(1,k)-r5(1,k-1)]/2;
end
auprBPNN=sum(area5);
save auprBPNN auprBPNN;
%%
p6=flip(prPred);
r6=flip(rePred);
area6(1,1)=p6(1,1)*r6(1,1)/2;
for k=2:10939
    area6(1,k)=[p6(1,k-1)+p6(1,k)]*[r6(1,k)-r6(1,k-1)]/2;
end
auprPred=sum(area6);
save auprPred auprPred;
%%
p7=flip(prLRLS);
r7=flip(reLRLS);
area7(1,1)=p7(1,1)*r7(1,1)/2;
for k=2:10939
    area7(1,k)=[p7(1,k-1)+p7(1,k)]*[r7(1,k)-r7(1,k-1)]/2;
end
auprLRLS=sum(area7);
save auprLRLS auprLRLS;

xlabel('Recall ');
ylabel('Precision');
title('10-fold-PR curve');
legend({'NTSMDA(aupr=0.6728)','KATZHMDA(aupr=0.6131)','NBLPIHMDA(aupr=0.6784)','BiRWMP(aupr=0.6557)','BPNNHMDA(aupr=0.7293)','HMDA-Pred(aupr=0.4808)','LRLSHMDA(aupr=0.6277)'},'FontSize',16,'Location','NorthEast')