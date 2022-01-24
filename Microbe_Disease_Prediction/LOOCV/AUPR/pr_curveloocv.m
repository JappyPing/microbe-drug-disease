load pr1KATZ.mat
load re1KATZ.mat
load pr1NBLPI.mat
load re1NBLPI.mat
load pr1NTS.mat
load re1NTS.mat
load pr1Bi.mat
load re1Bi.mat
load pr1BPNN.mat
load re1BPNN.mat
load pr1Pred.mat
load re1Pred.mat
load pr1LRLS.mat
load re1LRLS.mat
load interaction
[n,m] = size(interaction);
%t=m*n;
%%
plot(re1NTS,pr1NTS,':','LineWidth',1.4);
hold on;
plot(re1KATZ,pr1KATZ,':y','LineWidth',1.4);
hold on;
plot(re1NBLPI,pr1NBLPI,':r','LineWidth',1.4);
hold on;
plot(re1Bi,pr1Bi,':g','LineWidth',1.4);
hold on;
plot(re1BPNN,pr1BPNN,':m','LineWidth',1.4);
hold on;
plot(re1Pred,pr1Pred,':c','LineWidth',1.4);
hold on;
plot(re1LRLS,pr1LRLS,':k','LineWidth',1.4);
%%
p1=flip(pr1KATZ);
r1=flip(re1KATZ);
area1(1,1)=p1(1,1)*r1(1,1)/2;
for k=2:10939
    area1(1,k)=[p1(1,k-1)+p1(1,k)]*[r1(1,k)-r1(1,k-1)]/2;
end
aupr1KATZ=sum(area1);
save aupr1KATZ aupr1KATZ;
%%
p2=flip(pr1NBLPI);
r2=flip(re1NBLPI);
area2(1,1)=p2(1,1)*r2(1,1)/2;
for k=2:10939
    area2(1,k)=[p2(1,k-1)+p2(1,k)]*[r2(1,k)-r2(1,k-1)]/2;
end
aupr1NBLPI=sum(area2);
save aupr1NBLPI aupr1NBLPI;
%%
p3=flip(pr1NTS);
r3=flip(re1NTS);
area3(1,1)=p3(1,1)*r3(1,1)/2;
for k=2:10939
    area3(1,k)=[p3(1,k-1)+p3(1,k)]*[r3(1,k)-r3(1,k-1)]/2;
end
aupr1NTS=sum(area3);
save aupr1NTS aupr1NTS;
%%
p4=flip(pr1Bi);
r4=flip(re1Bi);
area4(1,1)=p4(1,1)*r4(1,1)/2;
for k=2:10939
    area4(1,k)=[p4(1,k-1)+p4(1,k)]*[r4(1,k)-r4(1,k-1)]/2;
end
aupr1Bi=sum(area4);
save aupr1Bi aupr1Bi;

%%
p5=flip(pr1BPNN);
r5=flip(re1BPNN);
area5(1,1)=p5(1,1)*r5(1,1)/2;
for k=2:10939
    area5(1,k)=[p5(1,k-1)+p5(1,k)]*[r5(1,k)-r5(1,k-1)]/2;
end
aupr1BPNN=sum(area5);
save aupr1BPNN aupr1BPNN;
%%
p6=flip(pr1Pred);
r6=flip(re1Pred);
area6(1,1)=p6(1,1)*r6(1,1)/2;
for k=2:10939
    area6(1,k)=[p6(1,k-1)+p6(1,k)]*[r6(1,k)-r6(1,k-1)]/2;
end
aupr1Pred=sum(area6);
save aupr1Pred aupr1Pred;
%%
p7=flip(pr1LRLS);
r7=flip(re1LRLS);
area7(1,1)=p7(1,1)*r7(1,1)/2;
for k=2:10939
    area7(1,k)=[p7(1,k-1)+p7(1,k)]*[r7(1,k)-r7(1,k-1)]/2;
end
aupr1LRLS=sum(area7);
save aupr1LRLS aupr1LRLS;


%aupr_MBiRW=trapz(r,p);
%save aupr_MBiRW aupr_MBiRW;
xlabel('Recall ');
ylabel('Precision');
title('PR curve');
legend({'NTSMDA(aupr=0.6947)','KATZHMDA(aupr=0.6178)','NBLPIHMDA(aupr=0.6885)','BiRWMP(aupr=0.6733)','BPNNHMDA(aupr=0.7301)','HMDA-Pred(aupr=0.8754)','LRLSHMDA(aupr=0.6847)'},'FontSize',12,'Location','NorthEast')