%KATZHMDA10cv(1,1,0.01,2);
%NBLPIHMDA10cv(1,1,0.2,5,5,0.75);
%BiRWMP10cv(1,1,0.01,2);
%BPNNHMDA10cv();
%NTSHMDA10cv(1,1,0.7,0.9,0.3,0.8,0.2);
%Pred10cv();
%LRLSHMDA10cv(1,1,0.5);
%NTSHMDA£¨blue£©
load positionNTS.mat;
load interaction;
[n,m]=size(interaction);
sID=textread('knowndiseasemicrobeinteraction.txt');
[pp,qq]=size(sID);
for i=1:pp
if positionNTS(i)>m*n-pp+1
positionNTS(i)=m*n-pp+1;
end
end
for k=1:m*n-pp+1
    tp=0;
   for t=1:pp
        if positionNTS(1,t)<=k
           tp=tp+1;
        end
    end
    tpr(1,k)=tp/pp; 
    fp=k*pp-tp;
    fpr(1,k)=fp/(pp*(m*n-pp)); 
end
%KATZHMDA (yellow)
load positionKATZ.mat;
load interaction;
[n1,m1]=size(interaction);
sID1=textread('knowndiseasemicrobeinteraction.txt');
[pp1,qq1]=size(sID1);


for i=1:pp1
if positionKATZ(i)>m1*n1-pp1+1
positionKATZ(i)=m1*n1-pp1+1;
end
end
for k=1:m1*n1-pp1+1
    tp1=0;
    for t=1:pp1
        if positionKATZ(1,t)<=k
            tp1=tp1+1;
        end
    end
    tpr1(1,k)=tp1/pp1; 
    fp1=k*pp1-tp1;
     fpr1(1,k)=fp1/(pp1*(m1*n1-pp1)); 
end
%NBLPIHMDA(red)
load positionNBLPI.mat;
load interaction;
[n2,m2]=size(interaction);
sID2=textread('knowndiseasemicrobeinteraction.txt');
[pp2,qq2]=size(sID2);


for i=1:pp2
if positionNBLPI(i)>m2*n2-pp2+1
positionNBLPI(i)=m2*n2-pp2+1;
end
end
for k=1:m2*n2-pp2+1
    tp2=0;
    for t=1:pp2
        if positionNBLPI(1,t)<=k
            tp2=tp2+1;
        end
    end
    tpr2(1,k)=tp2/pp2; 
    fp2=k*pp2-tp2;
     fpr2(1,k)=fp2/(pp2*(m2*n2-pp2)); 
end
%BiRWMP(green)
load positionBi.mat;
load interaction;
[n4,m4]=size(interaction);
sID4=textread('knowndiseasemicrobeinteraction.txt');
[pp4,qq4]=size(sID4);


for i=1:pp4
if positionBi(i)>m4*n4-pp4+1
positionBi(i)=m4*n4-pp4+1;
end
end
for k=1:m4*n4-pp4+1
    tp4=0;
    for t=1:pp4
        if positionBi(1,t)<=k
            tp4=tp4+1;
        end
    end
    tpr4(1,k)=tp4/pp4; 
    fp4=k*pp4-tp4;
     fpr4(1,k)=fp4/(pp4*(m4*n4-pp4)); 
end
%BPNNHMDA(fenghong)
load positionBPNN.mat;
load interaction;
[n5,m5]=size(interaction);
sID5=textread('knowndiseasemicrobeinteraction.txt');
[pp5,qq5]=size(sID5);


for i=1:pp5
if positionBPNN(i)>m5*n5-pp5+1
positionBPNN(i)=m5*n5-pp5+1;
end
end
for k=1:m5*n5-pp5+1
    tp5=0;
    for t=1:pp5
        if positionBPNN(1,t)<=k
            tp5=tp5+1;
        end
    end
    tpr5(1,k)=tp5/pp5; 
    fp5=k*pp5-tp5;
     fpr5(1,k)=fp5/(pp5*(m5*n5-pp5)); 
end
%HMDA_Pred(c)
load positionPred.mat;
load interaction;
[n6,m6]=size(interaction);
sID6=textread('knowndiseasemicrobeinteraction.txt');
[pp6,qq6]=size(sID6);


for i=1:pp6
if positionPred(i)>m6*n6-pp6+1
positionPred(i)=m6*n6-pp6+1;
end
end
for k=1:m6*n6-pp6+1
    tp6=0;
    for t=1:pp6
        if positionPred(1,t)<=k
            tp6=tp6+1;
        end
    end
    tpr6(1,k)=tp6/pp6; 
    fp6=k*pp6-tp6;
     fpr6(1,k)=fp6/(pp6*(m6*n6-pp6)); 
end
%LRLSHMDA(k:black)
load positionLRLS.mat;
load interaction;
[n7,m7]=size(interaction);
sID7=textread('knowndiseasemicrobeinteraction.txt');
[pp7,qq7]=size(sID7);


for i=1:pp7
if positionLRLS(i)>m7*n7-pp7+1
positionLRLS(i)=m7*n7-pp7+1;
end
end
for k=1:m7*n7-pp7+1
    tp7=0;
    for t=1:pp7
        if positionLRLS(1,t)<=k
            tp7=tp7+1;
        end
   end
    tpr7(1,k)=tp7/pp7; 
    fp7=k*pp7-tp7;
     fpr7(1,k)=fp7/(pp7*(m7*n7-pp7)); 
end
plot(fpr,tpr,'LineWidth',1.4);
hold on;
plot(fpr1,tpr1,'y','LineWidth',1.4);
hold on;
plot(fpr2,tpr2,'r','LineWidth',1.4);
hold on;
plot(fpr4,tpr4,'g','LineWidth',1.4);
hold on;
plot(fpr5,tpr5,'m','LineWidth',1.4);
hold on;
plot(fpr6,tpr6,'c','LineWidth',1.4);
hold on;
plot(fpr7,tpr7,'k','LineWidth',1.4);
xlabel('False Positive Rate','FontSize',16);
ylabel('True Positive Rate','FontSize',16);
legend({'NTSMDA(auc=0.8882)','KATZHMDA(auc=0.8354)','NBLPIHMDA(auc=0.9000)','BiRWMP(auc=0.8601)','BPNNHMDA(auc=0.9188)','HMDA-Pred(auc=0.8841)','LRLSHMDA(auc=0.8873)'},'FontSize',16,'Location','SouthEast');
title('10-fold-ROC curve')

