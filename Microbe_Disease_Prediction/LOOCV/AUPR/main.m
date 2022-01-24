%KATZHMDAcv(1,1,0.01,2);
%NBLPIHMDAcv(1,1,0.2,5,5,0.75);
%BiRWMPcv(1,1,0.01,2);
%BPNNHMDAcv();
%NTSHMDAcv(1,1,0.7,0.9,0.3,0.8,0.2);
%Predcv();
%LRLSHMDAcv(1,1,0.5);
%%
%NTSHMDA£¨blue£©
load positionNTS.mat;
load interaction;
[n,m] = size(interaction);
k=m*n;
%set threshold
for t=1:k
    [n1,m1]=size(positionNTS);
    tp=0;
    fp=0;
    tn=0;
    fn=0;
    for i=1:m1
        if positionNTS(1,i)<=t
            tp=tp+1;
        end
    end
    fp=m1-tp;
    fn=t-m1;
    tn=n*m-m1-fn;
    pr(t)=tp/(tp+fp);
    re(t)=tp/(tp+fn);
end
pr1NTS=pr(:,450:k);
re1NTS=re(:,450:k);
save pr1NTS pr1NTS
save re1NTS re1NTS
%%
%KATZHMDA (yellow)
load positionKATZ.mat;
load interaction;
[n,m] = size(interaction);
k=m*n;
%set threshold
for t=1:k
    [n1,m1]=size(positionKATZ);
    tp=0;
    fp=0;
    tn=0;
    fn=0;
    for i=1:m1
        if positionKATZ(1,i)<=t
            tp=tp+1;
        end
    end
    fp=m1-tp;
    fn=t-m1;
    tn=n*m-m1-fn;
    pr(t)=tp/(tp+fp);
    re(t)=tp/(tp+fn);
end
pr1KATZ=pr(:,450:k);
re1KATZ=re(:,450:k);
save pr1KATZ pr1KATZ
save re1KATZ re1KATZ
%%
%NBLPIHMDA(red)
load positionNBLPI.mat;
load interaction;
[n,m] = size(interaction);
k=m*n;
%set threshold
for t=1:k
    [n1,m1]=size(positionNBLPI);
    tp=0;
    fp=0;
    tn=0;
    fn=0;
    for i=1:m1
        if positionNBLPI(1,i)<=t
            tp=tp+1;
        end
    end
    fp=m1-tp;
    fn=t-m1;
    tn=n*m-m1-fn;
    pr(t)=tp/(tp+fp);
    re(t)=tp/(tp+fn);
end
pr1NBLPI=pr(:,450:k);
re1NBLPI=re(:,450:k);
save pr1NBLPI pr1NBLPI
save re1NBLPI re1NBLPI
%%
%BiRWMP(green)
load positionBi.mat;
load interaction;
[n,m] = size(interaction);
k=m*n;
%set threshold
for t=1:k
    [n1,m1]=size(positionBi);
    tp=0;
    fp=0;
    tn=0;
    fn=0;
    for i=1:m1
        if positionBi(1,i)<=t
            tp=tp+1;
        end
    end
    fp=m1-tp;
    fn=t-m1;
    tn=n*m-m1-fn;
    pr(t)=tp/(tp+fp);
    re(t)=tp/(tp+fn);
end
pr1Bi=pr(:,450:k);
re1Bi=re(:,450:k);
save pr1Bi pr1Bi
save re1Bi re1Bi
%%
%BPNNHMDA(fenghong)
load positionBPNN.mat;
load interaction;
[n,m] = size(interaction);
k=m*n;
%set threshold
for t=1:k
    [n1,m1]=size(positionBPNN);
    tp=0;
    fp=0;
    tn=0;
    fn=0;
    for i=1:m1
        if positionBPNN(1,i)<=t
            tp=tp+1;
        end
    end
    fp=m1-tp;
    fn=t-m1;
    tn=n*m-m1-fn;
    pr(t)=tp/(tp+fp);
    re(t)=tp/(tp+fn);
end
pr1BPNN=pr(:,450:k);
re1BPNN=re(:,450:k);
save pr1BPNN pr1BPNN
save re1BPNN re1BPNN
%%
%HMDA_Pred(c)
load positionPred.mat;
load interaction;
[n,m] = size(interaction);
k=m*n;
%set threshold
for t=1:k
    [n1,m1]=size(positionPred);
    tp=0;
    fp=0;
    tn=0;
    fn=0;
    for i=1:m1
        if positionPred(1,i)<=t
            tp=tp+1;
        end
    end
    fp=m1-tp;
    fn=t-m1;
    tn=n*m-m1-fn;
    pr(t)=tp/(tp+fp);
    re(t)=tp/(tp+fn);
end
pr1Pred=pr(:,450:k);
re1Pred=re(:,450:k);
save pr1Pred pr1Pred
save re1Pred re1Pred
%%
%LRLSHMDA(k:black)
load positionLRLS.mat;
load interaction;
[n,m] = size(interaction);
k=m*n;
%set threshold
for t=1:k
    [n1,m1]=size(positionLRLS);
    tp=0;
    fp=0;
    tn=0;
    fn=0;
    for i=1:m1
        if positionLRLS(1,i)<=t
            tp=tp+1;
        end
    end
    fp=m1-tp;
    fn=t-m1;
    tn=n*m-m1-fn;
    pr(t)=tp/(tp+fp);
    re(t)=tp/(tp+fn);
end
pr1LRLS=pr(:,450:k);
re1LRLS=re(:,450:k);
save pr1LRLS pr1LRLS
save re1LRLS re1LRLS
%%
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

xlabel('Recall ');
ylabel('Precision');
title('LOOCV-PR curve');
legend({'NTSMDA(aupr=0.6947)','KATZHMDA(aupr=0.6178)','NBLPIHMDA(aupr=0.6885)','BiRWMP(aupr=0.6733)','BPNNHMDA(aupr=0.7301)','HMDA-Pred(aupr=0.8416)','LRLSHMDA(aupr=0.6384)'},'FontSize',16,'Location','NorthEast')