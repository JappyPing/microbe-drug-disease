%%
load positionKATZ.mat
load interaction.mat
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
prKATZ=pr(:,450:k);
reKATZ=re(:,450:k);
save prKATZ prKATZ
save reKATZ reKATZ
%%
load positionNBLPI.mat
load interaction.mat
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
prNBLPI=pr(:,450:k);
reNBLPI=re(:,450:k);
save prNBLPI prNBLPI
save reNBLPI reNBLPI
%%
load positionBi.mat
load interaction.mat
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
prBi=pr(:,450:k);
reBi=re(:,450:k);
save prBi prBi
save reBi reBi
%%
load positionBPNN.mat
load interaction.mat
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
prBPNN=pr(:,450:k);
reBPNN=re(:,450:k);
save prBPNN prBPNN
save reBPNN reBPNN
%%
load positionNTS.mat
load interaction.mat
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
prNTS=pr(:,450:k);
reNTS=re(:,450:k);
save prNTS prNTS
save reNTS reNTS
%%
load positionPred.mat
load interaction.mat
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
prPred=pr(:,450:k);
rePred=re(:,450:k);
save prPred prPred
save rePred rePred
%%
load positionLRLS.mat
load interaction.mat
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
prLRLS=pr(:,450:k);
reLRLS=re(:,450:k);
save prLRLS prLRLS
save reLRLS reLRLS