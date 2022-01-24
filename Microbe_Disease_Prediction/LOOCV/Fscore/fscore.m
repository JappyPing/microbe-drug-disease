function fscore()
%%
load positionKATZ.mat;
load interaction.mat;
[n,m]=size(interaction);
k=3000;
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
    fp = t - tp;
    fn = m1 - tp;
    tn = n*m - m1 - fn - t;
    precision=tp/(tp+fp);
    recall=tp/(tp+fn);
    FscoreKATZ(1, t)=2*precision*recall/(precision+recall);
end
    save FscoreKATZ FscoreKATZ;
    %%
 load positionNBLPI.mat;
load interaction.mat;
[n,m]=size(interaction);
k=3000;
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
    fp = t - tp;
    fn = m1 - tp;
    tn = n*m - m1 - fn - t;
    precision=tp/(tp+fp);
    recall=tp/(tp+fn);
    FscoreNBLPI(1, t)=2*precision*recall/(precision+recall);
end
    save FscoreNBLPI FscoreNBLPI;
    %%
    load positionBi.mat;
load interaction.mat;
[n,m]=size(interaction);
k=3000;
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
    fp = t - tp;
    fn = m1 - tp;
    tn = n*m - m1 - fn - t;
    precision=tp/(tp+fp);
    recall=tp/(tp+fn);
    FscoreBi(1, t)=2*precision*recall/(precision+recall);
end
    save FscoreBi FscoreBi;
  %%
  load positionBPNN.mat;
load interaction.mat;
[n,m]=size(interaction);
k=3000;
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
    fp = t - tp;
    fn = m1 - tp;
    tn = n*m - m1 - fn - t;
    precision=tp/(tp+fp);
    recall=tp/(tp+fn);
    FscoreBPNN(1, t)=2*precision*recall/(precision+recall);
end
    save FscoreBPNN FscoreBPNN;
    %%
    load positionNTS.mat;
load interaction.mat;
[n,m]=size(interaction);
k=3000;
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
    fp = t - tp;
    fn = m1 - tp;
    tn = n*m - m1 - fn - t;
    precision=tp/(tp+fp);
    recall=tp/(tp+fn);
    FscoreNTS(1, t)=2*precision*recall/(precision+recall);
end
    save FscoreNTS FscoreNTS;
    %%
    load positionPred.mat;
load interaction.mat;
[n,m]=size(interaction);
k=3000;
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
    fp = t - tp;
    fn = m1 - tp;
    tn = n*m - m1 - fn - t;
    precision=tp/(tp+fp);
    recall=tp/(tp+fn);
    FscorePred(1, t)=2*precision*recall/(precision+recall);
end
    save FscorePred FscorePred;
    %%
    load positionLRLS.mat;
load interaction.mat;
[n,m]=size(interaction);
k=3000;
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
    fp = t - tp;
    fn = m1 - tp;
    tn = n*m - m1 - fn - t;
    precision=tp/(tp+fp);
    recall=tp/(tp+fn);
    FscoreLRLS(1, t)=2*precision*recall/(precision+recall);
end
    save FscoreLRLS FscoreLRLS;
end