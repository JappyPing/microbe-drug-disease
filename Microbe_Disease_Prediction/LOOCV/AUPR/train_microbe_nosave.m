function output = train_microbe_nosave(l,num,W1,W2,Sita1,Sita2,interaction,threshold,range,uplimit)
% load interaction;
[nd nm]=size(interaction);
interaction=interaction';

for i=1:nm
    [~,N]=size(find(interaction(i,:)));
    if N~=0
     Sign(i,:)=interaction(i,:)/N;
    else
     Sign(i,:)=interaction(i,:);
    end
end
% Sign=interaction;

for num=1:num

W1=normlize(-0.5,0.5,W1);
W2=normlize(-0.5,0.5,W2);
Sita1=normlize(-range,range,Sita1);
Sita2=normlize(-range,range,Sita2);
%        MinValue = min(min(W1(:,:)));
% MaxValue = max(max(W2(:,:))); 

%输入层赋值
S = Sign;


%隐藏层输入输出
Hi = zeros(nm,nd);
for i=1:nm
    for j=1:nm
        Hi(i,:) = Hi(i,:) + S(j,:)*W1(j,i);
    end
     Hi(i,:) = Hi(i,:) + Sita1(i);
    
end
Hi=normlize(-uplimit,uplimit,Hi);
Ho = TanH(Hi);
Ho=normlize(0,1,Ho);
% for i=1:39
%     Ho(:,i)=normlize(0,1,Ho(:,i));
% end

%输出层输入输出
Oi = zeros(nm,nd);
for i=1:nm
    for j=1:nm
        Oi(i,:) = Oi(i,:) + Ho(j,:)*W2(j,i);
    end
     Oi(i,:) = Oi(i,:) + Sita2(i);
    
end

Oi=normlize(-uplimit,uplimit,Oi);
Oo = TanH(Oi);
% for i=1:39
%     Oo(:,i)=normlize(0,1,Oo(:,i));
% end
Oo=normlize(0,1,Oo);

%计算误差，输出层
Err2 = zeros(1,nm);
Error = zeros(1,nm);
for i=1:nm
    aim = find(S(i,:));
    N = length(aim);
    if N == 0
        continue;
    end
    for j=1:N
    Err2(i) =Err2(i) + (1-Oo(i,aim(j))^2)*(1-Oo(i,aim(j)));
    Error(i) = Error(i) + (1-Oo(i,aim(j)));
    end
    Err2(i) = Err2(i)/N;
    Error(i) = Error(i)/N;
    Error(i) = 0.5*Error(i)^2;
end

%计算误差，隐藏层
Err1 = zeros(1,nm);
ErrW = zeros(1,nm);
for i=1:nm
    for j=1:nm
    ErrW(i) = ErrW(i) + W2(i,j)*Err2(j);
    end
end

for i=1:nm
    aim = find(S(i,:));
    N = length(aim);
    if N == 0
        continue;
    end
    for j=1:N
    Err1(i) = Err1(i) + (1-Ho(i,aim(j))^2);
    end
     Err1(i) = Err1(i)/N;
     Err1(i) = Err1(i)*ErrW(i);
end

% % 权重2更新
% for i=1:nm
%     meanHo=mean(Ho(i,:));
%     for j=1:nm
%         W2(i,j) = W2(i,j) + l * Err2(j) *  meanHo;
%     end
% end
% %权重1更新
% for i=1:nm
%     meanS=mean(S(i,:));
%     for j=1:nm
%         W1(i,j) = W1(i,j) + l * Err1(j) * meanS;
%     end
% end

%权重2更新
for i=1:nm
    MeanHo=0;
    aim = find(S(i,:));
    N = length(aim);
    if N == 0
        continue;
    end
        for j=1:N
    MeanHo = MeanHo + Ho(i,aim(j));
    end
%     MeanHo = MeanHo/N;
    
    for j=1:nm
        W2(i,j) = W2(i,j) + l * Err2(j) * MeanHo;
    end
end


%权重1更新

for i=1:nm
    MeanS=0;
    aim = find(S(i,:));
    N = length(aim);
    if N == 0
        continue;
    end
    for j=1:N
    MeanS = MeanS + S(i,aim(j));
    end
%     MeanS = MeanS/N;
    for j=1:nm
        W1(i,j) = W1(i,j) + l * Err1(j) * MeanS;
    end
end

%偏向更新
for i=1:nm
    Sita2(i) =Sita2(i) + l*Err2(i);
    Sita1(i) =Sita1(i) + l*Err1(i);
end


%Err = sum(abs(Err2));
Err(num) = mean(abs(Err2));
% Err(num) = mean(Error);
     if Err(num)<=threshold||(num>1&&Err(num)/Err(num-1)>1)
    break;

end


end
clear area;
% plot(Err);
output = Oo';
end



