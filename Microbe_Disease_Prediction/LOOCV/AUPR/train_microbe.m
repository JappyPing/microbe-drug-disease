function train_microbe(l,alpha,rho)
%train_microbe(0.1,1.6,6)

%interaction: adjacency matrix for the disease-microbe association network
%interaction(i,j)=1 means microbe j is related to disease i
load interaction;

% Gaussian kernel similarity between microbe: km
load km;
W1=km;
W2=km;

% nd:the number of diseases
% nm:the number of microbe
[nd nm]=size(interaction);
interaction=interaction';

% Sita1 is the bias vector of hidden layer 
Sita1 = rand(1,nm);

% Sita2 is the bias vector of hidden layer 
Sita2 = rand(1,nm);

% Across channel normalization scheme for input signals
for i=1:nm
    [~,N]=size(find(interaction(i,:)));
    if N~=0
     Sign(i,:)=interaction(i,:)/N;
    else
     Sign(i,:)=interaction(i,:);
    end
end



for num=1:1000

%data normalization scheme for weights and biases
W1=normlize(-0.5,0.5,W1);
W2=normlize(-0.5,0.5,W2);
Sita1=normlize(-rho,rho,Sita1);
Sita2=normlize(-rho,rho,Sita2);


%S is the input signals
S = Sign;


%input and output of hidden layer
Hi = zeros(nm,nd);
for i=1:nm
    for j=1:nm
        Hi(i,:) = Hi(i,:) + S(j,:)*W1(j,i);
    end
     Hi(i,:) = Hi(i,:) + Sita1(i);
    
end
Hi=normlize(-alpha,alpha,Hi);
Ho = TanH(Hi);
Ho=normlize(0,1,Ho);


%input and output of output layer
Oi = zeros(nm,nd);
for i=1:nm
    for j=1:nm
        Oi(i,:) = Oi(i,:) + Ho(j,:)*W2(j,i);
    end
     Oi(i,:) = Oi(i,:) + Sita2(i);
    
end

Oi=normlize(-alpha,alpha,Oi);
Oo = TanH(Oi);
Oo=normlize(0,1,Oo);

%calculate errors of output layer
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

%calculate errors of hidden layer
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


%update weights between hidden layer and output layer 
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


%update weights between input layer and hidden layer 
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

%update biases
for i=1:nm
    Sita2(i) =Sita2(i) + l*Err2(i);
    Sita1(i) =Sita1(i) + l*Err1(i);
end


Err(num) = mean(abs(Err2));
     if Err(num)<=0.001||(num>100&&Err(num)/Err(num-1)>1)
    break;

end


end

% plot(Err);
output = Oo;
 save('output.mat','output'); 
end



