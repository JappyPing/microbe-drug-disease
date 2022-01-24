function output = TanH(x)
% output = 2./(1+exp(-2*x))-1;
[a b]=size(x);
for i=1:a
    for j=1:b
output(i,j) = (exp(x(i,j))-exp(-x(i,j)))/(exp(x(i,j))+exp(-x(i,j)));
    end
end
end

