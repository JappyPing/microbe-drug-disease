function output = normlize(a,b,x)
[nd,nm] = size(x);
MinValue = min(min(x(:,:)));
MaxValue = max(max(x(:,:)));
if nd>1
    for i=1:nd
        for j=1:nm
            output(i,j) = ((b-a)*(x(i,j)-MinValue)/(MaxValue-MinValue))+a;
        end
    end
    else
        for i=1:nm
            if MinValue==0&&MaxValue==0
                output(i) = 0;
            else
        output(i) = ((b-a)*(x(i)-MinValue)/(MaxValue-MinValue))+a;
             end
        end

end

