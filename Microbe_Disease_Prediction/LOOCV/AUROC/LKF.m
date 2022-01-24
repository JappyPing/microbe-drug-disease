function [ sim ] = LKF( w,sim1,sim2 )
sim = w*sim1 + (1-w)*sim2;
end

