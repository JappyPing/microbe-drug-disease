for cv=1:10
KATZHMDA10cv(1,1,0.01,2);
overallauc(cv)=positiontooverallauc();
end
save overallauc overallauc
a=mean(overallauc)
b=std(overallauc)