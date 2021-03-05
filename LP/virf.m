function yf = virf(beta,beta0,h,shock)
%This proc calculates the impulse response function for a typical VAR
global T p vars

yf = zeros(h+p,vars);
fy = zeros(h+p,vars);
beta = beta(2:rows(beta),:);
fy[p+1,.] = shock'*inv(beta0);
yf[p+1,.] = fy[p+1,.];
i = 2;
while i<=h;
    fy[p+i,.] = (vecr(yf[i-1+p:i,.]))'*beta;
    yf[p+i,.] = fy[p+i,.];
    i = i+1;
endo;
yf = yf[p+1:rows(yf),.];
retp(yf);
end