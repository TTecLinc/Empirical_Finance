% Peilin Yang 3/5/2021

clear
global T p vars

vars = 6;		% Variables in the VAR
plg = 1;		% LAG LENGTH IN THE VAR. if plg = 0, lag length is set to maxp. if 1, lag length is chosen automatically 
maxp = 14;      % Max Lag for Info Criteria 
infoc = 0;      % 0 for AICc, -1 for AIC, 1 for SIC @
h = 24;			% Length of the Impulse Response @
shock = [0 0 0 0.5 0 0];	% Structural Shock for IRF Computation @
order = [1,2,3,4,5,6];	% Cholesky ordering for structural identification @
sig = 1;		% Standard error bands for IRF: 1 if 1 SD or 66% confidence level, 1.96 for 95% @
nw = 1;        		% Newey-West Lag length is nw+p+h @
z0=readtable('evnew.csv');
z = table2array(readtable('evnew.csv'));
T = size(z,1);
p=14;
%% ------------------------------------------------------------------------------
% [e beta beta0] = varols(z,p);
%	Generate the Dependent Variable Vector

y = z((p+1):T,:);

%	Generate the lags
x = ones(T-p,1);

i = 1;
if i > p;
    display('i>p');
else
	while i<=p;
		x = [x z((p+1-i):(T-i),:)];
		i = i+1;
	end
end


%	OLS on the VAR

beta = (x'*x)\(x'*y);
% A\=inv(A)

%	Get Structural residuals from reduced form

u = y - x*beta;
vcv = u'*u;

% Cholesky Decomposition
% eye: Identity matrix
%A = spdiags(B,d,A) replaces the diagonals specified by d with the columns of B. The output is sparse.
%A(logical(eye(size(A))))=B
%m_0=eye(vars);

a0 = chol(vcv)\diag(diag(chol(vcv)));
e = u*a0;

beta0=a0;
%% ------------------------------------------------------------------------------
%This proc calculates the impulse response function for a typical VAR

yf = zeros(h+p,vars);
fy = zeros(h+p,vars);
beta = beta(2:size(beta,1),:);
fy(p+1,:) = shock*inv(beta0);

yf(p+1,:) = fy(p+1,:);
i = 2;
while i<=h;

    m_0=flipud(yf(i:i-1+p,:));
    m_0=m_0';
    fy(p+i,:) =m_0(:)'*beta;
    yf(p+i,:) = fy(p+i,:);
    i = i+1;
end
yf = yf(p+1:size(yf,1),:);
var_f=yf;
% plot(yf(:,1));
% plot(yf(:,2));
% plot(yf(:,3));
%% ---------------------------------------------------------------------------------
%*******************************************************************
%Compute local-projection IRFs
%*******************************************************************/

w = shock/beta0;
a00 = (inv(beta0)').^2;

%	Generate the Dependent Variable Vector	***/

y = z((p+1):T,:);

%	Generate the lags
x = ones(T-p,1);

i = 1;
if i > p;
    display('i>p');
else
	while i<=p;
		x = [x z((p+1-i):(T-i),:)];
		i = i+1;
	end
end

yf = zeros(h+1,vars);
pf = zeros(h+1,vars);
nf = zeros(h+1,vars);
vf = zeros(h+1,vars);

yf(1,:) = w;
pf(1,:) = w;
nf(1,:) = w;
 

j = 1;
if j > h;
    display('j>h');
else
	while j<=h;
        
        yy = y(j:size(y,1),:);
        xx = x(1:size(x,1)-j+1,:);
        
        gamx = (xx'*xx)\(xx'*yy);
        u = yy - xx*gamx;
       
        if j==1;
            de = diag(beta0'*u'*u*beta0)';
        end;
        
        du = u'*u;
        [m0 p0]=max(shock);
        vf(j,:) = (a00(p0,:)).*de/du(p0,p0);
        vf(j,:) = vf(j,:)./(sum(vf(j,:)));
        
        i = 1;
        % Newey-West Standard Errors @
        l = 1;
        %var is (u-mean(u))'*(u-mean(u)');
        var=cov(u); %is different for freedom
%         var=zeros(vars,vars);
%         for i=1:size(u,2) 
%             for j=1:size(u,2) 
%                 var(i,j)=sum((u(:,i)-mean(u(:,i))).*(u(:,j)-mean(u(:,j))))/(size(u,1)-1);
%             end 
%         end
        
        if nw < 1;
            disp('nw');
        else
            while l<=nw;
                %gammak is (u(1+l:size(u,1),:)-mean(u(1+l:size(u,1),:))')'*(u(1:size(u,1)-l,:)-mean(u(1:size(u,1)-l,:))');
                gammak=cov(u(1+l:size(u,1),:));
                var = var + (1-(l/(nw+1)))*(gammak+gammak');
                l = l+1;
            end;
        end;
        
        while i<=vars;
            sd = (var(i,i)/(size(xx,1)-(vars*p))).*inv(xx'*xx);

            sd = w*sd(2:vars+1,2:vars+1)*w';
            sd = sd^0.5;
            pf(j+1,i) = sig*sd;
            nf(j+1,i) = -sig*sd;
            i = i+1;
        end

        a00 = a00 + (gamx(2:vars+1,:)'/(beta0)').^2;


        yf(j+1,:) = w*gamx(2:vars+1,:);
        pf(j+1,:) = pf(j+1,:)+yf(j+1,:);
        nf(j+1,:) = nf(j+1,:)+yf(j+1,:);
        j = j+1;
    end;
end

zlpf=yf(1:h,:);
pf_lp=pf(1:h,:);
nf_lp=nf(1:h,:);
vldf=vf(1:h,:);

%% Cubic Projection
%*******************************************************************
%Compute Cubic-projection IRFs
%*******************************************************************/

w = (shock/beta0);

%	Generate the Dependent Variable Vector	***/
y=bsxfun(@minus,z((p+1):T,:),mean(z((p+1):T,:)));

%	Generate the lags			***/
x = ones(T-p,1);

x = [x bsxfun(@minus,z(p:(T-1),:),mean(z(p:(T-1),:))) bsxfun(@minus,z(p:(T-1),:),mean(z(p:(T-1),:)).^2) bsxfun(@minus,z(p:(T-1),:),mean(z(p:(T-1),:)).^3)];

i = 2;
if i > p;
    display('i>p');
else
	while i<=p;
		x = [x bsxfun(@minus,z((p+1-i):(T-i),:),mean(z((p+1-i):(T-i),:)))];
		i = i+1;
    end
end

yf = zeros(h+1,vars);
pf = zeros(h+1,vars);
nf = zeros(h+1,vars);

yf(1,:) = w;
pf(1,:) = w;
nf(1,:) = w;



j = 1;
if j > h;
    display('j>k');
else
	while j<=h;
        yy = y(j:size(y,1),:);
        xx = x(1:size(x,1)-j+1,:);
%         save_xx=xx'*xx;
%         u=chol(save_xx);
%         tu=invutri(u);
%         local_B=tu*tu';
	    gamx = (xx'*xx)\(xx'*yy);
        u = yy - xx*gamx;
        
        i = 1;
        %     @ Newey-West Standard Errors @
        l = 1;
        var = cov(u);
        
        if nw < 1;
            disp('nw');
        else
            while l<=nw;
                %gammak is (u(1+l:size(u,1),:)-mean(u(1+l:size(u,1),:))')'*(u(1:size(u,1)-l,:)-mean(u(1:size(u,1)-l,:))');
                gammak=cov(u(1+l:size(u,1),:));
                var = var + (1-(l/(nw+1)))*(gammak+gammak');
                l = l+1;
            end;
        end;
    
    	while i<=vars;
    		sd = (var(i,i)/(size(xx,1)-(vars*p))).*inv(xx'*xx);
    		
    		sd = [w (w.^2) (w.^3)]*sd(2:3*vars+1,2:3*vars+1)*[w (w.^2) (w.^3)]';
    		sd = sd^0.5;
    		pf(j+1,i) = sig*sd;
    		nf(j+1,i) = -sig*sd;
    		i = i+1;
        end
        yf(j+1,:) = w*gamx(2:vars+1,:) + (w.^2)*(gamx((vars+2):(2*vars+1),:))+ (w.^3)*(gamx((2*vars + 2):(3*vars + 1),:));
        pf(j+1,:) = pf(j+1,:)+yf(j+1,:);
        nf(j+1,:) = nf(j+1,:)+yf(j+1,:);
    	j = j+1;
    end
end
zcpf=yf(1:h,:);
pf1=pf(1:h,:);
nf1= nf(1:h,:);
% 
% %% -------------------------------------------
% % Graph-
% 
figure(1);
plot(var_f(:,1));
hold on;
plot(zlpf(:,1));
hold on;
%plot(zcpf(:,1));
legend('VAR','Local Projection');

figure(2);
plot(var_f(:,2));
hold on;
plot(zlpf(:,2));
legend('VAR','Local Projection');


