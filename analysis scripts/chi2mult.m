function [p,chi2,chi2pair,pTadj,tbl] = chi2mult(X,n,g1,g2)
%% chi squared test for multiple groups
% test whether the proportions are indipendent

if size(n,1)>1
    n = reshape(n,1,length(n));
    X = X';
end

if nargin<3
   g1 = num2cell(1:size(X,1));
end

if nargin<4
   g2 = num2cell(1:size(X,2));
end

if sum(X(:,1))<n(1)
    X(end+1,:) = n-sum(X,1);
    g1(end+1) = {'Other'};
end
g1 = reshape(g1,length(g1),1);
%% expected frequencies
m = sum(X,2);
N = sum(n);
Xexp = m*n./N;

%% Calculated chi^2
chi2 = 0;
for i = 1:length(X(:))
    chi2 = chi2 + (X(i)-Xexp(i))^2/Xexp(i);
end

%% degress of freedom
df = (size(X,1)-1)*(size(X,2)-1); 
p = chi2cdf(chi2,df,'upper');

%%
tbl = {'Source' 'df' 'Chi^2' 'p';
    'Groups' df chi2 p};
%% Test statistics

% % Tstd = zeros(size(X));
Tadj = zeros(size(X));
pTadj = zeros(size(X));
sigTadj = zeros(size(X));
df = 1;
for i = 1:size(X,1)
    for j = 1:size(X,2)
        Tstd(i,j) = (X(i,j)-Xexp(i,j))/sqrt(Xexp(i,j));
        Tadj(i,j) = (X(i,j)-Xexp(i,j))/sqrt(Xexp(i,j)*(1-n(j)/N)*(1-m(i)/N));
        pTadj(i,j) = chi2cdf(Tadj(i,j)^2,df,'upper');        
    end
end

sigbound = zeros(size(X));
% sigbound(:) = .05/(.5*N*(N-1)); % compare all entries
sigbound(:) = .05/(.5*length(n)*(length(n)-1)); % compare one line
% sigbound(:) = .05/(length(n)-1); %compare first value to the others
% [~,b] = sort(pTadj,2);
% sigbound = sigbound.*b;
for i = 1:size(X,1)
    for j = 1:size(X,2)
        sigTadj(i,j) = pTadj(i,j)<sigbound(i,j);
    end
end
sigTadj = cat(3,Tadj,pTadj,sigTadj);

pTadj = num2cell(pTadj);
pTadj = cat(2,g1,pTadj);
g2 = cat(2,{'Ind p'},g2);
pTadj = cat(1,g2,pTadj);
%% Pairwise Bonferroni -  Method

df = 1;
k = 1;
testg = 1;
for i = 1:length(n)-1
    for j = i+1:length(n)
        chi2pair(k,1) = i;
        chi2pair(k,2) = j;
        chi2pair(k,3) = (X(testg,i)-Xexp(testg,i))^2/Xexp(testg,i)+(X(end,i)-Xexp(end,i))^2/Xexp(end,i)+(X(testg,j)-Xexp(testg,j))^2/Xexp(testg,j)+(X(end,j)-Xexp(end,j))^2/Xexp(end,j);
        chi2pair(k,4) = chi2cdf(chi2pair(k,3),df,'upper');
        k = k+1;
    end
end
sigbound = zeros(size(chi2pair,1),1);
sigbound(:) = .05/(.5*length(n)*(length(n)-1)); % compare one line
% sigbound(:) = .05/(length(n)-1); %compare first value to the others
% [~,b] = sort(chi2pair(:,4),1);
% sigbound = b.*sigbound;
k = 1;
for i = 1:length(n)-1
    for j = i+1:length(n)
        chi2pair(k,5) = chi2pair(k,4)<sigbound(k);
        k = k+1;
    end
end
     

%% P Value fisher exact test

% for i = 1:size(X,1)
%     for j = 1:size(X,2)
%         P(i,j) = factorial(m(i))*factorial(N-m(i))*factorial(n(j))*factorial(N-n(j))/((factorial(X(i,j))*factorial(n(j)-X(i,j))*factorial(m(i)-X(i,j))*factorial(N-m(i)-n(j)+X(i,j)))*factorial(N));
%     end
% end



end
