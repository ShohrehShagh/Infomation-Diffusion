function [term]=term_calculator(data,mu1,mu2,sigma1,sigma2)
T=size(data,2);
term=zeros(size(data));
for i=1:T+1
    term(:,i)=prod(normpdf(data(:,1:i-1),mu1,sigma1),2).*prod(normpdf(data(:,i:T),mu2,sigma2),2);
end
%term1=cumsum(((data-mu1).^2)/(2*sigma1^2),2);
%term2=diag(sum((data-mu2).^2,2)/(2*sigma2^2))*ones(size(data))-cumsum(((data-mu2).^2),2)/(2*sigma2^2);

% BT=size(data,2);
% term1=zeros(N,BT);
% term2=zeros(N,BT);
% term1(1:N,1)=((data(1:N,1)-mu1).^2)/(2*sigma1^2);
% term2(1:N,1)=sum((data(1:N,2:BT)-mu2).^2,2)/(2*sigma2^2);
% for t=2:BT
%     term1(1:N,t)=term1(1:N,t-1)+((data(1:N,t)-mu1).^2)/(2*sigma1^2);
%     term2(1:N,t)=term2(1:N,t-1)-((data(1:N,t)-mu2).^2)/(2*sigma2^2);
% end