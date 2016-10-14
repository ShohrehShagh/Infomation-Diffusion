function [t]=Sample_changepoint_online(prev_t,BL,p_b,t_ML,L0,t_start,N)
t=ones(1,N);
t(1:L0)=0;
t(t_ML==0 | t_ML==BL+1)=NaN;
t(~isnan(prev_t))=prev_t(~isnan(prev_t));
nodes=1:N;

for n1=nodes(t_ML~=0 & t_ML~=BL+1 & isnan(prev_t) & t>0)
    v=abs((1:BL)-t_ML(n1));
    FCP_C=0.5*p_b*(1-p_b).^v;
    FCP_C(t_ML(n1))=p_b;
    FCP_C(BL+1)=1-sum(FCP_C);
    index=find(rand <= cumsum(FCP_C)/sum(FCP_C), 1);
    if (index>BL)
        t(n1)=NaN;
    else
        t(n1)=index+t_start-1;
    end
end

