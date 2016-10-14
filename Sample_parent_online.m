function [z]=Sample_parent_online(t,prev_t,alphas,prev_z,pi,N,L0)
z=ones(1,N);
nodes=1:N;

z(1:L0)=0;

inds1=isnan(prev_t) & isnan(t) & z>0;
inds2=~isnan(prev_t) & ~isnan(t) & z>0;
inds3=isnan(prev_t) & ~isnan(t) & z>0;

z(inds1)=NaN;
z(inds2)=prev_z(inds2);

for n1=nodes(inds3)
    v=zeros(N,1);
    v(~isnan(t).* pi(n1,:)>0)=alphas(n1,~isnan(t).* pi(n1,:)>0)/sum(alphas(n1,:).* pi(n1,:));%%
    v(N+1)=1-sum(v);
    index=find(rand <= cumsum(v)/sum(v), 1);
    if(index>N)
        z(n1)=NaN;
    else
        z(n1)=index;
    end
end

