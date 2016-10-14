function [t]=Sample_changepoint_refine(x,prev_x,L0,term,N,t_start,BL,pi,m)
t=ones(1,N);
t(1:L0)=0;
t(~isnan(prev_x.t))=prev_x.t(~isnan(prev_x.t));
nodes=1:N;

for n1=nodes(isnan(prev_x.t) & t>0)
    t_n1=Sample_changepoint(x.t-t_start+1,x.z,n1,x.a,BL,term(n1,:),pi,m);
    t(n1)=t_n1+t_start-1;
end

