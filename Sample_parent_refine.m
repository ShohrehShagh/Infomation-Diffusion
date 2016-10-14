function [z]=Sample_parent_refine(x,prev_x,pi,N,L0,t_start_batch,BL)
z=ones(1,N);
nodes=1:N;

z(1:L0)=0;

inds1=isnan(prev_x.t) & isnan(x.t) & z>0;
inds2=~isnan(prev_x.t) & ~isnan(x.t) & z>0;
inds3=isnan(prev_x.t) & ~isnan(x.t) & z>0;

z(inds1)=NaN;
z(inds2)=prev_x.z(inds2);

for n1=nodes(inds3)
    z(n1)=Sample_parent(x.t,x.a,x.z,pi,n1,t_start_batch,BL);
end

