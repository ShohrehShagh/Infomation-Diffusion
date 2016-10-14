function [alphas]=Sample_alpha_refine(x,prev_x,a_mat,b_mat,pi,N,BL,pdf_val)
nodes=1:N;

inds1=(isnan(prev_x.t) & isnan(x.t))|(isnan(prev_x.t) & ~isnan(x.t));
inds2=~isnan(prev_x.t) & ~isnan(x.t);

for n1=nodes(inds1)
    alphas=Sample_alpha(x.z,n1,x.t,x.a,a_mat(n1,:),b_mat(n1,:),pi,N,BL,pdf_val);
end
for n1=nodes(inds2)
    inds2_prime=pi(:,n1)==1;
    alphas(n1,inds2_prime)=prev_x.a(n1,inds2_prime);
    alphas(inds2_prime,n1)=alphas(n1,inds2_prime);
end
