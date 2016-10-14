function [alphas]=Sample_alpha_online(t,prev_t,prev_alphas,a_mat,b_mat,pi,N)
alphas=zeros(N,N);
nodes=1:N;

inds1=(isnan(prev_t) & isnan(t))|(isnan(prev_t) & ~isnan(t));
inds2=~isnan(prev_t) & ~isnan(t);

for n1=nodes(inds1)
    inds1_prime=pi(:,n1)==1;
    alphas(n1,inds1_prime)=gamrnd(a_mat(n1,inds1_prime),b_mat(n1,inds1_prime));
    alphas(inds1_prime,n1)=alphas(n1,inds1_prime);
end
for n1=nodes(inds2)
    inds2_prime=pi(:,n1)==1;
    alphas(n1,inds2_prime)=prev_alphas(n1,inds2_prime);
    alphas(inds2_prime,n1)=alphas(n1,inds2_prime);
end
