function []=data_generation(N,mu1,mu2,sigma1,sigma2,a_mat,b_mat,pi,L0,data_num,graph,r)
% strengths f(\alpha)
alpha_R=zeros(N,N);
alpha_R(pi==1)=gamrnd(a_mat(pi==1),b_mat(pi==1));
alpha_R=tril(alpha_R,-1)'+alpha_R;
% parents f(z|\alpha)
z_R=zeros(1,N);
for n=1:N
    candidate_parents=find(pi(n,:));
    if(~isempty(candidate_parents))
        FCP_Z=alpha_R(n,candidate_parents);
        z_R(n)=candidate_parents(find(rand <= cumsum(FCP_Z)/sum(FCP_Z), 1));
    end
end
z_R=graph(:,1)';
% changepoints f(t_c|z,\alpha) & data f(d|t_c)
tc_R=zeros(1,N);
for i=1:N
    p=z_R(i);
    if(p==0)
        tc_R(i)=0;
    else
        p_alpha=1-exp(-alpha_R(i,p));
        tc_R(i)=tc_R(p)+geornd(p_alpha)+1;
    end
end

T=max(tc_R)+5;
all_data=zeros(N,T);
for i=L0+1:N
    Changepoint_index=tc_R(i);
    all_data(i,1:Changepoint_index-1)=normrnd(mu1,sigma1,[1,Changepoint_index-1]);
    all_data(i,Changepoint_index:T)=normrnd(mu2,sigma2,[1,T-Changepoint_index+1]);
end
Real.alpha=alpha_R;
Real.t=tc_R;
Real.z=z_R;
Real.T=T;
Real.d=all_data;
filename=sprintf('input_datas/Data_%d_%d',data_num,r);
save(filename,'Real');
Real