clc
clear all
Number_of_samples=10^5;
Number_of_runs=10;%00;
Burn_in=floor(Number_of_samples/100);
N_thin=10;
%% initialization
mu1=10;
mu2_vec=[100,100,11,11];
sigma1=1;
sigma2=1;
a1=1;
b1=0.5;
a2_vec=[40,2,40,2];
b2=0.5;
N=20;
L0=1;
p_b=0.5;
pi=zeros(N,N);
for i=2:N
    pi(i,1:i-1)=rand(1,i-1)<0.5;
    if(sum(pi(i,1:i-1))==0)
        pi(i,randsample(i-1,1))=1;
    end
end
%% implied links
links=zeros(N,N);
graph=zeros(N,2);
for i=1:N
    p_vector=find(pi(i,:)>0);
    P=ones(1,length(p_vector));
    if (isempty(P))
        graph(i,1)=0;
        graph(i,2)=i;
    else
        graph(i,1)=p_vector(find(rand <= cumsum(P)/sum(P), 1));
        graph(i,2)=i;
    end
end
for i=1:size(graph,1)
    if (graph(i,1)~=0)
        links(graph(i,1),graph(i,2))=-1;
        links(graph(i,2),graph(i,1))=1;
    end
end
save(sprintf('outputs/Init_Value%d.mat',N),'-v7.3')

N_thin=100;
Burn_in=floor(Number_of_samples/100);
for r=1:Number_of_runs
    %% Data Generation
    b_mat=b1*pi;
    b_mat=tril(b_mat,-1)'+b_mat;
    b_mat(links~=0)=b2;
    LD=length(mu2_vec);
    for mu_index=1:LD
        mu_index
        %% alpha required value  
        candidate_alphas1=linspace(0.01,5,100);
        candidate_alphas2=linspace(10,30,100);
        pdf_val=zeros(2,length(candidate_alphas1));
        pdf_val(1,:)=gampdf(candidate_alphas1,a1,b1);
        pdf_val(2,:)=gampdf(candidate_alphas2,a2_vec(mu_index),b2);
        %%
        mu2=mu2_vec(mu_index);
        a2=a2_vec(mu_index);
        a_mat=a1*pi;
        a_mat=tril(a_mat,-1)'+a_mat;
        a_mat(links~=0)=a2;
        data_generation(N,mu1,mu2,sigma1,sigma2,a_mat,b_mat,pi,L0,mu_index,graph,r);
        %% Algos
        Number_of_batches=4;
        [~]=Gibbs_tza(N,pi,L0,a_mat,b_mat,mu1,mu2,sigma1,sigma2,Number_of_samples,mu_index,1,Burn_in,N_thin,r,pdf_val);
        online_tza(N,pi,L0,a_mat,b_mat,mu1,mu2,sigma1,sigma2,mu_index,Number_of_batches,Number_of_samples,r,p_b,Burn_in,N_thin,pdf_val);
        Gibbs_za(N,pi,L0,a_mat,b_mat,mu1,mu2,sigma1,sigma2,Number_of_samples,mu_index,1,Burn_in,N_thin,r,pdf_val);
        online_za(N,pi,L0,a_mat,b_mat,mu1,mu2,sigma1,sigma2,mu_index,Number_of_batches,Number_of_samples,r,p_b,Burn_in,N_thin,pdf_val)
    end
end
