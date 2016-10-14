function [S_alphas]=Sample_alpha(S_parents,n1,S_changepoints,S_alphas,a_vec,b_vec,pi,N,BL,pdf_val)
study_inds=find(pi(n1,:)>0);
l_s_i=length(study_inds);
candidate_alphas1=linspace(0.01,5,100);
candidate_alphas2=linspace(10,30,100);

for i=1:l_s_i
    n2=study_inds(i);    
    if(isnan(S_changepoints(n1)) || isnan(S_parents(n1)))
        S_alphas_n1_n2=gamrnd(a_vec(n2),b_vec(n2));
    else
        denom=sum(S_alphas(n1,:).*pi(n1,:))-S_alphas(n1,n2);
        pow=S_changepoints(n1)-S_changepoints(S_parents(n1))-1;
        if((a_vec(n2)*b_vec(n2)==0.5) || (a_vec(n2)*b_vec(n2)==1))
            coeff1=pdf_val(1,:);
            candidate_alphas=candidate_alphas1;
            exp_candidate_alphas=exp(-candidate_alphas1);
            term2=exp_candidate_alphas.^pow;
        elseif(a_vec(n2)*b_vec(n2)==20)
            coeff1=pdf_val(2,:);
            candidate_alphas=candidate_alphas2;
            exp_candidate_alphas=exp(-candidate_alphas2);
            term2=exp_candidate_alphas.^pow;
        end
        if(S_parents(n1)==n2)
            coeff2=(1-exp_candidate_alphas).*term2.*candidate_alphas;
            FCP_A=coeff1.*coeff2./(denom+candidate_alphas);
        else
            FCP_A=coeff1./(denom+candidate_alphas);
        end        
        S_alphas_n1_n2=candidate_alphas(find(rand <= cumsum(FCP_A)./sum(FCP_A), 1,'first'));
    end
    S_alphas(n1,n2)=S_alphas_n1_n2;
    S_alphas(n2,n1)=S_alphas(n1,n2);
end
