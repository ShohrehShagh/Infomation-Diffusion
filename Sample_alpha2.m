function [S_alphas]=Sample_alpha(S_parents,n1,S_changepoints,S_alphas,a_vec,b_vec,pi,N,BL,pdf_val)
study_inds=find(pi(n1,:)>0);
l_s_i=length(study_inds);
candidate_alphas=1:0.5:100;%%

for i=1:l_s_i
    n2=study_inds(i);
    if(a_vec(n2)*b_vec(n2)==4.5)%%
        coeff1=pdf_val(1,:);
    else
        coeff1=pdf_val(2,:);
    end
    if(isnan(S_changepoints(n1)) || isnan(S_parents(n1)))
        FCP_A=coeff1;
    else
        pow=S_changepoints(n1)-S_changepoints(S_parents(n1));
        coeff3=coeff1./((sum(S_alphas(n1,:))).^(pow+2));%%
        if(S_parents(n1)==n2)
            FCP_A=coeff3.*(candidate_alphas.^2);
        else
            FCP_A=coeff3.*((sum(S_alphas(n1,:))-S_alphas(n1,S_parents(n1))).^(pow));%%
        end
    end
    S_alphas_n1_n2=candidate_alphas(find(rand <= cumsum(FCP_A)./sum(FCP_A), 1,'first'));
    S_alphas(n1,n2)=S_alphas_n1_n2;
    S_alphas(n2,n1)=S_alphas(n1,n2);
end