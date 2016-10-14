function [S_parents_n1]=Sample_parent(S_changepoints,S_alphas,S_parents,pi,n1,t_start_batch,BL)
if(S_changepoints(n1)<t_start_batch)
    S_parents_n1=S_parents(n1);
else
    N=size(S_alphas,1);
    p_mat1=1-exp(-S_alphas.*pi);
    p_mat2=S_alphas./repmat(sum(S_alphas.*pi,2),1,N);
    if (isnan(S_changepoints(n1)))
        candidate_parents=find(pi(n1,:)>0);
        ft=ones(1,length(candidate_parents));
        ft(~isnan(S_changepoints(candidate_parents)))=(1-p_mat1(n1,candidate_parents(~isnan(S_changepoints(candidate_parents))))).^(BL+1);
        FCP_Z=ft.*p_mat2(n1,candidate_parents);%%%
    else
        candidate_parents=find((pi(n1,:)>0).*(S_changepoints<S_changepoints(n1)));
        FCP_Z=p_mat2(n1,candidate_parents).*p_mat1(n1,candidate_parents).*...
            ((1-p_mat1(n1,candidate_parents)).^(S_changepoints(n1)-S_changepoints(candidate_parents)-1));
    end
    if (~isempty(candidate_parents))
        S_parents_n1=candidate_parents(find(rand <= cumsum(FCP_Z)/sum(FCP_Z), 1));
        if (isempty(S_parents_n1))
            S_parents_n1=NaN;
        end
    elseif (isempty(candidate_parents))
        S_parents_n1=NaN;
    end
end