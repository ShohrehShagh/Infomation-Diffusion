function [S_changepoints_n1]=Sample_changepoint(S_changepoints,S_parents,n1,S_alphas,BL,term_n1,pi,m)
z_n1=S_parents(n1);
if(isnan(z_n1)||isnan(S_changepoints(z_n1)))
    S_changepoints_n1=NaN;
else
    p_mat=1-exp(-S_alphas.*pi);
    p_i_z=p_mat(n1,z_n1)*(1-p_mat(n1,z_n1)).^((1:BL)-S_changepoints(z_n1)-1);
    p_i_z(1:S_changepoints(z_n1))=10^(-30);
    children=find(S_parents==n1 & ~isnan(S_changepoints));%%
    if(isempty(children) || m==1)
        p_vec=term_n1(1:BL).*p_i_z;
        child_term=1;
    else
        p_c_i=zeros(length(children),BL);
        child_term=1;
        for i=1:length(children)
            c_i=children(i);
            p_c_i(i,1:S_changepoints(c_i)-1)=p_mat(c_i,n1)*(1-p_mat(c_i,n1)).^((S_changepoints(c_i)-2):-1:0);
            p_c_i(i,S_changepoints(c_i):BL)=10^(-30);
            if(~isnan(S_changepoints(c_i)))
                child_term=0;
            else
                child_term=child_term*((1-p_mat(c_i,n1))^(BL+1));
            end
        end
        p_vec=term_n1(1:BL).*prod(p_c_i,1).*p_i_z;
    end
    FCP_C=p_vec;
    FCP_C(BL+1)=term_n1(BL+1)*((1-p_mat(n1,z_n1))^(BL+1))*child_term;
    S_changepoints_n1=find(rand <= cumsum(FCP_C)/sum(FCP_C), 1);
    if(isempty(S_changepoints_n1))
        S_changepoints_n1=NaN;
    elseif(S_changepoints_n1==BL+1)
        S_changepoints_n1=NaN;
    end
end

