function [x_star]=refine(x,prev_x,pi,L0,N,a_mat,b_mat,term,t_start,pdf_val,m,type,t_Real)
    BL=size(term,2)-1;
    x_star=x;
    % infection times
    if(type==1)
        x_star.t=Sample_changepoint_refine(x_star,prev_x,L0,term,N,t_start,BL,pi,m);
    else
        x_star.t=t_Real;
    end
    % \alphas
    x_star.a=Sample_alpha_refine(x_star,prev_x,a_mat,b_mat,pi,N,BL,pdf_val); 
    % parents
    x_star.z=Sample_parent_refine(x_star,prev_x,pi,N,L0,t_start,BL);
end
