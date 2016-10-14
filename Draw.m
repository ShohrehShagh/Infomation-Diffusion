function [x]=Draw(prev_x,N,BL,p_b,t_ML,L0,t_start_batch,a_mat,b_mat,pi,type,t_Real)
% infection times
if(type==1)
    t=Sample_changepoint_online(prev_x.t,BL,p_b,t_ML,L0,t_start_batch,N);
else
    t=t_Real;
end
% \alphas
alphas=Sample_alpha_online(t,prev_x.t,prev_x.a,a_mat,b_mat,pi,N);
% parents
z=Sample_parent_online(t,prev_x.t,alphas,prev_x.z,pi,N,L0);
%t(isnan(z))=NaN;%%

x.t=t;
x.z=z;
x.a=alphas;

