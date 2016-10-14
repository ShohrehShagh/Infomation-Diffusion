function [x_new,prev_x_new]=Accept(x,prev_x,x_star,prev_x_star,term,N,t_ML,t_start,p_b,pi,BL)
     R1=R_calculator(x_star,prev_x_star,N,t_ML',p_b,pi,BL);
     R2=R_calculator(x,prev_x,N,t_ML',p_b,pi,BL);
     R_data=data_likelihood_ratio(x_star,x,term,N,t_start);
     rho= prod(R1./R2)*R_data;
    if (rand <rho)
        x_new=x_star;
        prev_x_new=prev_x_star;
    else
        x_new=x;
        prev_x_new=prev_x;
    end
end
function [Ratio]=R_calculator(x,prev_x,N,t_ML,p_b,pi,BL)
    nodes=1:N;
    R=ones(N,1);
    %%
    inds1=~isnan(prev_x.t) & ~isnan(x.t);
    R(inds1)=1;
    %%
    inds2=isnan(prev_x.t) & isnan(x.t);
    A=zeros(N,1);
    B=zeros(N,1);
    if(any(inds2))
        for n1=nodes(inds2)
            candidates=(pi(n1,:).*~isnan(x.t))>0;
            AA=zeros(N,1);
            for l=find(candidates)
                AA(l)=(x.a(n1,l)/sum(x.a(n1,pi(n1,:)>0)))*(1-exp(-x.a(n1,l)))*sum((exp(-x.a(n1,l))).^((x.t(l)+1:BL)-x.t(l)));
            end
            A(n1)=sum(AA);
            B(n1)=0.5*sum(p_b*((1-p_b).^(abs((1:BL)-t_ML(n1)))))*sum(x.a(n1,candidates)/sum(x.a(n1,pi(n1,:)>0)));
        end
        R(inds2)=(1-A(inds2))./(1-B(inds2));
    end
    %%
    inds3=isnan(prev_x.t) & ~isnan(x.t) & ~isnan(x.z);%%%
    p_r2=zeros(N,1);
    pow=zeros(N,1);
    R2_num=zeros(N,1);
    R2_denom=zeros(N,1);
    if(any(inds3))
        for n1=nodes(inds3)
            p_r2(n1)=1-exp(-x.a(n1,x.z(n1)));
        end
        pow(inds3)=x.t(inds3)-x.t(x.z(inds3))-1;
        R2_num(inds3)=p_r2(inds3).*((1-p_r2(inds3)).^pow(inds3));
        R2_denom(inds3)=0.5*p_b*((1-p_b).^(abs(x.t(inds3)-t_ML(inds3))));
        R(inds3)=R2_num(inds3)./R2_denom(inds3);
        R(pow(inds3)<0 | isnan(pow(inds3)))=0;
    end
    %%
    inds4=isnan(prev_x.t) & ~isnan(x.t) & isnan(x.z);
    R(inds4)=0;
    %%
    Ratio=R;
end
function [data_likelihood]=data_likelihood_ratio(x1,x2,term,N,t_start)
    tempt1=x1.t(2:N)-t_start+1;
    tempt1(tempt1<1)=1;
    tempt1(isnan(tempt1))=size(term,2);
    %
    tempt2=x2.t(2:N)-t_start+1;
    tempt2(tempt2<1)=1;
    tempt2(isnan(tempt2))=size(term,2);
    %
    prod_term=term(sub2ind(size(term),2:N,tempt1))./term(sub2ind(size(term),2:N,tempt2));    
    if (any(prod_term==0))
        data_likelihood=0;
    else
        data_likelihood=prod(prod_term(~isnan(prod_term)));
    end
end