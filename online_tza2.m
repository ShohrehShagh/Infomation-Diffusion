function []=online_tza2(N,pi,L0,a_mat,b_mat,mu1,mu2,sigma1,sigma2,data_num,Number_of_batches,Number_of_samples,r,p_b,Burn_in,N_thin,pdf_val)
load(sprintf('input_datas/Data_%d_%d',data_num,r));
all_data=Real.d;
T=Real.T;
BL=floor(T/Number_of_batches);
%% first block
batch=1
just_t_hat=Gibbs_tza(N,pi,L0,a_mat,b_mat,mu1,mu2,sigma1,sigma2,Number_of_samples,data_num,Number_of_batches,Burn_in,N_thin,r,pdf_val);
load(sprintf('Samples/Samples_tza_%d_%d',data_num,r));
particle_size=size(Samples_To_Save.parents,1);
%% next blocks
Samples.parents=NaN(Number_of_samples-Burn_in,N);
Samples.changepoints=NaN(Number_of_samples-Burn_in,N);
Samples.alphas=cell(Number_of_samples-Burn_in,1);
for batch=2:Number_of_batches
    batch
    %% init
    particle_set=Samples_To_Save;
    t_start_batch=BL*(batch-1)+1;
    if(batch==Number_of_batches)
        t_end_batch=T;
    else
        t_end_batch=min(BL*batch);
    end
    data=all_data(:,t_start_batch:t_end_batch);
    term=term_calculator(data,mu1,mu2,sigma1,sigma2);
    [~,t_ML]=max(term,[],2);
    for i=1:N
        if(term(i,t_ML(i))==0)
            t_ML(i)=0;
        end
    end
    just_t_hat=Just_t(N,batch,term,L0,BL,just_t_hat);
    %% Samples
    for m=1:Number_of_samples
        rand_num=ceil(particle_size*rand);
        prev_x_star.t=particle_set.changepoints(rand_num,:);
        prev_x_star.z=particle_set.parents(rand_num,:);
        prev_x_star.a=particle_set.alphas{rand_num};
        x_star=Draw(prev_x_star,N,BL,p_b,t_ML',L0,t_start_batch,a_mat,b_mat,pi);
        if(m==1)
            while (any(isnan(x_star.z) & ~isnan(x_star.t)))%%%
                rand_num=ceil(particle_size*rand);
                prev_x_star.t=particle_set.changepoints(rand_num,:);
                prev_x_star.z=particle_set.parents(rand_num,:);
                prev_x_star.a=particle_set.alphas{rand_num};
                x_star=Draw(prev_x_star,N,BL,p_b,t_ML',L0,t_start_batch,a_mat,b_mat,pi);
            end
            prev_x=prev_x_star;
            x=x_star;
        else
            [x,prev_x]=Accept(x,prev_x,x_star,prev_x_star,term,N,t_ML,t_start_batch,p_b,pi,BL);
            x=refine(x,prev_x,pi,L0,N,a_mat,b_mat,term,t_start_batch,pdf_val,m);
        end
        if (m>Burn_in)
            Samples.parents(m-Burn_in,:)=x.z;
            Samples.changepoints(m-Burn_in,:)=x.t;
            Samples.alphas{m-Burn_in}=x.a;
            if(mod(m,N_thin)==0)
                Samples_To_Save.parents((m-Burn_in)/N_thin,:)=x.z;
                Samples_To_Save.changepoints((m-Burn_in)/N_thin,:)=x.t;
                Samples_To_Save.alphas{(m-Burn_in)/N_thin}=x.a;
            end
        end
    end
    save(sprintf('Samples/Samples_online_tza_%d_%d_%d',data_num,batch,r),'Samples_To_Save','-v7.3')
end
