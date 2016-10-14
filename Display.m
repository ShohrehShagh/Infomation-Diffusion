clc
clear all
load outputs/Init_Value20.mat
% a_mat=a1*pi;
% a_mat=tril(a_mat,-1)'+a_mat;
% a_mat(links~=0)=a2;
warning off
for i=1:4
    i
    mu2=mu2_vec(i);
%     b2=b2_vec(i);
%     b_mat=b1*pi;
%     b_mat=tril(b_mat,-1)'+b_mat;
%     b_mat(links~=0)=b2;
    %r=0;
    for r=1:Number_of_runs
        %if (exist(sprintf('Final_Results/tza_%d_%d.mat',i,rr), 'file') && ...
                %exist(sprintf('Final_Results/za_%d_%d.mat',i,rr), 'file') && ...
                %exist(sprintf('Samples/Results_t_%d_%d.mat',i,rr), 'file'))
            %r=r+1;
            %% real values
            load(sprintf('input_datas/Data_%d_%d',i,r));
            alpha_R=Real.alpha;
            tc_R=Real.t;
            z_R=Real.z;
            sample_distance=Real.s_dis;
            
            N=size(tc_R,2);
            %% gibbs tza
            load(sprintf('Final_Results/tza_%d_%d',i,r));
            %% just t
            load(sprintf('Samples/Results_t_%d_%d',i,r));
            tc_hat_prime=S_changepoints;
            %% gibbs za
            load(sprintf('Final_Results/za_%d_%d',i,r));
            %% display parameters
            D1t(i,r)=mean(floor(abs(hats.t-tc_R)/sample_distance));
            D2t(i,r)=mean(floor(abs(tc_hat_prime-tc_R)/sample_distance));
            %
            D1z(i,r)=sum(abs(hats.z-z_R)~=0);
            D2z(i,r)=sum(abs(hats_prime.z-z_R)~=0);
            %
            D1a(i,r)=mean(mean(abs(hats.a-alpha_R)));
            D2a(i,r)=mean(mean(abs(hats_prime.a-alpha_R)));
        %end
    end
    cl=ceil(0.05*r);
    cu=floor(0.95*r);
    s_D1t=sort(D1t(i,1:r),'ascend');s_D2t=sort(D2t(i,1:r),'ascend');
    s_D1z=sort(D1z(i,1:r),'ascend');s_D2z=sort(D2z(i,1:r),'ascend');
    s_D1a=sort(D1a(i,1:r),'ascend');s_D2a=sort(D2a(i,1:r),'ascend');
    A_T(i,:)=[mean(D1t(i,1:r)),nanmean(D2t(i,1:r))];L_T(i,:)=[s_D1t(cl),s_D2t(cl)];U_T(i,:)=[s_D1t(cu),s_D2t(cu)];
    A_Z(i,:)=[mean(D1z(i,1:r)),mean(D2z(i,1:r))];L_Z(i,:)=[s_D1z(cl),s_D2z(cl)];U_Z(i,:)=[s_D1z(cu),s_D2z(cu)];
    A_A(i,:)=[mean(D1a(i,1:r)),mean(D2a(i,1:r))];L_A(i,:)=[s_D1a(cl),s_D2a(cl)];U_A(i,:)=[s_D1a(cu),s_D2a(cu)];
end
warning on
%% graphs
numgroups = size(A_T, 1);
numbars = size(A_T, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));

h=figure(1);
bar(A_T,'grouped','BarWidth',1)
set(gca,'XTickLabel',{'A','B','C','D'},'FontSize',18);
legend({'$D_t(\hat{\mathbf{t}},\mathbf{t}^R)$','$D_t(\hat{\mathbf{t}}^\prime,\mathbf{t}^R)$'}, 'Interpreter', 'latex');
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_T(:,i), A_T(:,i)-L_T(:,i), U_T(:,i)-A_T(:,i),'k', 'linestyle', 'none');
end
hold off

h=figure(2);
subplot(1,2,1)
bar(A_Z,'grouped','BarWidth',1)
set(gca,'XTickLabel',{'A','B','C','D'},'FontSize',18);
legend({'$D_z(\hat{\mathbf{z}},\mathbf{z}^R)$','$D_z(\hat{\mathbf{z}}^\prime,\mathbf{z}^R)$'}, 'Interpreter', 'latex');
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_Z(:,i), A_Z(:,i)-L_Z(:,i), U_Z(:,i)-A_Z(:,i),'k', 'linestyle', 'none');
end
hold off

subplot(1,2,2)
bar(A_A,'grouped','BarWidth',1)
set(gca,'XTickLabel',{'A','B','C','D'},'FontSize',18);
legend({'$D_{\alpha}(\hat{\mathbf{\alpha}},\mathbf{\alpha}^R)$','$D_{\alpha}(\hat{\mathbf{\alpha}}^\prime,\mathbf{\alpha}^R)$'}, 'Interpreter', 'latex');
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_A(:,i), A_A(:,i)-L_A(:,i), U_A(:,i)-A_A(:,i),'k', 'linestyle', 'none');
end
hold off

figure(3)
subplot(2,2,1)
hist(D1t(3,:))
title({'$D_{t}(\hat{\mathbf{t}},\mathbf{t}^R)$ , Scenario C'}, 'Interpreter', 'latex');
subplot(2,2,2)
hist(D1t(4,:))
title({'$D_{t}(\hat{\mathbf{t}},\mathbf{t}^R)$ , Scenario D'}, 'Interpreter', 'latex');
subplot(2,2,3)
hist(D2t(3,:))
title({'$D_{t}({\hat{\mathbf{t}}}^\prime,\mathbf{t}^R)$ , Scenario C'}, 'Interpreter', 'latex');
subplot(2,2,4)
hist(D2t(4,:))
title({'$D_{t}({\hat{\mathbf{t}}}^\prime,\mathbf{t}^R)$ , Scenario D'}, 'Interpreter', 'latex');

figure(4)
subplot(2,4,1)
hist(D1z(1,:))
title({'$D_{z}(\hat{\mathbf{z}},\mathbf{z}^R)$ , Scenario A'}, 'Interpreter', 'latex');
subplot(2,4,2)
hist(D1z(2,:))
title({'$D_{z}(\hat{\mathbf{z}},\mathbf{z}^R)$ , Scenario B'}, 'Interpreter', 'latex');
subplot(2,4,3)
hist(D1z(3,:))
title({'$D_{z}(\hat{\mathbf{z}},\mathbf{z}^R)$ , Scenario C'}, 'Interpreter', 'latex');
subplot(2,4,4)
hist(D1z(4,:))
title({'$D_{z}(\hat{\mathbf{z}},\mathbf{z}^R)$ , Scenario D'}, 'Interpreter', 'latex');
subplot(2,4,5)
hist(D2z(1,:))
title({'$D_{z}({\hat{\mathbf{z}}}^{\prime},\mathbf{z}^R)$ , Scenario A'}, 'Interpreter', 'latex');
subplot(2,4,6)
hist(D2z(2,:))
title({'$D_{z}({\hat{\mathbf{z}}}^{\prime},\mathbf{z}^R)$ , Scenario B'}, 'Interpreter', 'latex');
subplot(2,4,7)
hist(D2z(3,:))
title({'$D_{z}({\hat{\mathbf{z}}}^{\prime},\mathbf{z}^R)$ , Scenario C'}, 'Interpreter', 'latex');
subplot(2,4,8)
hist(D2z(4,:))
title({'$D_{z}({\hat{\mathbf{z}}}^{\prime},\mathbf{z}^R)$ , Scenario D'}, 'Interpreter', 'latex');

figure(5)
subplot(2,4,1)
hist(D1a(1,:))
title({'$D_{\alpha}(\hat{\mathbf{\alpha}},\mathbf{\alpha}^R)$ , Scenario A'}, 'Interpreter', 'latex');
subplot(2,4,2)
hist(D1a(2,:))
title({'$D_{\alpha}(\hat{\mathbf{\alpha}},\mathbf{\alpha}^R)$ , Scenario B'}, 'Interpreter', 'latex');
subplot(2,4,3)
hist(D1a(3,:))
title({'$D_{\alpha}(\hat{\mathbf{\alpha}},\mathbf{\alpha}^R)$ , Scenario C'}, 'Interpreter', 'latex');
subplot(2,4,4)
hist(D1a(4,:))
title({'$D_{\alpha}(\hat{\mathbf{\alpha}},\mathbf{\alpha}^R)$ , Scenario D'}, 'Interpreter', 'latex');
subplot(2,4,5)
hist(D2a(1,:))
title({'$D_{\alpha}({\hat{\mathbf{\alpha}}}^{\prime},\mathbf{\alpha}^R)$ , Scenario A'}, 'Interpreter', 'latex');
subplot(2,4,6)
hist(D2a(2,:))
title({'$D_{\alpha}({\hat{\mathbf{\alpha}}}^{\prime},\mathbf{\alpha}^R)$ , Scenario B'}, 'Interpreter', 'latex');
subplot(2,4,7)
hist(D2a(3,:))
title({'$D_{\alpha}({\hat{\mathbf{\alpha}}}^{\prime},\mathbf{\alpha}^R)$ , Scenario C'}, 'Interpreter', 'latex');
subplot(2,4,8)
hist(D2a(4,:))
title({'$D_{\alpha}({\hat{\mathbf{\alpha}}}^{\prime},\mathbf{\alpha}^R)$ , Scenario D'}, 'Interpreter', 'latex');