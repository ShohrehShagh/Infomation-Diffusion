clear all
clc
Number_of_batches=4;
Number_of_runs=10;
for mu_index=1:4
    for r=1:Number_of_runs
        load(sprintf('input_datas/Data_%d_%d',mu_index,r));
        T=Real.T;
        inds=Real.alpha>0;
        BL=floor(T/Number_of_batches);
        load(sprintf('Final_Results/tza_%d_%d',mu_index,r));
        Real_t=Real.t;
        Real_z=Real.z;
        Real_t(Real_t>BL)=BL+1;
        Real_z(Real_t>BL)=Inf;
        hats.z(isinf(hats.t))=Inf;
        hats.marg_t(isinf(hats.marg_t))=BL+1;
        hats.just_t(isinf(hats.just_t))=BL+1;
        Dt(1,mu_index,r)=sum(hats.t~=Real_t);
        Dt_marg(1,mu_index,r)=mean(abs(hats.marg_t-Real_t));
        Dt_just(1,mu_index,r)=mean(abs(hats.just_t-Real_t));
        Dz(1,mu_index,r)=sum(hats.z~=Real_z);
        Da(1,mu_index,r)=mean(nanmean(abs(hats.a-Real.alpha)));%./Real.alpha        
        load(sprintf('Final_Results/za_%d_%d',mu_index,r));
        hats.z(Real_t>BL)=Inf;
        Dz2(1,mu_index,r)=sum(hats.z~=Real_z);
        Da2(1,mu_index,r)=mean(nanmean(abs(hats.a-Real.alpha)));
        for i=2:Number_of_batches
            load(sprintf('Final_Results/online_tza_%d_%d_%d',mu_index,i,r))
            Real_t=Real.t;
            Real_z=Real.z;
            if(i==Number_of_batches)
                Real_t(Real_t>T)=T+1;
                Real_z(Real_t>T)=Inf;
                hats.marg_t(isinf(hats.marg_t))=T+1;
                hats.just_t(isinf(hats.just_t))=T+1;
            else
                Real_t(Real_t>i*BL)=i*BL+1;
                Real_z(Real_t>i*BL)=Inf;
                hats.marg_t(isinf(hats.marg_t))=i*BL+1;
                hats.just_t(isinf(hats.just_t))=i*BL+1;
            end
            hats.t(1)=0;
            hats.marg_t(isinf(hats.marg_t))=BL+1;
            hats.just_t(isinf(hats.just_t))=BL+1;
            Dt(i,mu_index,r)=sum(hats.t~=Real_t);
            Dz(i,mu_index,r)=sum(hats.z~=Real_z);
            Da(i,mu_index,r)=mean(nanmean(abs(hats.a-Real.alpha)));
            Dt_marg(i,mu_index,r)=mean(abs(hats.marg_t-Real_t));
            Dt_just(i,mu_index,r)=mean(abs(hats.just_t-Real_t));
            load(sprintf('Final_Results/online_za_%d_%d_%d',mu_index,i,r));
            Dz2(i,mu_index,r)=sum(hats.z~=Real_z);
            Da2(i,mu_index,r)=mean(nanmean(abs(hats.a-Real.alpha)));
        end
        Real_t=Real.t;
        Real_t(Real_t>T)=T;%Inf;
        Real_z=Real.z;
        load(sprintf('Final_Results/total/tza_%d_%d',mu_index,r))
        Dt(5,mu_index,r)=sum(hats.t~=Real_t);
        hats.marg_t(isinf(hats.marg_t))=T+1;
        hats.just_t(isinf(hats.just_t))=T+1;
        Dt_marg(5,mu_index,r)=mean(abs(hats.marg_t-Real_t));
        Dt_just(5,mu_index,r)=mean(abs(hats.just_t-Real_t));
        Dz(5,mu_index,r)=sum(hats.z~=Real_z);
        Da(5,mu_index,r)=mean(nanmean(abs(hats.a-Real.alpha)));
        load(sprintf('Final_Results/total/za_%d_%d',mu_index,r));
        Dz2(5,mu_index,r)=sum(hats.z~=Real_z);
        Da2(5,mu_index,r)=mean(nanmean(abs(hats.a-Real.alpha)));
    end
    %%
    cl=ceil(0.05*r);
    cu=floor(0.95*r);
    for j=1:Number_of_batches+1
        s_Dt=sort(Dt(j,mu_index,:),'ascend');
        s_Dz=sort(Dz(j,mu_index,:),'ascend');
        s_Da=sort(Da(j,mu_index,:),'ascend');
        A_T(j,mu_index)=mean(s_Dt);L_T(j,mu_index)=s_Dt(cl);U_T(j,mu_index)=s_Dt(cu);
        A_Z(j,mu_index)=mean(s_Dz);L_Z(j,mu_index)=s_Dz(cl);U_Z(j,mu_index)=s_Dz(cu);
        A_A(j,mu_index)=mean(s_Da);L_A(j,mu_index)=s_Da(cl);U_A(j,mu_index)=s_Da(cu);
    end
end
for j=1:Number_of_batches+1
    s_Dt_marg1=sort(Dt_marg(j,1,:),'ascend');s_Dz_marg1=sort(Dz(j,1,:),'ascend');s_Da_marg1=sort(Da(j,1,:),'ascend');
    s_Dt_just1=sort(Dt_just(j,1,:),'ascend');s_Dz_just1=sort(Dz2(j,1,:),'ascend');s_Da_just1=sort(Da2(j,1,:),'ascend');
    s_Dt_marg2=sort(Dt_marg(j,2,:),'ascend');s_Dz_marg2=sort(Dz(j,2,:),'ascend');s_Da_marg2=sort(Da(j,2,:),'ascend');
    s_Dt_just2=sort(Dt_just(j,2,:),'ascend');s_Dz_just2=sort(Dz2(j,2,:),'ascend');s_Da_just2=sort(Da2(j,2,:),'ascend');
    s_Dt_marg3=sort(Dt_marg(j,3,:),'ascend');s_Dz_marg3=sort(Dz(j,3,:),'ascend');s_Da_marg3=sort(Da(j,3,:),'ascend');
    s_Dt_just3=sort(Dt_just(j,3,:),'ascend');s_Dz_just3=sort(Dz2(j,3,:),'ascend');s_Da_just3=sort(Da2(j,3,:),'ascend');
    s_Dt_marg4=sort(Dt_marg(j,4,:),'ascend');s_Dz_marg4=sort(Dz(j,4,:),'ascend');s_Da_marg4=sort(Da(j,4,:),'ascend');
    s_Dt_just4=sort(Dt_just(j,4,:),'ascend');s_Dz_just4=sort(Dz2(j,4,:),'ascend');s_Da_just4=sort(Da2(j,4,:),'ascend');
    %%
    A_T1(j,1)=mean(s_Dt_marg1);L_T1(j,1)=s_Dt_marg1(cl);U_T1(j,1)=s_Dt_marg1(cu);
    A_T1(j,2)=mean(s_Dt_just1);L_T1(j,2)=s_Dt_just1(cl);U_T1(j,2)=s_Dt_just1(cu);
    A_Z1(j,1)=mean(s_Dz_marg1);L_Z1(j,1)=s_Dz_marg1(cl);U_Z1(j,1)=s_Dz_marg1(cu);
    A_Z1(j,2)=mean(s_Dz_just1);L_Z1(j,2)=s_Dz_just1(cl);U_Z1(j,2)=s_Dz_just1(cu);
    A_A1(j,1)=mean(s_Da_marg1);L_A1(j,1)=s_Da_marg1(cl);U_A1(j,1)=s_Da_marg1(cu);
    A_A1(j,2)=mean(s_Da_just1);L_A1(j,2)=s_Da_just1(cl);U_A1(j,2)=s_Da_just1(cu);
    %%
    A_T2(j,1)=mean(s_Dt_marg2);L_T2(j,1)=s_Dt_marg2(cl);U_T2(j,1)=s_Dt_marg2(cu);
    A_T2(j,2)=mean(s_Dt_just2);L_T2(j,2)=s_Dt_just2(cl);U_T2(j,2)=s_Dt_just2(cu);
    A_Z2(j,1)=mean(s_Dz_marg2);L_Z2(j,1)=s_Dz_marg2(cl);U_Z2(j,1)=s_Dz_marg2(cu);
    A_Z2(j,2)=mean(s_Dz_just2);L_Z2(j,2)=s_Dz_just2(cl);U_Z2(j,2)=s_Dz_just2(cu);
    A_A2(j,1)=mean(s_Da_marg2);L_A2(j,1)=s_Da_marg2(cl);U_A2(j,1)=s_Da_marg2(cu);
    A_A2(j,2)=mean(s_Da_just2);L_A2(j,2)=s_Da_just2(cl);U_A2(j,2)=s_Da_just2(cu);
    %%
    A_T3(j,1)=mean(s_Dt_marg3);L_T3(j,1)=s_Dt_marg3(cl);U_T3(j,1)=s_Dt_marg3(cu);
    A_T3(j,2)=mean(s_Dt_just3);L_T3(j,2)=s_Dt_just3(cl);U_T3(j,2)=s_Dt_just3(cu);
    A_Z3(j,1)=mean(s_Dz_marg3);L_Z3(j,1)=s_Dz_marg3(cl);U_Z3(j,1)=s_Dz_marg3(cu);
    A_Z3(j,2)=mean(s_Dz_just3);L_Z3(j,2)=s_Dz_just3(cl);U_Z3(j,2)=s_Dz_just3(cu);
    A_A3(j,1)=mean(s_Da_marg3);L_A3(j,1)=s_Da_marg3(cl);U_A3(j,1)=s_Da_marg3(cu);
    A_A3(j,2)=mean(s_Da_just3);L_A3(j,2)=s_Da_just3(cl);U_A3(j,2)=s_Da_just3(cu);
    %%
    A_T4(j,1)=mean(s_Dt_marg4);L_T4(j,1)=s_Dt_marg4(cl);U_T4(j,1)=s_Dt_marg4(cu);
    A_T4(j,2)=mean(s_Dt_just4);L_T4(j,2)=s_Dt_just4(cl);U_T4(j,2)=s_Dt_just4(cu);
    A_Z4(j,1)=mean(s_Dz_marg4);L_Z4(j,1)=s_Dz_marg4(cl);U_Z4(j,1)=s_Dz_marg4(cu);
    A_Z4(j,2)=mean(s_Dz_just4);L_Z4(j,2)=s_Dz_just4(cl);U_Z4(j,2)=s_Dz_just4(cu);
    A_A4(j,1)=mean(s_Da_marg4);L_A4(j,1)=s_Da_marg4(cl);U_A4(j,1)=s_Da_marg4(cu);
    A_A4(j,2)=mean(s_Da_just4);L_A4(j,2)=s_Da_just4(cl);U_A4(j,2)=s_Da_just4(cu);
end
%% graphs
numgroups = size(A_T1, 1);
numbars = size(A_T1, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));

h=figure(1);
%subplot(2,2,1)
b=bar(A_T1,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
%legend('marginal','just t');
legend({'$D_t(\hat{\mathbf{t}},\mathbf{t}^R)$','$D_t(\hat{\mathbf{t}^\prime},\mathbf{t}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
%title({'$D_t(\hat{\mathbf{t}},\mathbf{t}^R)$'}, 'Interpreter', 'latex');
title({'S1'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_T1(:,i), A_T1(:,i)-L_T1(:,i), U_T1(:,i)-A_T1(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off

h=figure(2);
%subplot(2,2,2)
b=bar(A_T2,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
%legend('marginal','just t');
legend({'$D_t(\hat{\mathbf{t}},\mathbf{t}^R)$','$D_t(\hat{\mathbf{t}^\prime},\mathbf{t}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
%title({'$D_t(\hat{\mathbf{t}},\mathbf{t}^R)$'}, 'Interpreter', 'latex');
title({'S2'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_T2(:,i), A_T2(:,i)-L_T2(:,i), U_T2(:,i)-A_T2(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off

h=figure(3);
%subplot(2,2,3)
b=bar(A_T3,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
%legend('marginal','just t');
legend({'$D_t(\hat{\mathbf{t}},\mathbf{t}^R)$','$D_t(\hat{\mathbf{t}}^\prime,\mathbf{t}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
%title({'$D_t(\hat{\mathbf{t}},\mathbf{t}^R)$'}, 'Interpreter', 'latex');
title({'S3'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_T3(:,i), A_T3(:,i)-L_T3(:,i), U_T3(:,i)-A_T3(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off

h=figure(4);
%subplot(2,2,4)
b=bar(A_T4,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
%legend('marginal','just t');
legend({'$D_t(\hat{\mathbf{t}},\mathbf{t}^R)$','$D_t(\hat{\mathbf{t}^\prime},\mathbf{t}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
%title({'$D_t(\hat{\mathbf{t}},\mathbf{t}^R)$'}, 'Interpreter', 'latex');
title({'S4'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_T4(:,i), A_T4(:,i)-L_T4(:,i), U_T4(:,i)-A_T4(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off
%%
numgroups = size(A_Z1, 1);
numbars = size(A_Z1, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));

h=figure(11);
%subplot(2,2,1)
b=bar(A_Z1,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
legend({'$D_z(\hat{\mathbf{z}},\mathbf{z}^R)$','$D_z(\hat{\mathbf{z}^\prime},\mathbf{z}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
title({'S1'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_Z1(:,i), A_Z1(:,i)-L_Z1(:,i), U_Z1(:,i)-A_Z1(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off

h=figure(12);
%subplot(2,2,2)
b=bar(A_Z2,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
legend({'$D_z(\hat{\mathbf{z}},\mathbf{z}^R)$','$D_z(\hat{\mathbf{z}^\prime},\mathbf{z}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
title({'S2'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_Z2(:,i), A_Z2(:,i)-L_Z2(:,i), U_Z2(:,i)-A_Z2(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off

h=figure(13);
%subplot(2,2,3)
b=bar(A_Z3,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
legend({'$D_z(\hat{\mathbf{z}},\mathbf{z}^R)$','$D_z(\hat{\mathbf{z}^\prime},\mathbf{z}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
title({'S3'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_Z3(:,i), A_Z3(:,i)-L_Z3(:,i), U_Z3(:,i)-A_Z3(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off

h=figure(14);
%subplot(2,2,4)
b=bar(A_Z4,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
legend({'$D_z(\hat{\mathbf{z}},\mathbf{z}^R)$','$D_z(\hat{\mathbf{z}^\prime},\mathbf{z}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
title({'S4'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_Z4(:,i), A_Z4(:,i)-L_Z4(:,i), U_Z4(:,i)-A_Z4(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off
%%
numgroups = size(A_A1, 1);
numbars = size(A_A1, 2);
groupwidth = min(0.8, numbars/(numbars+1.5));

h=figure(21);
%subplot(2,2,1)
b=bar(A_A1,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
%legend('tza','za');
legend({'$D_\alpha(\hat{\mathbf{\alpha}},\mathbf{\alpha}^R)$','$D_\alpha(\hat{\mathbf{\alpha}^\prime},\mathbf{\alpha}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
title({'S1'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_A1(:,i), A_A1(:,i)-L_A1(:,i), U_A1(:,i)-A_A1(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off

h=figure(22);
%subplot(2,2,2)
b=bar(A_A2,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
legend({'$D_\alpha(\hat{\mathbf{\alpha}},\mathbf{\alpha}^R)$','$D_\alpha(\hat{\mathbf{\alpha}^\prime},\mathbf{\alpha}^R)$'}, 'Interpreter', 'latex');
b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
title({'S2'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_A2(:,i), A_A2(:,i)-L_A2(:,i), U_A2(:,i)-A_A2(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off

h=figure(23);
%subplot(2,2,3)
b=bar(A_A3,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
legend({'$D_\alpha(\hat{\mathbf{\alpha}},\mathbf{\alpha}^R)$','$D_\alpha(\hat{\mathbf{\alpha}^\prime},\mathbf{\alpha}^R)$'}, 'Interpreter', 'latex');b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
title({'S3'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_A3(:,i), A_A3(:,i)-L_A3(:,i), U_A3(:,i)-A_A3(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off

h=figure(24);
%subplot(2,2,4)
b=bar(A_A4,'grouped','BarWidth',1);
set(gca,'XTickLabel',{'B1','B2','B3','B4','Batch'});%,'FontSize',18);
legend({'$D_\alpha(\hat{\mathbf{\alpha}},\mathbf{\alpha}^R)$','$D_\alpha(\hat{\mathbf{\alpha}^\prime},\mathbf{\alpha}^R)$'}, 'Interpreter', 'latex');b(1).EdgeColor = 'blue';
b(2).EdgeColor = 'red';
b(1).FaceColor = 'blue';
b(2).FaceColor = 'red';
title({'S4'});
grid on
hold on;
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
    errorbar(x, A_A4(:,i), A_A4(:,i)-L_A4(:,i), U_A4(:,i)-A_A4(:,i),'k', 'linestyle', 'none','LineWidth',3);
end
hold off


