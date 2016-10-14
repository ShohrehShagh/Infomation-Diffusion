clc
clear all
N=20;
mu_index=2;
load(sprintf('input_datas/Data_%d_1',mu_index));
t_R=Real.t;
T=Real.T;
BL=floor(T/4);
%
Interval=[1,BL]
load(sprintf('Samples/Samples_tza_%d_1',mu_index));
t1=mode(Samples_To_Save.changepoints);
z1=mode(Samples_To_Save.parents);
[Real.t;t1]
[Real.z;z1]
%
Interval=[BL+1,2*BL]
load(sprintf('Samples/Samples_online_tza_%d_2_1',mu_index));
t2=mode(Samples_To_Save.changepoints);
z2=mode(Samples_To_Save.parents);
[Real.t;t2]
[Real.z;z2]
%
Interval=[2*BL+1,3*BL]
load(sprintf('Samples/Samples_online_tza_%d_3_1',mu_index));
t3=mode(Samples_To_Save.changepoints);
[Real.t;t3]
%
Interval=[3*BL+1,4*BL]
load(sprintf('Samples/Samples_online_tza_%d_4_1',mu_index));
t4=mode(Samples_To_Save.changepoints);
[Real.t;t4]
%
load(sprintf('Samples/total/Samples_tza_%d_1',mu_index));
t5=mode(Samples_To_Save.changepoints);
[Real.t;t5]

% NN=200;
% figure
% %subplot(3,1,1)
% a=1;
% b=0.5;
% aa=gamrnd(a,b,[1,10000]);
% [counts,centers] = hist(aa,NN);
% h = counts/sum(counts);
% bar(centers,h,'blue');
% mean(aa)
% %%
% hold on
% %subplot(3,1,2)
% a=40;
% b=0.5;
% aa=gamrnd(a,b,[1,10000]);
% [counts,centers] = hist(aa,NN);
% h = counts/sum(counts);
% bar(centers,h,'green');
% mean(aa)
% %%
% hold on
% %subplot(3,1,3)
% a=2;
% b=0.5;
% aa=gamrnd(a,b,[1,10000]);
% [counts,centers] = hist(aa,NN);
% h = counts/sum(counts);
% bar(centers,h,'red');
% mean(aa)
