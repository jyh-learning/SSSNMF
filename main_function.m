addpath(genpath('.\'))
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% load data
load('seeds.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=3; % number of class
n=210; % number of data samples

%%%  load the affinity matrix, all the affinity matrix construction methods
%%%  can be used here
load seeds_affinity_matrix
%% Progressive Ensemble
addpath('\\144.214.36.164\jyh\ProEnSymNMF')
nClass=C;
iter=10;
num_en=20;
mode=2; % 2 stands for the hard manner, while 1 stands for the soft manner

%%% init the input matrix
rng(102)
V=rand(n,C,num_en,iter);


[acc_PEsymNMF,nmi_PEsymNMF,aveNMI,PairwiseNMI,SS,HH,res_PEsymNMF]=...
    progressive_EnSymNMF_v3_outputS_fixdV(W,nClass,iter,num_en,mode,gnd,V)
%%
figure
for i=1:10
    aa(i)=mean(cell2mat(acc_PEsymNMF{i}));
end
% plot(aa)
yyaxis left
plot(aa);
%
% hold on
% figure
for i=1:10
    bb(i)=mean(PairwiseNMI{i});

end

yyaxis right
plot(bb);