function [ACCres,NMIres,aveNMI,PairwiseNMI,SS,HH,ClusterRes]=progressive_EnSymNMF_v3_outputS_fixdV(W,nClass,iter,num_en,mode,gnd,V)
%%%% in this version, we introduce an alpha to contronal each simple 


% SymNMFpara.lambda=lambda; % mu has value under GNMF
SymNMFpara.k=nClass;   % set the embedding dimension = nClass
SymNMFpara.maxiter=500;
SymNMFpara.num_en=num_en;

n=length(W);
S=W;



for kk=1:iter
    SS{kk}=S;
    [H,a]=progressive_EnSymNMF_v3_subV_fixedH(S,V(:,:,:,kk),SymNMFpara);
    HH{kk}=H;
   for j=1:num_en
        [ACC{j},NMI{j},Cres{j}]=cal_ACC_NMF_symNMF_v3(H{j},gnd);
    end
    [aveNMI(kk),PairwiseNMI{kk}]=cal_aveNMI_symNMF(H);
    if mode==1 % soft model
        S=zeros(n);
        for i=1:num_en
            S=S+a(i)*H{i}*H{i}';
        end
%         S=S/num_en;
    end
    if mode==2 % soft hard model
        for j=1:num_en
            AA=repmat(max(H{j}')',1,nClass);
            H{j}(H{j}<AA)=0;
            H{j}(H{j}>0)=1;
        end
        S=zeros(n);
        for i=1:num_en
            S=S+a(i)*H{i}*H{i}';
        end
        S=S/max(max(S));
%         S=S+W;
    end
    ACCres{kk}=ACC;
    NMIres{kk}=NMI;
    ClusterRes{kk}=Cres;
end