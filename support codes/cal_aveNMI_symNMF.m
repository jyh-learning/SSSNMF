function [aveNMI,PairwiseNMI]=cal_aveNMI_symNMF(H)
for i=1:length(H)
    [~,res{i}]=max(H{i}');
end
% PairwiseNMI=zeros(length(H))*(zeros(length(H))-1)/2;
aveNMI=0;
PairwiseNMI=0;

n=length(H);

for i=1:length(H)
    for j=1:length(H)
        try 
            temp= MutualInfo(res{i},res{j});
        catch
            temp=mean(PairwiseNMI);
        end
        PairwiseNMI((i-1)*n+j)=temp;
        aveNMI=temp+aveNMI;
    end
end
