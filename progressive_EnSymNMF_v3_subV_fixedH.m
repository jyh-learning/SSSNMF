function [H,a,obj]=progressive_EnSymNMF_v3_subV_fixedH(S,V,para)%L=D-O;
% work for both NMF and GNMF
% lambda=para.lambda;
k=para.k;
% p=para.p;
maxiter=para.maxiter;
num_en=para.num_en;


% d=size(S,1);
n=size(S,2);

for i=1:num_en
    H{i}=V(:,:,i);
end
a=rand(num_en,1);

for i=1:num_en
    obj_temp(i)=a(i)*a(i)*sum(sum((S-H{i}*H{i}').^2));
end
obj(1)=sum(obj_temp);



for iter=1:maxiter
    parfor kk=1:num_en
        NUM{kk}=S*H{kk};
        DEN{kk}=H{kk}*H{kk}'*H{kk}+eps;
        E{kk}=(NUM{kk}./DEN{kk});
        H{kk}=H{kk}.*(E{kk}.^(1/4));
        h(kk)=norm(S-H{kk}*H{kk}').^2;
%         obj_temp(kk)=a(kk)*h(kk);
    end
    sumH=sum(1./h);
    a=(1./h')/sumH;
    obj(iter+1)=h*(a.^2);


    disp(['the ', num2str(iter), ' obj is ', num2str(obj(iter))]);
    if (iter>300 && abs(obj(iter+1)-obj(iter))<10^-3)
        break;
    end
end