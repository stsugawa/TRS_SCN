function [net,ct_res]=map_ct(ct,age,gender,icv)

%with global CT regressed out
%net=atanh(partialcorr(ct,[age,age.^2,gender,mean(ct,2)]));
%without global CT regressed out
net=atanh(partialcorr(ct,[age,gender]));
ct_res=[];
% 
% %Number of subjects
% S=size(ct,1); 
% %Number of regions
% N=size(ct,2); 
% 
% %age=age-mean(age); 
% %icv=icv-mean(icv); 
% %Regress confounds for each region
% ct_res=zeros(S,N); 
% for i=1:N
%     %[b,~,stats]=glmfit([age,age.^2,gender,icv,mean(ct,2)],ct(:,i)); %,
%     %[b,~,stats]=glmfit([age,age.^2,gender,icv],ct(:,i)); %,
%     [b,~,stats]=glmfit([age,age.^2,gender,mean(ct,2)],ct(:,i)); %,
%     %[b,~,stats]=glmfit([age,gender,icv],ct(:,i)); %,
%     %[b,~,stats]=glmfit([age,age.^2,gender],ct(:,i)); %,
%     %[b,~,stats]=glmfit([age,age.^2],ct(:,i)); %,
%     ct_res(:,i)=stats.resid+b(1); %b(1) is the constant term
%     
% end
% net=atanh(corr(ct_res)); %r to z 
% %net=(corr(ct_res)); %r to z 
% %net=cov(ct_res); 
