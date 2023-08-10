clear all
close all

load tmp_mini %results from main.m
load c    %coordinates
N=size(c,1); 

regions={'CaudalMiddleFrontal.L',...
'EntorhinalCortex.L',...
'PostcentralGyrus.L',...
'LateralFrontalTriangularis.L',...
'SupramarginalGyrus.L',...
'Insula.L',...
'LateralOrbitofrontal.L',...
'LateralFrontalOrbitalis.L',...
'MiddleTemporal.L',...
'Pericalcarine.L',...
'Parahippocampal.L',...
'ParacentralGyrus.L',...
'MedialOrbitofrontal.L',...
'Cuneus.L',...
'InferiorTemporal.L',...
'RostralMiddleFrontal.L',...
'RostralAnteriorCingulate.L',...
'IsthmusCingulateGyrus.L',...
'InferiorOccipitalCortex.L',...
'LingualGyrus.L',...
'SuperiorParietal.L',...
'LateralFrontalOpercularis.L',...
'FusiformGyrus.L',...
'CaudalAnteriorCingulate.L',...
'SuperiorFrontalGyrus.L',...
'Precuneus.L',...
'TransverseTemporal.L',...
'PrecentralGyrus.L',...
'InferiorParietal.L',...
'PosteriorCingulate.L',...
'SuperiorTemporal.L',...
'CaudalMiddleFrontal.R',...
'EntorhinalCortex.R',...
'PostcentralGyrus.R',...
'LateralFrontalTriangularis.R',...
'SupramarginalGyrus.R',...
'Insula.R',...
'LateralOrbitofrontal.R',...
'LateralFrontalOrbitalis.R',...
'MiddleTemporal.R',...
'Pericalcarine.R',...
'Parahippocampal.R',...
'ParacentralGyrus.R',...
'MedialOrbitofrontal.R',...
'Cuneus.R',...
'InferiorTemporal.R',...
'RostralMiddleFrontal.R',...
'RostralAnteriorCingulate.R',...
'IsthmusCingulateGyrus.R',...
'InferiorOccipitalCortex.R',...
'LingualGyrus.R',...
'SuperiorParietal.R',...
'LateralFrontalOpercularis.R',...
'FusiformGyrus.R',...
'CaudalAnteriorCingulate.R',...
'SuperiorFrontalGyrus.R',...
'Precuneus.R',...
'TransverseTemporal.R',...
'PrecentralGyrus.R',...
'InferiorParietal.R',...
'PosteriorCingulate.R',...
'SuperiorTemporal.R',}; 


x=con_mat{1}; 
y=zeros(N,N); 
y(ind_upper)=delta; 
x=x.*y; 
x=x+x'; x=full(x); 
deg=sum(~~x); 
dlmwrite('edges.txt',x,'delimiter',' ','precision','%.3f')


bb=[c,ones(62,1),ones(62,1)]; 
fid=fopen('nodes.txt','w'); 
for i=1:62
    if i<62
        fprintf(fid,'%0.2f %0.2f %0.2f 1 %d %s\n',c(i,1),c(i,2),c(i,3),deg(i),regions{i}(1:end-9)); 
    else
        fprintf(fid,'%0.2f %0.2f %0.2f 1 %d %s',c(i,1),c(i,2),c(i,3),deg(i),regions{i}(1:end-9)); 
    end
end
fclose(fid);


% c=c-repmat(min(c),size(c,1),1);
% c=c./repmat(max(c),size(c,1),1); 
% r=c(:,[1,2]);
% 
% for i=1:length(ii)
%     x=ii(i); y=jj(i);
%     plot([r(x,1);r(y,1)],[r(x,2);r(y,2)],'k. ','MarkerSize',20,...
%           'MarkerEdgeColor','black');
%       if i==1
%           hold on;
%       end
%      plot([r(x,1);r(y,1)],[r(x,2);r(y,2)],'k-','Color','k','LineWidth',2);     
% end