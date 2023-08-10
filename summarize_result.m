for Group = {'All_TnRS0HC1'} 
for Hypothesis = {'C-P' 'P-C'}
    clear T
    Hyp=cell2mat(Hypothesis);
    load(strcat(cell2mat(Group),'_',Hyp))
     [ii,jj]=ind2sub([N,N],ind_upper(ind));
    for i=1:length(ii)
        x=ii(i); y=jj(i);
        T(i,:)={regions{x}, regions{y}, net_c(x,y) ,net_p(x,y)};
    end
    if length(ii)==0
        T={'None','Result'}
    else
    T(:,5)=num2cell(delta(ind)')
    A=cell2table(T)
    deltatable=table(delta)
    writetable(A,strcat(cell2mat(Group),'_',Hyp,'.txt'))
    %writetable(deltatable,strcat('Delta_',cell2mat(Group),'_',Hyp,'_',num2str(REP),'.txt'))
    end
end
end