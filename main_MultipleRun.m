
for Dataset = {'All_TRS0TnRS1' 'All_TRS0THC1' 'All_TnRS0HC1'}
    for Hypothesis = {'C-P' 'P-C'}
        REP="default";
        rng(REP);
        [ndata,text,alldata]=xlsread(strcat(cell2mat(Dataset),'.xlsx'));
        %[ndata,text,alldata]=xlsread('Destrieux_SCZ_Subtype.xlsx');
        
        %Select either NBS or FDR by commenting out one of the two lines:
        Method='NBS';
        %Method='FDR';
        
        Thresh=3; %NBS test statistic threshold (not used for FDR)
        
        %Select an alternative hypothesis by commenting out one of the two lines:
        %Hyp='P>C'; %test for patients greater than controls
        %Hyp='P<C'; %test for patients less than controls
        Hyp=cell2mat(Hypothesis);
        %With bootrapping, inference is based on a genuine t-statistic, rather than
        %a difference. The disadvantage of boostrapping is increased running time.
        Bootstrap=1; %1=yes, 0=no
        B=1000; %number of bootstraps  #before 202106, B value was 100
        
        J=50000; %number of permutations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(Method,'FDR') && J<50000
            fprintf('Number of permutations forced to J=50000. For FDR, J=%d is insufficient.\n',J);
            J=50000;
        end
        
        fprintf('Method: %s, Perms: %d, Boots: %d, Alternative Hypothesis: %s\n',...
            Method,J,B,Hyp);
        
        if strcmp(Hyp,'C-P') %strcmp(Hyp,'P<C')
            %difference if C-P
            %row numbers of controls
            ind_c=find(ndata(:,1)==0);
            %row numbers of trs
            ind_p=find(ndata(:,1)==1);
        elseif strcmp(Hyp,'P-C') %strcmp(Hyp,'P>C')
            %swap the indexes so that the difference is P-C
            %row numbers of controls
            ind_p=find(ndata(:,1)==0);
            %row numbers of trs
            ind_c=find(ndata(:,1)==1);
        end
        %gender
        gender=ndata(:,2);
        %age
        age=ndata(:,3);
        
        %regional cortical thickness estimates
        %ct=ndata(:,175:210);
        ct=ndata(:,5:length(text));
        regions=text(5:length(text));
        
        S=size(ct,1); %number of subjects
        N=size(ct,2); %number of regions
        ind_upper=find(triu(ones(N,N),1));
        
        %observed data
        net_p=map_ct(ct(ind_p,:),age(ind_p),gender(ind_p));
        net_c=map_ct(ct(ind_c,:),age(ind_c),gender(ind_c));
        
        % ind_thresh=find(abs(net_all(ind_upper))>0.2);
        % ind_upper=ind_upper(ind_thresh);
        % fprintf('Density: %0.2f, connections: %d\n',length(ind_upper)/(N*(N-1)/2),length(ind_upper));
        
        delta=net_c(ind_upper)-net_p(ind_upper); delta=delta';
        if Bootstrap
            deltaB=zeros(B,length(ind_upper));
            for i=1:B
                ind_pB=ind_p(ceil(rand(1,length(ind_p))*length(ind_p)));
                ind_cB=ind_c(ceil(rand(1,length(ind_c))*length(ind_c)));
                net_p=map_ct(ct(ind_pB,:),age(ind_pB),gender(ind_pB),icv(ind_pB));
                net_c=map_ct(ct(ind_cB,:),age(ind_cB),gender(ind_cB),icv(ind_cB));
                deltaB(i,:)=net_c(ind_upper)-net_p(ind_upper);
            end
            delta=delta./std(deltaB); %t-statistic
        end
        
        deltaX=zeros(J,length(ind_upper)); %permutations by connections
        max_sz=zeros(J,1);
        frst=0;
        for i=1:J
            ind=randperm(S);
            ind_cX=ind(1:length(ind_c));
            ind_pX=ind(length(ind_c)+1:length(ind_c)+length(ind_p));
            net_pX=map_ct(ct(ind_pX,:),age(ind_pX),gender(ind_pX),icv(ind_pX));
            net_cX=map_ct(ct(ind_cX,:),age(ind_cX),gender(ind_cX),icv(ind_cX));
            deltaX(i,:)=net_cX(ind_upper)-net_pX(ind_upper);
            if Bootstrap
                deltaB=zeros(B,length(ind_upper));
                for j=1:B
                    ind_pXB=ind_pX(ceil(rand(1,length(ind_pX))*length(ind_pX)));
                    ind_cXB=ind_cX(ceil(rand(1,length(ind_cX))*length(ind_cX)));
                    net_p=map_ct(ct(ind_pXB,:),age(ind_pXB),gender(ind_pXB),icv(ind_pXB));
                    net_c=map_ct(ct(ind_cXB,:),age(ind_cXB),gender(ind_cXB),icv(ind_cXB));
                    deltaB(j,:)=net_c(ind_upper)-net_p(ind_upper);
                end
                deltaX(i,:)=deltaX(i,:)./std(deltaB); %t-statistic
            end
            show_progress(i,J,frst);
            frst=1;
        end
        
        if strcmp(Method,'NBS')
            %NBS
            Stats.thresh=Thresh;
            Stats.alpha=0.05;
            Stats.N=N;
            Stats.test_stat=[delta;deltaX];
            Stats.size='intensity';
            [n_cnt,con_mat,pval]=NBSstats(Stats);
            ind=[];
            for i=1:n_cnt
                ind=[ind;find(con_mat{i}(ind_upper))];
            end
        elseif strcmp(Method,'FDR')
            %FDR
            %uncorrected p-values
            for i=1:length(ind_upper)
                p(i)=sum(delta(i)<=deltaX(:,i))/J;
            end
            p_uncorr=zeros(N,N); p_uncorr(ind_upper)=p;
            p_uncorr=p_uncorr+p_uncorr';
            %FDR correction
            [p_srt,ind_srt]=sort(p);
            ind=find(p_srt<0.05*[1:length(ind_upper)]/length(ind_upper));
        end
        
        if isempty(ind)
            fprintf('Null hypothesis is true\n');
        else
            %observed data
            net_p=map_ct(ct(ind_p,:),age(ind_p),gender(ind_p));
            net_c=map_ct(ct(ind_c,:),age(ind_c),gender(ind_c));
            fprintf('Null hypothesis rejected for:\n');
            if strcmp(Method,'FDR')
                [ii,jj]=ind2sub([N,N],ind_upper(ind_srt(1:ind(end))));
                for i=1:length(ii)
                    x=ii(i); y=jj(i);
                    fprintf('%s and %s | group1 %0.2f | group2 %0.2f | pval: %0.6f\n',...
                        regions{x},regions{y},net_c(x,y),net_p(x,y),p_uncorr(x,y));
                end
            elseif strcmp(Method,'NBS')
                [ii,jj]=ind2sub([N,N],ind_upper(ind));
                for i=1:length(ii)
                    x=ii(i); y=jj(i);
                    fprintf('%s and %s | group1 %0.8f | group2 %0.2f\n',...
                        regions{x},regions{y},net_c(x,y),net_p(x,y));
                end
                save tmp_mini con_mat delta ind_upper
            end
        end
        save(strcat(cell2mat(Dataset),'_',Hyp))
    end
end