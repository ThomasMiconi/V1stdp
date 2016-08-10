addpath(genpath('/home/tm150/MATLAB/freezeColors'));
addpath(genpath('/home/tm150/MATLAB/subplot_tight'));

% Before doing anything, you need to run ./stdp learn , and let it run for a while (should be OK when it gets to ~300K evals, 
% which should take less than a day - check the output files it generates every 100K evals )

% Then you need to run ./stdp test to generate test data.

% Each section / figure is guarded by an "if 0". Change to "if 1" to generate this particular figure. Be sure to run the commands 
% specified in the code for each figure to generate the appropriate data.


load('w.txt'); w = w(1:100, 1:100);
wlat = w;
%load ('resps_test.txt'); resps = resps_test;

NBCELLS_S1 = size(w, 1);




% w = zscore(wlat_S1);
wt = w'; w = w(~eye(size(w))); wt = wt(~eye(size(wt)));
disp('Correlation between the lateral weight matrix and its transpose (i.e. symmmetry of lateral connections):');
corrsym=corr(w, wt); disp(corrsym);
%figure; plot(w, wt, '*'); xlabel('Weight A->B'); ylabel('Weight B->A'); title('Lateral connection symmetry');
%set(gcf, 'Position', [100 100 300 300]);
%print(gcf, '-dpng', '-r0', ['connectionsymmetry_' FILENAME '.png'])

if  0
    % Is there a correlation between connection pattern similarity and mutual connection weight?
    corrs=[];
    clat = corr(wlat_S1);
    for nc1=1:NBCELLS_S1
        for nc2=nc1+1:NBCELLS_S1
            corrs = [corrs; wlat_S1(nc1, nc2) + wlat_S1(nc2, nc1)  clat(nc1, nc2) ];
        end
    end
    disp('Correlation between connection pattern similarity and mutual connection weight:');
    corrwinwmut =  (corr(corrs(:,1), corrs(:,2)));
    disp(corrwinwmut);
    %figure; plot(corrs(:,1), corrs(:,2), '*'); ylabel('Correlation between A and B''s weight vectors'); xlabel('Weight of A <-> B connection');  title('Mutual connection weight vs. overall connection similarity');
    %set(gcf, 'Position', [100 100 300 300]);
    %print(gcf, '-dpng', '-r0', ['mutualvsother_' FILENAME '.png'])
    
    
    
    
    % Is there a correlation between firing similarity and mutual connection weights?
    corrfirings = corr(rasterS); %, 'type', 'spearman');
    corrfirings = corrfirings(~eye(size(corrfirings)));
    disp('Correlation between firing similarity and mutual connection weight:');
    disp(corr(w, corrfirings)); %, 'type', 'spearman'));
end


%return;

if  0
    
    fid=fopen('w_500000.dat'); w1=fread(fid, 'double'); fclose(fid);
    fid=fopen('w_1000000.dat'); w2=fread(fid, 'double'); fclose(fid);
    disp('Corr b/w lateral weight vectors at 500K/1M: ');
    disp(corr(w1(:), w2(:)));
    
end



if  0
    % CLUSTERS
    % Are there clusters of firing patterns?
    
    
    
    if 1
        disp('Loading data (clusters)...');
        load ('resps_test.txt');
        ll = load('lastnv_test.txt');ll = ll(1:100,:);
        
        %myresps = resps_test(1:100, end-999:end);  % last 1000.
        myresps = resps_test(1:100, 1:1000);  % first 1000.
        
        sr = myresps';
        %sr=sr(1:10:1000,:);
        
        %nums = 1:100;
        
        corrsrRP = corr(sr'); % Notice the transpose - columns are image presentations
        corrsrRP(isnan(corrsrRP(:))) = 0;
        
        
        %lnk = linkage(double(corrsrRP), 'average', 'correlation');
        lnk = linkage(double(corrsrRP), 'single', 'correlation');
        %lnk = linkage(sr, 'single', 'correlation');
        
        
        % You just need to have the endrogram window pop up - there's apparently no way to recover the leaf-order from linkage  without that (optimalleaforder does something else...)
        figure; [h t oRP] = dendrogram(lnk,size(corrsrRP,1));
        clustmat = corrsrRP(oRP',oRP');
        title(pwd);
        
        % Same thing for the BAthellier data:
        load('./CC_And_Vector_Values.mat');
        pm =  mean(P, 3); % averaging over all presentations of the same stimulus...
        corrsrRP_b = corr(pm);
        corrsrRP_b(isnan(corrsrRP_b(:))) = 0;
        
        %lnk_b = linkage(double(corrsrRP_b), 'single', 'correlation');
        %lnk_b = linkage(pm, 'single', 'correlation');
        lnk_b = linkage(corrsrRP_b, 'single', 'correlation');
        %lnk_b = linkage(corrsrRP_b, 'average', 'correlation');
        figure; [h t oRP_b] = dendrogram(lnk_b,size(corrsrRP_b,1));
        clustmat_b = corrsrRP_b(oRP_b',oRP_b');
        
        
    end
    figclust = figure;
    
    
    %subplot_tight(2,2, 1, [.1 .03]);
    subplot_tight(3,2, 1, .1);
    imagesc(sr'); colormap(hot);
    set(gca,'xtick',[1 size(sr, 1)], 'ytick',[1 100]);
    %ylabel('Cell #'); xlabel('Stimulus #'); % Too much space
    text(-55, 66, 'Cell #', 'rotation', 90);
    text(350, 115, 'Stimulus #');
    cb_r = colorbar;
    limz = [min(sr(:)) max(sr(:))];
    set(cb_r,  'ytick', limz); %, 'yticklabel'); %, sprintf('%1.2f|',limz) );
    title('\bf{A} \rm{}- Population responses')
    
    subplot_tight(3,2, 2, .1);
    imagesc(sr(oRP,:)');
    set(gca,'xtick',[1 size(sr, 1)], 'ytick',[1 100]);
    cb_ss = colorbar;
    limz = [min(sr(:)) max(sr(:))];
    set(cb_ss,  'ytick', limz); %, 'yticklabel', sprintf('%1.2f|',limz) );
    %title('\textbf{B} Sorted by similarity', 'Interpreter', 'Latex')
    title('\bf{B} \rm{}- Sorted by similarity')
    
    minval=min(min(clustmat(:)), min(clustmat_b(:)));
    maxval=max(max(clustmat(:)), max(clustmat_b(:)));
    disp('Min and max values for the correlation matrices in Cluster figure:');
    disp([minval maxval]);
    
    
    subplot_tight(3,2, 3, .1);
    imagesc(clustmat);
    %axis square;
    colormap(hot); %caxis([minvalspont maxvalspont]);
    set(gca,'xtick',[]);  set(gca,'xticklabel',[])
    set(gca,'ytick',[]);  set(gca,'yticklabel',[])
    %axes('Position', [0.01 0.01 0.99 0.99], 'Visible', 'off');
    cb1 = colorbar;
    limz = [min(clustmat(:)) max(clustmat(:))];
    set(cb1,  'ytick', limz, 'yticklabel', sprintf('%1.2f|',limz) );
    %set(cb1, 'Position', [0.44 .1 .025 .35]);
    title('\bf{C} \rm{} - Correlation matrix of \bf{B}');
    
    
    %subplot_tight(1,3, 3, .0751);
    subplot_tight(3,2, 4, .1);
    %imagesc(clustmatspont);
    imagesc(clustmat_b);
    %axis square;
    colormap(hot); %caxis([minvalspont maxvalspont]);
    set(gca,'xtick',[]);  set(gca,'xticklabel',[])
    set(gca,'ytick',[]);  set(gca,'yticklabel',[])
    %axes('Position', [0.01 0.01 0.99 0.99], 'Visible', 'off');
    cb1_b = colorbar;
    limz = [min(clustmat_b(:)) max(clustmat_b(:))];
    set(cb1_b,  'ytick', limz, 'yticklabel', sprintf('%1.2f|',limz) );
    title('\bf{D} \rm{-} Mouse A1');
    %set(cb1_b, 'Position', [0.89 .1 .025 .35]);
    
    
    if  1
        % V TRACES
        %figtraces = figure;
        %imagesc(ll(:,4001:8000)) %  <--- To see some presentations with  gamma and others without...
        
        LW = 1;
        
        %STARTTIME = 8050 - 350 + 1;
        %cellz = [8 2 39];  % First 2 cells are same cluster, 3rd cell is different cluster.
        STARTTIME = 7001;
        cellz = [99 96 5];  % First 2 cells are same cluster, 3rd cell is different cluster.
        
        
        subplot_tight(3,2, 5:6, .1);
        
        
        plot(ll(cellz(1), STARTTIME:STARTTIME+2*350-1), 'r', 'linewidth', LW)
        hold on;
        plot(ll(cellz(2), STARTTIME:STARTTIME+2*350-1), 'b', 'linewidth', LW)
        plot(ll(cellz(3), STARTTIME:STARTTIME+2*350-1), 'g', 'linewidth', LW)
        legend(['Cell ' num2str(cellz(1)) ], ['Cell ' num2str(cellz(2)) ], ['Cell ' num2str(cellz(3)) ], 'orientation', 'horizontal');
        legend(gca, 'boxoff');
        hold off
        xlabel('Time (ms)');
        ylabel('V (mV)');
        axis([1 700 -100 50])
        line([350 350], [-100 50], 'linestyle', '--', 'color', 'k');
        line([250 250], [-100 50], 'linestyle', '--', 'color', 'k');
        line([600 600], [-100 50], 'linestyle', '--', 'color', 'k');
        
        %set(figtraces, 'Units', 'Inches', 'Resize', 'Off', 'Position', [10 10 6 3], 'PaperPosition', [0 0 6 3]);
        %saveas(figtraces, 'figuretraces.eps', 'epsc');
        
        title('\bf{E} \rm{-} Voltage traces');
        
    end
    
    set(findall(figclust,'-property','FontSize'),'FontSize',8)
    
    set(figclust, 'Units', 'Centimeters', 'Resize', 'Off', 'Position', [10 10 6 6], 'PaperPosition', [0 0 6 6]);
    %print(gcf, '-dtiff', '-r300', ['test.tiff'])
    print(figclust, '-depsc', '-r300', 'figurecluster.eps')
    saveas(figclust, 'figurecluster.png');
end



if 0
        % SPONTANEOUS
        % You need to run CLUSTERS first!
        
        %./stdp spont latconnmult 5.0 wie .5 wpenscale .33 altpmult .75 delayparam 5.0 

    
    if 1
        disp('Loading data (spontaneous activity)...');
        lls = load('lastnspikes_spont.txt');
        
        lls = lls(1:100,:);
        lls = lls(:, 100001:300000); BINSIZE = 50;  PROPORT = 30/100.0;
        
        
        
        lc=zeros(size(lls));
        for nn=1:BINSIZE
            lc = lc + circshift(lls, [0 nn-1]);
        end;
        lc=lc(:, 1:BINSIZE:end);
        disp(['Proportion of bins with >1 spike: ' num2str(mean(lc(:)>1))]);
        lc(lc > 1) = 1;  % Some bins may contains more than one spike.
        
        % Selecting the PROPORT % most active time bins...
        [tmp idz] = sort(sum(lc), 'descend');  idz = idz(1:ceil(PROPORT * length(idz)));
        srspont = lc(:, sort(idz))';  % The 'sort' here is to prevent the srspont to be sorted by descending activity - we want them in their original order!
        
        disp(['Number bins before culling to 1000:' num2str(size(srspont, 1))]);
        if (size(srspont, 1) > 1000)
            %    srspont = srspont(1:1000,:);
        end
        
        %nums = 1:100;
        
        disp('Computing correlation (spontaneous activity)...');
        corrsrspontRP = corr(srspont'); % Notice the transpose - columns are image presentations
        corrsrspontRP(isnan(corrsrspontRP(:))) = 0;
        
        
        %lnk = linkage(double(corrsrspontRP), 'average', 'correlation');
        lnk = linkage(double(corrsrspontRP), 'single', 'correlation');
        %lnk = linkage(srspont, 'single', 'correlation');
        
        disp('done!');
        
        
        % You just need to have the dendrogram window pop up - there's apparently no way to recover the leaf-order (and thus the similarity grouping) from linkage  without that (optimalleaforder does something else...)
        figure; [h t oRPspont] = dendrogram(lnk,size(corrsrspontRP,1));
        clustmatspont = corrsrspontRP(oRPspont',oRPspont');
        title(pwd);
        
        
        disp('Computing cross-correlation evoked/spontaneous....');
        xcorrspontevoked = corr(srspont(oRPspont,:)', sr(oRP,:)');
        disp('done!');
        
        
    end
    
    
    
    
    figspont = figure;
    
    MARGIN=.1;
    
    if 1
        %figrasterspont = figure;
        subplot_tight(2,3, 1, MARGIN);
        %imagesc( 1 - ll(:,501:1000)); colormap(gray);
        rast = lls(:,501:1000);
        indicez = find(rast > 0); [neurz, timez] = ind2sub(size(rast), indicez);
        plot(timez, neurz, '.k');
        %axis square;
        axis([1 500 1 100]);
        set(gca, 'xtick', [1 500], 'ytick', [1 100],  'ydir', 'reverse');
        %ylabel('Cell #'); xlabel('Time (ms)');  % Too far from plot...
        text(-45, 66, 'Cell #', 'rotation', 90);
        text(155, 115, 'Time (ms)');
        title('\bf{A} \rm{-} Spontaneous firing');
        
        %set(figrasterspont, 'Units', 'Inches', 'Resize', 'Off', 'Position', [10 10 2 2 ], 'PaperPosition', [0 0 2 2 ]);
        %print(figrasterspont, '-deps', '-r300', 'figurerasterspont.eps')
        %saveas(figrasterspont, 'figurerasterspont.png');
    end
    
    
    %subplot_tight(1,3, 1, .0751);
    subplot_tight(2,3, 2, MARGIN);
    imagesc(max(srspont(:)) - srspont'); colormap(gray);
    %set(gca,'xtick',[1 size(srspont,1)] ,'ytick',[1 100], 'XAxisLocation','top');
    set(gca,'xtick',[] ,'ytick',[1 100], 'XAxisLocation','top');
    text(size(srspont, 1)*1.05, 1, num2str(size(srspont, 1)));
    
    title(['\bf{B} \rm{-} ' num2str(PROPORT * 100) '% most active ' num2str(BINSIZE) 'ms time-bins']);
    %title('50ms timebins (1/3 most active)');
    
    
    
    %subplot_tight(1,3, 2, .0751);
    subplot_tight(2,3, 3, MARGIN);
    imagesc(1.0 - srspont(oRPspont,:)');
    %set(gca,'xtick',[1 size(srspont,1)],'ytick',[1 100], 'XAxisLocation','top');
    set(gca,'xtick',[],'ytick',[1 100], 'XAxisLocation','top');
    text(size(srspont, 1)*1.05, 1, num2str(size(srspont, 1)));
    title('\bf{C} \rm{-} Sorted by similarity');
    
    minvalspont=min(min(clustmatspont(:)), min(xcorrspontevoked(:)));
    maxvalspont=max(max(clustmatspont(:)), max(xcorrspontevoked(:)));
    disp('Min and max values for the correlation matrices in Spontaneous activity figure:');
    disp([minvalspont maxvalspont]);
    
    
    %subplot_tight(2,2, 3, .1);
    subplot_tight(2,3, 5, MARGIN);
    imagesc(clustmatspont);
    %axis square;
    colormap(hot); %caxis([minvalspont maxvalspont]);
    set(gca,'xtick',[],'xticklabel',[],'ytick' ,[]);
    cb1spont = colorbar;
    limz = [min(clustmatspont(:)) max(clustmatspont(:))];
    set(cb1spont,  'ytick', limz, 'yticklabel', sprintf('%1.2f|',limz) );
    %set(cb1spont, 'Position', [0.47 .1 .025 .35]);
    set(cb1spont, 'Position', [0.61 .1 .02 .35]);
    title('\bf{D} \rm{-} Correlation matrix');
    
    
    subplot_tight(2,3, 6, MARGIN);
    %imagesc(clustmatspont);
    imagesc(xcorrspontevoked');
    %axis square;
    colormap(hot); %caxis([minvalspont maxvalspont]);
    set(gca,'xtick',[], 'xticklabel',[] ,'ytick',[],'yticklabel',[]);
    %axes('Position', [0.01 0.01 0.99 0.99], 'Visible', 'off');
    cb1xcorr = colorbar;
    limz = [min(xcorrspontevoked(:)) max(xcorrspontevoked(:))];
    set(cb1xcorr,  'ytick', limz, 'yticklabel', sprintf('%1.2f|',limz) );
    %set(cb1xcorr, 'Position', [0.92 .1 .025 .35]);
    set(cb1xcorr, 'Position', [0.91 .1 .02 .35]);
    title('\bf{E} \rm{-} Cross-corr. spont./evoked');
    
    
    %caxis([minvalspont maxvalspont]);
    %set(cb1spont, 'xticklabel',  [minvalspont 0 maxvalspont], 'xtick', [minvalspont 0 maxvalspont])
    %set(cb1spont, 'Position', [0.95 .1 .025 .35]);
    
    set(findall(figspont,'-property','FontSize'),'FontSize',8)
    
    set(figspont,  'Units', 'Inches', 'Resize', 'Off', 'Position', [10 10 7.2 4], 'PaperPosition', [0 0 7.2 4]);
    %print(gcf, '-dtiff', '-r300', ['test.tiff'])
    print(figspont, '-depsc', '-r300', 'figurespont.eps')
    print(figspont, '-dpng', '-r300', 'figurespont.png')
end






if  0
    % RANDOM
    disp('Clusters for randomized-patch learning...')
    % Are there clusters of firing patterns in the randomized-patches network? (No!)
    
    
    if 1
        load ('../testRandom/resps_test.txt');
        ll_r = load('../testRandom/lastnv_test.txt');ll_r = ll_r(1:100,:);
        
        %myresps = resps_test(1:100, end-999:end);  % last 1000.
        myresps_r = resps_test(1:100, 1:1000);  % first 1000.
        
        sr_r = myresps_r';
        %sr=sr(1:10:1000,:);
        
        %nums = 1:100;
        
        corrsrRP_r = corr(sr_r'); % Notice the transpose - columns are image presentations
        corrsrRP_r(isnan(corrsrRP_r(:))) = 0;
        
        
        %lnk = linkage(double(corrsrRP), 'average', 'correlation');
        lnk_r = linkage(double(corrsrRP_r), 'single', 'correlation');
        %lnk = linkage(sr, 'single', 'correlation');
        
        
        % You just need to have the endrogram window pop up - there's apparently no way to recover the leaf-order from linkage  without that (optimalleaforder does something else...)
        figure; [h t oRP_r] = dendrogram(lnk_r,size(corrsrRP_r,1));
        clustmat_r = corrsrRP_r(oRP_r',oRP_r');
        title(pwd);
        
        
    end
    figclust_r = figure;
        
    subplot_tight(3,2, 4, .1);
    showw('../testRandom/wff.txt');
    axis normal;
    freezeColors;
    title('\bf{D} \rm{}- Receptive fields')

    
    %subplot_tight(2,2, 1, [.1 .03]);
    subplot_tight(3,2, 1, .1);
    imagesc(sr_r'); colormap(hot);
    set(gca,'xtick',[1 size(sr_r, 1)], 'ytick',[1 100]);
    %ylabel('Cell #'); xlabel('Stimulus #'); % Too much space
    text(-55, 66, 'Cell #', 'rotation', 90);
    text(350, 115, 'Stimulus #');
    cb_r_r = colorbar;
    limz = [min(sr_r(:)) max(sr_r(:))];
    set(cb_r_r,  'ytick', limz); %, 'yticklabel'); %, sprintf('%1.2f|',limz) );
    title('\bf{A} \rm{}- Population responses')
    
    subplot_tight(3,2, 2, .1);
    imagesc(sr_r(oRP_r,:)');
    set(gca,'xtick',[1 size(sr_r, 1)], 'ytick',[1 100]);
    cb_ss_r = colorbar;
    limz = [min(sr_r(:)) max(sr_r(:))];
    set(cb_ss_r,  'ytick', limz); %, 'yticklabel', sprintf('%1.2f|',limz) );
    %title('\textbf{B} Sorted by similarity', 'Interpreter', 'Latex')
    title('\bf{B} \rm{}- Sorted by similarity')
    
    minval=min(min(clustmat_r(:)));
    maxval=max(max(clustmat_r(:)));
    disp('Min and max values for the correlation matrices in Cluster figure:');
    disp([minval maxval]);
    
    
    subplot_tight(3,2, 3, .1);
    imagesc(clustmat_r);
    %axis square;
    colormap(hot); %caxis([minvalspont maxvalspont]);
    set(gca,'xtick',[]);  set(gca,'xticklabel',[])
    set(gca,'ytick',[]);  set(gca,'yticklabel',[])
    %axes('Position', [0.01 0.01 0.99 0.99], 'Visible', 'off');
    cb1 = colorbar;
    limz = [min(clustmat_r(:)) max(clustmat_r(:))];
    set(cb1,  'ytick', limz, 'yticklabel', sprintf('%1.2f|',limz) );
    %set(cb1, 'Position', [0.44 .1 .025 .35]);
    title('\bf{C} \rm{} - Correlation matrix of \bf{B}');
    
    
    
    if  1
        % V TRACES
        %figtraces = figure;
        % imagesc(ll(:,5001:8000))   %  <--- To see some presentations with  gamma and others without...
        
        LW = 1;
        
        STARTTIME = 8050 - 350 + 1;
        cellz = [8 2 39];  % First 2 cells are same cluster, 3rd cell is different cluster.
        
        %STARTTIME = 10500 - 2*350 + 1;
        %cellz = [11 17 9];  % First 2 cells are same cluster, 3rd cell is different cluster.
        
        %STARTTIME =  2*350 + 1;
        %cellz = [10 18 9];  % First 2 cells are same cluster, 3rd cell is different cluster.
        
        
        subplot_tight(3,2, 5:6, .1);
        
        
        plot(ll_r(cellz(1), STARTTIME:STARTTIME+2*350-1), 'r', 'linewidth', LW)
        hold on;
        plot(ll_r(cellz(2), STARTTIME:STARTTIME+2*350-1), 'b', 'linewidth', LW)
        plot(ll_r(cellz(3), STARTTIME:STARTTIME+2*350-1), 'g', 'linewidth', LW)
        legend(['Cell ' num2str(cellz(1)) ], ['Cell ' num2str(cellz(2)) ], ['Cell ' num2str(cellz(3)) ], 'orientation', 'horizontal');
        legend(gca, 'boxoff');
        hold off
        xlabel('Time (ms)');
        ylabel('V (mV)');
        axis([1 700 -100 50])
        line([350 350], [-100 50], 'linestyle', '--', 'color', 'k');
        line([250 250], [-100 50], 'linestyle', '--', 'color', 'k');
        line([600 600], [-100 50], 'linestyle', '--', 'color', 'k');
        
        %set(figtraces, 'Units', 'Inches', 'Resize', 'Off', 'Position', [10 10 6 3], 'PaperPosition', [0 0 6 3]);
        %saveas(figtraces, 'figuretraces.eps', 'epsc');
        
        title('\bf{E} \rm{-} Voltage traces');
        
    end
    
    
    set(findall(figclust_r,'-property','FontSize'),'FontSize',9)
    
    %set(figclust, 'Units', 'Centimeters', 'Resize', 'Off', 'Position', [10 10 12 11], 'PaperPosition', [0 0 12 11]);
        set(figclust_r, 'Units', 'Centimeters', 'Resize', 'Off', 'Position', [10 10 6 6], 'PaperPosition', [0 0 6 6]);

    %print(gcf, '-dtiff', '-r300', ['test.tiff'])
    print(figclust_r, '-depsc', '-r300', 'figurecluster_r.eps')
    saveas(figclust_r, 'figurecluster_r.png');
end







% Varying params, plus or minus 20%
if 0
    
    dirnames = {'../testPlus20pc' ; '../testMinus20pc'};
    
    % The two-step process is required because I still don't know how to get the order information from dendrogram without actually displaying the dendrogram.
    
    if 1
    for nn = 1:2
        dirname = dirnames{nn};
        disp('Loading data (varying parameters)...');
        myresps = load ([dirname '/resps_test.txt']);
        
        srV{nn} = myresps(1:100,1:1000)';
        
        corrsrRPV{nn} = corr(srV{nn}'); % Notice the transpose - columns are image presentations
        corrsrRPV{nn}(isnan(corrsrRPV{nn}(:))) = 0;
        lnk = linkage(double(corrsrRPV{nn}), 'single', 'correlation');
        figure; [h t oRPV{nn}] = dendrogram(lnk,size(corrsrRPV{nn},1));
                        sortedV{nn} = corrsrRPV{nn}(oRPV{nn}',oRPV{nn}');

    end
    end
    
    figvary = figure;
    for nn = 1:2
        dirname = dirnames{nn};

                subplot_tight(2, 3, (nn-1)*3 + 1, .1);
        showw([dirname '/wff.txt']);
        if (nn == 1)
            ylabel('Parameters +20%');
        else
            ylabel('Parameters -20%');
        end
        title('Receptive fields');
        freezeColors;
        
        subplot_tight(2, 3, (nn-1)*3 + 3, .1);
        imagesc(sortedV{nn})
        %set(gca, 'xtick', [1 100]); set(gca, 'ytick', [1 100]);
        title({'Correlation matrix' ; 'of response patterns'});
        axis square;
        colormap(hot);
        freezeColors;
        
        subplot_tight(2, 3, (nn-1)*3 + 2, .1);
        %imagesc(zscore(sr(oRP, :), 0, 2));
        imagesc((srV{nn}(oRPV{nn}, :)));
        %set(gca, 'xtick', [1 100]); set(gca, 'ytick', [1 100]);
        %set(gca,'yaxislocation','right');
        xlabel('Cell #'); ylabel('Trial #');
        title({'Response patterns'; '(Sorted by similarity)'});
        axis square;
        colormap(hot);
        cbV{nn} = colorbar;
        limz = [min(srV{nn}(:)) max(srV{nn}(:))];
        set(cbV{nn},  'ytick', limz);  % Doesn't work ??...
            set(cbV{nn}, 'Position', [0.61 (1.025 - (nn) * .45) .01 .3]);

        freezeColors; 
        

        
        set(figvary,  'Units', 'Inches', 'Resize', 'Off', 'Position', [10 10 7.2 5], 'PaperPosition', [0 0 7.2 5]);
        print(figvary, '-depsc', '-r300', 'figurevary.eps')
        print(figvary, '-dpng', '-r300', 'figurevary.png')
        
    end
end





% 400 Neurons
if 0

    if 1
        disp('Loading data (400 Neurons)...');
        myresps = load ('../test400/resps_test.txt');
        
        sr400 = myresps(1:400, :)';
        
        corrsrRP400 = corr(sr400'); % Notice the transpose - columns are image presentations
        corrsrRP400(isnan(corrsrRP400(:))) = 0;
        lnk = linkage(double(corrsrRP400), 'single', 'correlation');
        figure; [h t oRP400] = dendrogram(lnk,size(corrsrRP400,1));
    end
    
    sorted400 = corrsrRP400(oRP400',oRP400');
    
    fig400 = figure; 

    
        subplot_tight(5, 2, 5:10, .1);
    showw('../test400/wff.txt');
        ylabel('400 E Neurons');
    
    title('Receptive fields');

    freezeColors; 
    

    
    
    subplot_tight(5, 2, [2 4], .1);
    imagesc(sorted400)
    %set(gca, 'xtick', [1 100]); set(gca, 'ytick', [1 100]);
    title({'Correlation matrix' ; 'of response patterns'});
    axis square;
    colormap(hot);
    freezeColors;
    
    subplot_tight(5, 2, [1 3], .1);
    %imagesc(zscore(sr(oRP, :), 0, 2));
    imagesc((sr400(oRP400, :)));
    xlabel('Cell #'); ylabel('Trial #');
    title({'Response patterns'; '(Sorted by similarity)'});
        cb400 = colorbar;
        limz = [min(sr400(:)) max(sr400(:))];
        set(cb400,  'ytick', limz);  % Doesn't work ??...
            set(cb400, 'Position', [0.4 .65 .015 .25]);
    
    axis square;
    colormap(hot);
 

    set(fig400,  'Units', 'Inches', 'Resize', 'Off', 'Position', [10 10 7.2 6], 'PaperPosition', [0 0 7.2 6]);
        print(fig400, '-depsc', '-r300', 'figure400.eps')
        print(fig400, '-dpng', '-r300', 'figure400.png')
end





if 0
    
    % Receptive Fields  + Lateral Weights
    figRF = figure;
    MARGIN=.1;
    
    subplot_tight(1,3,1, MARGIN);
    
    showw;colormap(gray);
    title('\bf{A} \rm{}- Receptive fields');
    
    
    
    subplot_tight(1,3,2, MARGIN);
    
    load('w.txt');
    
    w = w(1:100, 1:100);
    
    imagesc(w); cbw =  colorbar(); %colormap(hot);
    set(cbw, 'Position', [0.61 .25 .01 .5]);
    
    axis square;
    set(gca,'xtick',[],'ytick',[]);
    title('\bf{B} \rm{}- Lateral weights');
    
    subplot_tight(1,3,3, MARGIN);
    xx = 0:.15:4.5;
    yy = histc(w(:), xx) ./ numel(w);  % The distribution of (lateral E-E) weights.
    bar(xx + .05, min(yy, .023));
    axis([0 4.5 0 .025]);
    set(gca, 'ytick', [0 .01 .023], 'yticklabel', {'0' '0.01' '0.89'});
    axis square;
    rectangle('Position', [-.2 .020 .5 .0003], 'clipping', 'off', 'edgecolor', 'k', 'facecolor', 'k');
    rectangle('Position', [-.2 .0203 .5 .0005], 'clipping', 'off', 'edgecolor', 'w', 'facecolor', 'w');
    rectangle('Position', [-.2 .02068 .5 .0003], 'clipping', 'off', 'edgecolor', 'k', 'facecolor', 'k');
    title('\bf{C} \rm{}- Weight distribution');
    
    
    set(findall(figRF,'-property','FontSize'),'FontSize',8)
    
    set(figRF,  'Units', 'Inches', 'Position', [10 10 7.2 3], 'PaperPosition', [0 0 7.2 3]);
    print(figRF, '-depsc', '-r300', 'figureRF.eps')
    print(figRF, '-dpng', '-r300', 'figureRF.png')
    
    
end



if 0
    
    % This portion is just here for easy copy-pasting. Not used to make the actual figures.
    %Must define STUDIEDCELL, TIMEBIN
    
    %No smoothing, just time-binned:
    ll=load('lastnspikes_pulse_nolat_noinh.txt'); tt = reshape(ll(STUDIEDCELL, :), TIMEBIN, 350/TIMEBIN, 50); z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2); hold off;
    plot(1:TIMEBIN:350, z, 'b', 'Linewidth', 2)
    ll=load('lastnspikes_pulse.txt'); tt = reshape(ll(STUDIEDCELL, :), TIMEBIN, 350/TIMEBIN, 50); z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2); hold on;
    plot(1:TIMEBIN:350, z, 'r', 'Linewidth', 2); axis([0  150 0 200]);
    
    %No Time-binning, just smoothing
    TIMEBIN=1; ll=load('lastnspikes_pulse_nolat_noinh.txt'); tt = reshape(ll(STUDIEDCELL, :), TIMEBIN, 350/TIMEBIN, 50); z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2); hold off;
    plot(1:TIMEBIN:350, smooth(z, 10), 'b', 'Linewidth', 2)
    ll=load('lastnspikes_pulse.txt'); tt = reshape(ll(STUDIEDCELL, :), TIMEBIN, 350/TIMEBIN, 50); z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2); hold on;
    plot(1:TIMEBIN:350, smooth(z, 10), 'r', 'Linewidth', 2); axis([0  150 0 200]);
    
    
end


if 1
    
    %PULSE
    figpulse = figure;
    
            % NOTE: you need to run:
        %./stdp spont latconnmult 5.0 wie .5 wpenscale .33 altpmult .75 delayparam 5.0 noinh
        %./stdp spont latconnmult 5.0 wie .5 wpenscale .33 altpmult .75 delayparam 5.0 noinh nolat
     
    % The stimuli for which we plot pulse responses (i.e. the strongest stimulus for each cluster, as returned by analyzewlat.m).
    stimlist = [120 30 50 873]; % 625   193   927    ]; %  241   629   466];    % We only use the first few ones.
    
    % For each of these, you need to run ./stdp pulse [stim #]  latconnmult 5.0 wie .5 wpenscale .33 altpmult .75 delayparam 5.0 
    % AND ./stdp pulse [stim #]  latconnmult 5.0 wie .5 wpenscale .33 altpmult .75 delayparam 5.0  noinh nolat
    
    
    
    stimlist = stimlist - 1; % Convert from Matlab to (zero-indexed) C++. 
    

    % First we plot the spontaneous/divergence
    
    %subplot_tight(1,3,1);
    subplot(1,3,1);
    
    if 1
        ll=load('lastnspikes_spont.txt');  %s=sum(ll);z=zeros(1000, 50); z(:)=s(:); s1=ll(5,:); z1=z; z1(:) = s1(:);
        ll1 = ll(1:100, 1:100000);
        ll=load('lastnspikes_spont_noinh.txt');
        ll2 = ll(1:100, 1:100000);
        ll=load('lastnspikes_spont_nolat_noinh.txt');
        ll3 = ll(1:100, 1:100000);
    end
    
    STUDIEDCELL = 10; % Must be a cell that is part of one of the clusters
    
    TIMEBIN=  8;
    
    tt = reshape(ll1(STUDIEDCELL, :), TIMEBIN, 1000/TIMEBIN, 100); % 10ms time bins, 100 of them, for 100 trials
    z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2); % sum within each time bin, divide by 10ms, multiply by 1000Hz sampling rate to get firing rate in Hz. Then take the mean along the second (trials) dimension.
    hold off;
    plot(0:TIMEBIN:59*TIMEBIN,z(1:60), 'r', 'Linewidth', 2);
    
    tt = reshape(ll2(STUDIEDCELL, :), TIMEBIN, 1000/TIMEBIN, 100); % 10ms time bins, 100 of them, for 100 trials
    z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2); % sum within each time bin, divide by 10ms, multiply by 1000Hz sampling rate to get firing rate in Hz. Then take the mean along the second (trials) dimension.
    hold on;
    plot(0:TIMEBIN:59*TIMEBIN,z(1:60), 'g', 'Linewidth', 2);
    
    tt = reshape(ll3(STUDIEDCELL, :), TIMEBIN, 1000/TIMEBIN, 100); % 10ms time bins, 100 of them, for 100 trials
    z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2); % sum within each time bin, divide by 10ms, multiply by 1000Hz sampling rate to get firing rate in Hz. Then take the mean along the second (trials) dimension.
    hold on;
    plot(0:TIMEBIN:59*TIMEBIN,z(1:60), 'b', 'Linewidth', 2);
    
    hold off;
    xlabel('Time (ms)');
    ylabel('Spikes (average 100 trials)');
    set(gca,'xtick',[0 100 200], 'ytick',[0 50 100 150 200]);
    title('Spontaneous activity');
    axis([0  200 0 200]);
    legend('Normal', 'No inhib', 'No inhib/No lat', 'location', 'East'); %, 'orientation', 'horizontal');
    legend(gca, 'boxoff');
    
    
       NBPLOTS = length(stimlist);
    
    for numstim = 1:NBPLOTS
        zestim = stimlist(numstim);
        
        %subplot_tight(NBPLOTS, 1, numstim);
        %subplot_tight(NBPLOTS/2, 3, 1+numstim + floor((numstim-1) / 2));
        subplot(NBPLOTS/2, 3, 1+numstim + floor((numstim-1) / 2));
        llNINL=load(['lastnspikes_pulse_' num2str(zestim) '_nolat_noinh.txt']);
        ll=load(['lastnspikes_pulse_' num2str(zestim) '.txt']);
        ll = ll(1:100,:);
        llNINL = llNINL(1:100,:);
        
        sumf = sum(ll,2)';
        clustercells = find(sumf>300); disp(clustercells);
        
        allspikes = [];
        allspikesNINL = [];
        
        
        
        TIMEBIN = 6;
        SMOOTHING = 0;
        
        for STUDIEDCELL = clustercells
            tt = reshape(ll(STUDIEDCELL, :), 1, 350, 50);
            tt = tt(:,1:330, :);
            
            tt = reshape(tt, TIMEBIN, 330/TIMEBIN, 50); z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2);
            if SMOOTHING
                %z = smooth(z, 10);
                z = filter(ones(1,10)/10.0, 1, z);
                
            end
            allspikes = [allspikes z];
            
            
            tt = reshape(llNINL(STUDIEDCELL, :), 1, 350, 50);
            tt = tt(:,1:330, :);
            
            tt = reshape(tt, TIMEBIN, 330/TIMEBIN, 50); z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2);
            if SMOOTHING
                %z = smooth(z, 10);
                z = filter(ones(1,10)/10.0, 1, z);
            end
            allspikesNINL = [allspikesNINL z];
        end
        
        
        
        avgNINL = mean(allspikesNINL);
        % hold off; plot(mean(allspikes ./ (1.0 + repmat(avgNINL, 350, 1)), 2), 'r') ; hold on; plot(mean(allspikesNINL ./ (1.0 + repmat(avgNINL, 350, 1)), 2), 'b'); axis([0  150 0 10]); % Not much difference...
        hold off; plot(1:TIMEBIN:330, mean(allspikes, 2), 'r', 'Linewidth', 2) ; hold on; plot(1:TIMEBIN:330, mean(allspikesNINL, 2), 'b', 'Linewidth', 2); axis([0  150 0 200]);
        hold off;
        
        line([100 100], [0 150], 'linestyle', '--', 'color', 'k');
        set(gca,'xtick',[0 100 120], 'ytick',[0 50 100 150]);
        if(numstim >2)
            xlabel('Time (ms)');
        end
        if(mod(numstim, 2) == 1)
            ylabel('Firing rate');
        end
        title(['Cluster ' num2str(numstim)]);
        axis([1 120 0 150]);
    end
    
    set(findall(figpulse,'-property','FontSize'),'FontSize',9)
    
    %set(figpulse, 'Units', 'Inches', 'Resize', 'Off', 'Position', [10 10 5 3], 'PaperPosition', [0 0 5 3]);
    set(figpulse, 'Units', 'Inches',  'Position', [10 10 7.2 3], 'PaperPosition', [0 0 7.2 3]);
    %saveas(figpulse, 'figurepulse.png');
    print(figpulse, '-depsc', '-r300', 'figurepulse.eps')
    print(figpulse, '-dpng', '-r300', 'figurepulse.png')
    
    
end


     %3     1     8     4     6     9     2     5     7
  % 903   120    10   873    30   859   633   216    50
  %  36    29    34    31    31    32    34    31    28

% MIXING

%You must run these:
% ./stdp mix 633 30 latconnmult 5.0 wie .5 wpenscale .33 altpmult .75 delayparam 5.0
% ./stdp mix 633 30 latconnmult 5.0 wie .5 wpenscale .33 altpmult .75 delayparam 5.0 nolat
% ./stdp mix 633 30 latconnmult 5.0 wie .5 wpenscale .33 altpmult .75 delayparam 5.0 nospike
% ./stdp mix 633 30 latconnmult 5.0 wie .5 wpenscale .33 altpmult .75 delayparam 5.0 nospike nolat


    %The stimuli we mix are the stimuli evoking the largest responses of two different (arbitrarily picked) clusters
    NUM1 = 633; NUM2 = 30;  
    



    
if  0

    % Note : noelat = no E-E connections. nolat = no connections at all
    % (including inhibitory). Both make the curves smoother (reduce the
    % sharpness of the transition), but nolat even more (unsurprisingly).
    
    disp('Loading data (Mixing experiments...)');
    rl = load (['resps_mix_' num2str(NUM1 - 1) '_' num2str(NUM2 - 1) '.txt']);   % -1 converts from Matlab 1-indexing to C++ 0-indexing
    rnl = load (['resps_mix_' num2str(NUM1 - 1) '_' num2str(NUM2 - 1) '_nolat.txt']); 
    rl = rl(1:100, :); rnl = rnl(1:100, :);
    %rl = load ('respssumv_mix.txt'); rnl = load ('respssumv_mix_nolat.txt');
    %rl = load ('respssumv_mix_nospike.txt'); rnl = load ('respssumv_mix_nolat_nospike.txt');
    
    NBMIXES = size(rl, 2)/3;
    
    MARGIN = .1;
    
    figmix=figure;
    
    colormap(hot);
    ylimz = [0 max(rl(:))];
    
    subplot_tight(1,4,1, MARGIN);
    imagesc(rl(:,1:NBMIXES)); set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'100%S2', '100%S1'}, 'ytick',[1 100]); ylabel('Cell #');
    title('\bf{A} \rm{}- Mix Stim.1 / Stim.2');
    caxis(ylimz);
    
    
    subplot_tight(1,4,2, MARGIN); imagesc(rl(:, (NBMIXES + 1):2*NBMIXES)); set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'0%S1', '100%S1'}, 'ytick',[]);
    title('\bf{B} \rm{}- Stim.1 only');
    caxis(ylimz);
    
    subplot_tight(1,4,3, MARGIN); imagesc(rl(:, (2*NBMIXES + 1):end)); set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'100%S2', '0%S2'}, 'ytick',[]);
    title('\bf{C} \rm{}- Stim.2 only');
    caxis(ylimz);
    
    
    cbr = colorbar;
    set(cbr,  'ytick', ylimz); %, 'yticklabel', sprintf('%1.2f|',limz) );
    set(cbr, 'Position', [0.47 .15 .020 .70]);
    
    
    subplot_tight(1,4,4, MARGIN);
    
    bs=[]; bs_nl=[];
    % We regress the response to each mix upon the responses at 100% Stim1 and 100% Stim2, and store the resulting regression weights - both for fulll-network and without-lat. conn.
    for nn=1:NBMIXES
        [b bint] = regress(rl(:, nn)./norm(rl(:,nn)), [rl(:,1)./norm(rl(:,1)), rl(:, NBMIXES)./norm(rl(:,NBMIXES))]);
        [bnl bint] = regress(rnl(:, nn)./norm(rnl(:,nn)), [rnl(:,1)./norm(rnl(:,1)), rnl(:, NBMIXES)./norm(rnl(:,NBMIXES))]);
        bs = [bs b]; bs_nl = [bs_nl bnl];
    end
    set(figmix,'defaultAxesColorOrder',[1 0 0;0 0 1]);
    hold off; plot([1:NBMIXES], bs', 'LineWidth', 2);
    axis([1 NBMIXES -.05 1.05]); hold on;
    plot([1:NBMIXES], bs_nl', ':', 'LineWidth', 2);
    axis([1 NBMIXES -.05 1.05]);
    set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'100%S2' '100% S1'}, 'ytick',[]);
    title('\bf{D} \rm{}- Regression');
    legend('Beta S2', 'Beta S1', 'Beta S2 (No lat.)', 'Beta S1 (No lat.)',  'location', 'southoutside');
    legend boxoff;
    
    
    set(figmix, 'Units', 'Inches',  'Position', [10 10 7.2 3], 'PaperPosition', [0 0 7.2 3]);
    
    print(figmix, '-depsc', '-r300', 'figuremix.eps')
    print(figmix, '-dpng', '-r300', 'figuremix.png')
    
    
end


% MIXING - potentials

if  1
    disp('Loading data (Mixing experiments with potentials...)');


    
    rl = load (['respssumv_mix' num2str(NUM1 - 1) '_' num2str(NUM2 - 1) '.txt']);   % -1 converts from Matlab 1-indexing to C++ 0-indexing
    rnl = load (['respssumv_mix' num2str(NUM1 - 1) '_' num2str(NUM2 - 1) '_nolat.txt']); 

    rl = rl(1:100, :) ./ 350; rnl = rnl(1:100, :) ./350;  % We divide by 350 ms (the duration of each presentation) to get the average V for each presentation
    
    ylimz = fix( [min(rl(:)) max(rl(:))] );
    
    NBMIXES = size(rl, 2)/3;
    figmix_v=figure;
    
    colormap(hot);
    
    subplot_tight(2,4,1, MARGIN);
    imagesc(rl(:,1:NBMIXES)); set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'100%S2', '100%S1'}, 'ytick',[1 100]); ylabel('Cell #');
    title('\bf{A} \rm{}- Mix Stim.1 / Stim.2');
    caxis(ylimz);
    
    subplot_tight(2,4,2, MARGIN); imagesc(rl(:, (NBMIXES + 1):2*NBMIXES)); set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'0%S1', '100%S1'}, 'ytick',[]);
    title('\bf{B} \rm{}- Stim.1 only');
    caxis(ylimz);
    
    subplot_tight(2,4,3, MARGIN); imagesc(rl(:, (2*NBMIXES + 1):end)); set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'100%S2', '0%S2'}, 'ytick',[]);
    title('\bf{C} \rm{}- Stim.2 only');
    caxis(ylimz);
    
    cbv = colorbar;
    set(cbv,  'ytick', ylimz) ; %, 'yticklabel', sprintf('%1.2f|', ylimz) );
    %set(cbv, 'Position', [0.47 .15 .020 .70]);
    set(cbv, 'Position', [0.47 .60 .020 .25]);
    
    
    subplot_tight(2,4,4, MARGIN);
    bs=[]; bs_nl=[];
    for nn=1:NBMIXES
        [b bint] = regress(rl(:, nn)./norm(rl(:,nn)), [rl(:,1)./norm(rl(:,1)), rl(:, NBMIXES)./norm(rl(:,NBMIXES))]);
        [bnl bint] = regress(rnl(:, nn)./norm(rnl(:,nn)), [rnl(:,1)./norm(rnl(:,1)), rnl(:, NBMIXES)./norm(rnl(:,NBMIXES))]);
        bs = [bs b]; bs_nl = [bs_nl bnl];
    end
    set(figmix_v,'defaultAxesColorOrder',[1 0 0;0 0 1]);
    hold off;
    plot([1:NBMIXES], bs', 'LineWidth', 2); axis([1 NBMIXES -.05 1.05]);
    %hold on; plot([1:NBMIXES], bs_nl', ':', 'LineWidth', 2);
    axis([1 NBMIXES -.05 1.05]); set(gca,'xtick',[],'ytick',[]);
    
    
    title('\bf{D} \rm{}- Regression');
    %legend('Beta S2', 'Beta S1', 'Beta S2 (No lat.)', 'Beta S1 (No lat.)',  'location', 'southoutside');
    legend('Beta S2', 'Beta S1',  'location', 'southoutside');
    legend boxoff;
    
    
    %set(figmix_v, 'Units', 'Inches',  'Position', [10 10 7.2 2.4], 'PaperPosition', [0 0 7.2 2.4]);
    %print(figmix_v, '-depsc', '-r300', 'figuremix_v.eps')
    %print(figmix_v, '-dpng', '-r300', 'figuremix_v.png')
    
    
    
    
    % Same, with spiking blocked / no spikes:
    
    rl = load (['respssumv_mix' num2str(NUM1 - 1) '_' num2str(NUM2 - 1) '_nospike.txt']);   % -1 converts from Matlab 1-indexing to C++ 0-indexing
    %rnl = load (['respssumv_mix_' num2str(NUM1 - 1) '_' num2str(NUM2 - 1) '_nolat_nospike.txt']); 
    
    rl = rl(1:100, :) ./ 350; %rnl = rnl(1:100, :) ./350;  % We divide by 350 ms (the duration of each presentation) to get the average V for each presentation
    
    ylimz = fix( [min(rl(:)) max(rl(:))] );
    NBMIXES = size(rl, 2)/3;
    
    %figmix_v_nospike=figure;
    
    colormap(hot);
    
    subplot_tight(2,4,5, MARGIN);
    imagesc(rl(:,1:NBMIXES)); set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'100%S2', '100%S1'}, 'ytick',[1 100]); ylabel('Cell #');
    title('\bf{E} \rm{}- Mix Stim.1 / Stim.2');
    caxis(ylimz);
    
    subplot_tight(2,4,6, MARGIN); imagesc(rl(:, (NBMIXES + 1):2*NBMIXES)); set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'0%S1', '100%S1'}, 'ytick',[]);
    title('\bf{F} \rm{}- Stim.1 only');
    caxis(ylimz);
    
    subplot_tight(2,4,7, MARGIN); imagesc(rl(:, (2*NBMIXES + 1):end)); set(gca,'xtick',[1 NBMIXES], 'xticklabel', {'100%S2', '0%S2'}, 'ytick',[]);
    title('\bf{G} \rm{}- Stim.2 only');
    caxis(ylimz);
    
    cbvns = colorbar;
    set(cbvns,  'ytick', ylimz) ; %, 'yticklabel', sprintf('%1.2f|', ylimz) );
    %set(cbvns, 'Position', [0.47 .15 .020 .70]);
    set(cbvns, 'Position', [0.47 .15 .020 .25]);
    
    subplot_tight(2,4,8, MARGIN);
    
    bs=[]; bs_nl=[];
    for nn=1:NBMIXES
        [b bint] = regress(rl(:, nn)./norm(rl(:,nn)), [rl(:,1)./norm(rl(:,1)), rl(:, NBMIXES)./norm(rl(:,NBMIXES))]);
        %[bnl bint] = regress(rnl(:, nn)./norm(rnl(:,nn)), [rnl(:,1)./norm(rnl(:,1)), rnl(:, NBMIXES)./norm(rnl(:,NBMIXES))]);
        bs = [bs b]; bs_nl = [bs_nl bnl];
    end
    %set(figmix_v_nospike,'defaultAxesColorOrder',[1 0 0;0 0 1]);
    hold off; plot([1:NBMIXES], bs', 'LineWidth', 2); axis([1 NBMIXES -.05 1.05]);
    %hold on; plot([1:NBMIXES], bs_nl', ':', 'LineWidth', 2);
    axis([1 NBMIXES -.05 1.05]); set(gca,'xtick',[],'ytick',[]);
    
    title('\bf{H} \rm{}- Regression');
    %legend('Beta S2', 'Beta S1', 'Beta S2 (No lat.)', 'Beta S1 (No lat.)',  'location', 'southoutside');
    legend('Beta S2', 'Beta S1',   'location', 'southoutside');
    legend boxoff;
    
   set(figmix_v , 'Units', 'Inches',  'Position', [10 10 7.2 4.8], 'PaperPosition', [0 0 7.2 4.8]);
    %print(figmix_v_nospike, '-depsc', '-r300', 'figuremix_v_nospike.eps')
    %print(figmix_v_nospike, '-dpng', '-r300', 'figuremix_v_nospike.png')
    print(figmix_v , '-depsc', '-r300', 'figuremix_v.eps')
    print(figmix_v , '-dpng', '-r300', 'figuremix_v.png')
    
end


return;

% Previous version:
load resps_mix.txt; r= resps_mix; cc = corr(r); v = cc(:, [60 61]);
load resps_mix_nolat.txt; rnolat= resps_mix_nolat; ccnl = corr(rnolat); vnl = ccnl(:, [60 61]);
figmix = figure;
%set(0,'DefaultAxesFontSize',14); subplot_tight(141); imagesc(r(:,1:30)); title('Mixed stimuli,  0% / 100% to 100% / 0%'); subplot_tight(142); imagesc(r(:,31:60)); title('Stimulus 1 alone, 0% to 100%'); subplot_tight(143); imagesc(r(:,61:end)); title('Stimulus 2 alone,  100% to 0%'); subplot_tight(144); plot(1:30,v(1:30,1), 'r', 1:30, v(1:30,2), 'b', 'Linewidth', 2); axis([1 30 -.3 1.05]);
subplot_tight(141); imagesc(r(:,1:30)); set(gca,'xtick',[]);  set(gca,'xticklabel',[]); set(gca,'ytick',[]);  set(gca,'yticklabel',[]);
subplot_tight(142); imagesc(r(:,61:end)); set(gca,'xtick',[]);  set(gca,'xticklabel',[]); set(gca,'ytick',[]);  set(gca,'yticklabel',[]);
subplot_tight(143); imagesc(r(:,31:60)); set(gca,'xtick',[]);  set(gca,'xticklabel',[]); set(gca,'ytick',[]);  set(gca,'yticklabel',[]);
subplot_tight(144); plot(1:30,v(1:30,1), 'r', 1:30, v(1:30,2), 'b', 'Linewidth', 2); axis([1 30 -.2 1.05]); set(gca,'xtick',[]);  set(gca,'xticklabel',[]); set(gca,'ytick',[]);  set(gca,'yticklabel',[]);

figmix.Units = 'inches';
figmix.Position = [1 1 7 3];

%figmixnolat = figure;
%set(0,'DefaultAxesFontSize',14); subplot_tight(141); imagesc(rnolat(:,1:30)); title('Mixed stimuli,  0% / 100% to 100% / 0%'); subplot_tight(142); imagesc(rnolat(:,31:60)); title('Stimulus 1 alone, 0% to 100%'); subplot_tight(143); imagesc(rnolat(:,61:end)); title('Stimulus 2 alone,  100% to 0%'); subplot_tight(144); plot(1:30,vnl(1:30,1), 'm', 1:30, vnl(1:30,2), 'c', 'Linewidth', 2); axis([1 30 -.3 1.05]);







