addpath(genpath('/home/tm150/MATLAB/freezeColors'));

% You need to run './stdp test' first !

NBCELLS_S1 = 100;



if 1 == 1
    % First, look at clusters of correlation patterns between response vectors -
    % cluster together the response vectors (to many different images) that look
    % similar.  
    
%    pt = load('patterns_test.txt')'; % input patterns that gave rise to this resp, in the same format as sr (columns are cells)

    nums = 1:NBCELLS_S1;
    %sr = sr(:, sum(myresps') > 1); % We ignore columns (i.e. cells) that don't fire enough. NO!
    %nums = nums(sum(myresps') > 1);
	%sr = sr(sum(sr,2)>1e-6,:); % We only take the presentations that did elicit some response. WELL, NO.

    % Useful to save time, but also to see which dimension is the cells vs the patterns
    %sr = sr(1:2:end, :);
    %pt = pt(1:2:end, :);

    %SUMRANGE = 0;
    %for n=1:SUMRANGE
    %    sr = sr + circshift(rasterS, n);
    %end
    %if SUMRANGE>0; sr = sr(1:SUMRANGE:end,:); end
 
    %sr(:,81)=[]; sr(:,24)=[];   
    %sr = sr + 1e-5*rand(size(sr)); % To address the problem of images evoking no response at all - but then we can't use Spearman...


    %sr(sr>0) = 1;

	%sr=sr(1:3000,:);
    
    
    
    if 0
        w = load('w.txt'); 
load('wff.txt');
load ('resps_test.txt'); resps = resps_test; 
    
    sr = myresps';

        corrsrRP = corr(sr'); % Notice the transpose - columns are image presentations
        corrsrRP(isnan(corrsrRP(:))) = 0;
        lnk = linkage(double(corrsrRP), 'single', 'correlation');
        figure; [h t oRP] = dendrogram(lnk,size(corrsrRP,1));
    end
    
    
    
    
    
    sorted = corrsrRP(oRP',oRP');

    figanalyzesmall = figure; 
    
    subplot_tight(1, 3, 3, .1);
    imagesc(sorted)
    %set(gca, 'xtick', [1 100]); set(gca, 'ytick', [1 100]);
    title({'Correlation matrix' ; 'of response patterns'});
    axis square;
    colormap(hot);
    freezeColors;
    
    subplot_tight(1, 3, 2, .1);
    %imagesc(zscore(sr(oRP, :), 0, 2));
    imagesc((sr(oRP, :)));
    %set(gca, 'xtick', [1 100]); set(gca, 'ytick', [1 100]);
    %set(gca,'yaxislocation','right');
    xlabel('Cell #'); ylabel('Trial #');
    title({'Response patterns'; '(Sorted by similarity)'});
    axis square;
    colormap(hot);
 
    freezeColors; % Otherwise everythong goes grayscale !
    subplot_tight(1, 3, 1, .1);
    showw;
    title('Receptive fields');

    set(figanalyzesmall,  'Units', 'Inches', 'Resize', 'Off', 'Position', [10 10 7.2 3], 'PaperPosition', [0 0 7.2 3]);



end
