addpath(genpath('/home/tm150/MATLAB/freezeColors'));

% This analyzes the output of the network and generate various figures. It needs shoow.m

% First you need to run './stdp learn' to build the network. Then you need to run './stdp test' to generate the output data (with frozen weights).


NBCELLS_S1 = 100;

w = load('w.txt'); 
load('wff.txt');
load ('resps_test.txt'); resps = resps_test; 
%load ('resps.txt'); 



% If you want to only use the excitatory subnetwork - which you should ! 
resps = resps(1:NBCELLS_S1, :);
w = w(1:NBCELLS_S1, 1:NBCELLS_S1);
wff = wff(1:NBCELLS_S1, :);


wlat = w;
lastresp = max(find(sum(resps)>0));
if lastresp > 1000
    myresps = resps(:, lastresp-1000:lastresp);
else
    myresps = resps;
end

r = load('resps.txt'); re = r(1:100,:); ri = r(101:end,:);
disp(['Max R-Exc: ' num2str(max(re(:))) ', median-max: ' num2str(median(max(re))) ', sum R-Exc: ' num2str(sum(re(:)))])


% Cossell et al. 2015: 7% most correlated firings have 50% of the total weight. 12% most correlated RFs have 50% of the total weight.
% What is the % of cells that represent 50% of all the weights?
% NOTE : This should not include inhibitory connections !
wvect = w(eye(size(w)) ~= 1);
wvect_s = sort(wvect, 1, 'descend'); % Ignore the diagonal (self) weights (which are forced to zero anyway). Also, sort from highest weight to lowest weight!
csw = cumsum(wvect_s); csw = csw / max(csw);
loc50w = min(find(csw > .5)); prop50w = loc50w / numel(csw); disp([num2str(prop50w) ' is the proportion of pairs that contain 50% of the total weight.']);
disp( [ num2str(csw(floor(.05*numel(csw)))) ' is the proportion of weights contained in the 5% strongest connections.']);
disp( [ num2str(csw(floor(.1*numel(csw)))) ' is the proportion of weights contained in the 10% strongest connections.']);

r = myresps;
cc = corr(r'); cc = cc - eye(size(cc));
ccvect = (cc(eye(size(cc)) ~= 1));
[x ccorder] = sort(ccvect, 1, 'descend'); % ccorder contains the order of the pairs from most-correlated responses to least-correlated response
wvect_cc = wvect(ccorder); csw = cumsum(wvect_cc); csw = csw / max(csw);
loc50w_cc = min(find(csw > .5)); prop50w_cc = loc50w_cc / numel(csw); disp([num2str(prop50w_cc) ' is the proportion of most correlated-response pairs that contain 50% of the total weight.']);

ccwff = corr(wff');  ccwff = ccwff - eye(size(ccwff));
ccwffvect = (ccwff(eye(size(ccwff)) ~= 1));
ccwff_spearman = corr(wff', 'type', 'spearman');  ccwff_spearman = ccwff_spearman - eye(size(ccwff_spearman));
ccwffvect_spearman = (ccwff_spearman(eye(size(ccwff_spearman)) ~= 1));
[x ccwfforder] = sort(ccwffvect, 1, 'descend'); % ccwfforder contains the order of the pairs from most-correlated responses to least-correlated response
wvect_ccwff = wvect(ccwfforder); csw = cumsum(wvect_ccwff); csw = csw / max(csw);
loc50w_ccwff = min(find(csw > .5)); prop50w_ccwff = loc50w_ccwff / numel(csw); disp([num2str(prop50w_ccwff) ' is the proportion of most correlated-RF pairs that contain 50% of the total weight.']);

corrsmostconn10pc = ccwffvect_spearman(wvect >= wvect_s(floor(.1*numel(wvect_s))));
disp(['Spearman correlations among the 10% most strongly connected (which account for ' num2str( csw(numel(corrsmostconn10pc))) ' of the total connection weight): Mean ' num2str(mean(corrsmostconn10pc)) ', median ' num2str(median(corrsmostconn10pc)) ]);






% w = zscore(wlat_S1);
wt = w'; w = w(~eye(size(w))); wt = wt(~eye(size(wt)));
disp('Correlation between the lateral weight matrix and its transpose (i.e. symmmetry of lateral connections):');
corrsym=corr(w, wt); disp(corrsym);
%figure; plot(w, wt, '*'); xlabel('Weight A->B'); ylabel('Weight B->A'); title('Lateral connection symmetry'); 
%set(gcf, 'Position', [100 100 300 300]);
%print(gcf, '-dpng', '-r0', ['connectionsymmetry_' FILENAME '.png'])   

if 1 == 0
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

% Are there clusters of firing patterns?



if 1 == 1
    % First, look at clusters of correlation patterns between response vectors -
    % cluster together the response vectors (to many different images) that look
    % similar.  
    
%    pt = load('patterns_test.txt')'; % input patterns that gave rise to this resp, in the same format as sr (columns are cells)
    
    sr = myresps';
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
    
    
    corrsrRP = corr(sr'); % Notice the transpose - columns are image presentations
    corrsrRP(isnan(corrsrRP(:))) = 0;
    lnk = linkage(double(corrsrRP), 'single', 'correlation');


    figure; [h t oRP] = dendrogram(lnk,size(corrsrRP,1));
    sorted = corrsrRP(oRP',oRP');
    %showw; % Should appear in the dendrogram window
    title(pwd);
    figcorrpres = figure; 
    
    zz=4:6; tt = zz; for n=2:3; zz = zz + 9; tt = [tt zz ]; end;
    subplot(6, 9, tt);
    imagesc(sorted)
    %set(gca, 'xtick', [1 100]); set(gca, 'ytick', [1 100]);
    title({'Correlation matrix' ; 'of response patterns'});
    freezeColors;
    zz=7:9; tt = zz; for n=2:3; zz = zz + 9; tt = [tt zz ]; end;
    subplot(6, 9, tt);
    %imagesc(zscore(sr(oRP, :), 0, 2));
    imagesc((sr(oRP, :)));
    %set(gca, 'xtick', [1 100]); set(gca, 'ytick', [1 100]);
    set(gca,'yaxislocation','right');
    xlabel('Cell #'); ylabel('Trial #');
    title({'Response patterns'; ' (not contrast-enhanced)'});
    freezeColors;
    zz=1:3; tt = zz; for n=2:3; zz = zz + 9; tt = [tt zz ]; end;
    subplot(6, 9, tt);
    imagesc(wlat);
    zz = pwd; title(zz(end-5:end));
%    imagesc(pt(oRP, :));
%    set(gca, 'xtick', [1 100]); set(gca, 'ytick', [1 100]);
%    xlabel(pwd);
%    zz = 63:68;  tt = zz; for n=2:2; zz = zz + 10; tt = [tt zz ]; end;  
    freezeColors;
    tt = 31:33;
    subplot(6, 9, tt);
    imagesc(sum(sr(oRP, :), 2)');
    xlabel('Total activity for each pattern');
    freezeColors;
    zz = 28:30; tt = zz; for n=2:3; zz = zz + 9; tt = [tt zz ]; end;
    subplot(6, 9, tt);
    showw;




    % We identify a group of responses as belonging to specific clusters. Note that some of the parameters may need adjustment - e.g. the distance cutoff for clustering, or the minimum total firing limit for selection
    % _s = "sorted" (by the previous sorting process); _f = "filtered" (by being in a cluster, AND having high response)
    
    r_s = sr(oRP,:);
    nums_s = oRP;
    c_s = corrsrRP(oRP,:);
    s_s = sum(r_s,2);
    inclust = sum(c_s>.9,2) > 40; % It's "in a cluster" if it correlates >.9 with more than 40 other response patterns

    selcrit = inclust == 1 & s_s > median(s_s); % Only those of firing above median are selected.

    r_f = r_s(selcrit, :);
    nums_f = nums_s(selcrit);
    clustnums  = clusterdata(r_f,'criterion', 'distance', 'distance', 'correlation', 'cutoff', .3);
        
    maxresp_clusts = [];
    nummaxresp_clusts = [];
    for nn=unique(clustnums, 'stable')'
        allthisclust=find(clustnums == nn); nums_allthisclust = nums_f(allthisclust);
        %[rmax_thisclust i] = max(sum(r_f(allthisclust, :), 2));
        [rmax_thisclust i] = max(max(r_f(allthisclust, :), [], 2));
        nummax_thisclust = nums_allthisclust(i); 
        maxresp_clusts = [maxresp_clusts rmax_thisclust];
        nummaxresp_clusts = [nummaxresp_clusts nummax_thisclust];
    end
    disp([unique(clustnums, 'stable')'; lastresp - 1000 + nummaxresp_clusts; maxresp_clusts]); disp(['Last Resp: ' num2str(lastresp)]);

    % Now r_f contains the selected response patterns, nums_f contains their appearance number (1 to 1001), and clustnums the cluster they've been assigned to.




    % Summing the test patterns, and the responses to them by clusters.
    %pts = pt(oRP,:); 
    %srs = sr(oRP,:);
    %pl = load('patterns_learn.txt');
    %sumpts = ([(sum(pts(1:30,:).^4)); sum(pts(43:55,:).^4) ; sum(pts(64:73, :).^4)]);
    %sumsrs = ([(sum(srs(1:30,:).^1)); sum(srs(43:55,:).^1) ; sum(srs(64:73, :).^1)]);
    % imagesc(zscore(corr((pts').^4, pl, 'type', 'pearson'), 0,2)) % correlation between *enhanced* input patterns and learning patterns shows the expected zones corresponding to the clusters (especially with zscore)
    % imagesc(zscore(corr((srs').^1, pl, 'type', 'pearson'), 0,2)) % same thing for response patterns with learning patterns
    % imagesc(zscore(corr((srs').^1, (pts').^4, 'type', 'pearson'), 0,2)) % however, the input (test) patterns do not really correlate with the corresponding response patterns (except along the diagonal)...
    % imagesc(corr(sumsrs',sumpts.^4'));   % HOWEVER, if you sum over clusters AND enhance the input patterns, then each  summed-enhanced-input does correlated with the corresponding summed-response (alog the diagonal), although very weakly(r~=.1/.2, not significant)

    %print(gcf, '-dpng', '-r0', ['clusters_presentations_' FILENAME '.png'])   
end

if 1 == 0
    % Look at clusters between the firing patterns of different cells across
    % images, rather than the full response vectors to different images
    % across all cells.
    sr = rasterS;


    SUMRANGE = 0;
    for n=1:SUMRANGE
        sr = sr + circshift(rasterS, n);
    end
	if SUMRANGE > 0;  sr = sr(1:SUMRANGE:end,:); end
    
    nums=1:NBCELLS_S1;
    %sr = sr(:, sum(rasterS) > 1e-10);
    %nums = nums(sum(rasterS) > 1e-10);
    sr = sr(:, sum(rasterS) > 1);
    nums = nums(sum(rasterS) > 1);
    
    %sr = sr + 1e-5*rand(size(sr)); % To address the problem of images evoking no response at all - but then we can't use Spearman...
    %sr=sr(sum(sr,2)>1e-10,:);

    %sr=sr(sum(sr,2)>0, :);
    %sr(sr>0) = 1;

    %sr=sr(:,mean(sr)>0.0070);

    corrsr = corr(sr);  % Columns are cells
    % Clustering cells according to the similarity of their "firing correlation with other cells" pattern 
    %lnk = linkage(double(corrsr), 'single', 'correlation');
    %figure; [h t o] = dendrogram(lnk,size(corrsr,1));
    %% Clustering cells according to the similarity of their firing pattern 
    lnk = linkage(double(corrsr'), 'single', 'correlation');
    figure; [h t o] = dendrogram(lnk,size(corrsr',1));
    sorted = corrsr(o',o');
	sorted = sorted - diag(diag(sorted)); % We don't care that cells correlate with themselves!
    figure; imagesc(sorted)
    %print(gcf, '-dpng', '-r0', ['clusters_neurons_' FILENAME '.png'])   
end

%figure; imagesc(rasterS(3001:3501,:)');
%print(gcf, '-dpng', '-r0', ['rasterS_' FILENAME '.png'])   
%showS1;
%print(gcf, '-dpng', '-r0', ['RF_S1_' FILENAME '.png'])   
    

