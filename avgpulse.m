ll=load('lastnspikes_pulse.txt'); 
ll = ll(1:100,:);
sumf = sum(ll,2)';
clustercells = find(sumf>300)

allspikes = [];
allspikesNINL = [];

llNINL=load('lastnspikes_pulse_nolat_noinh.txt'); 
ll=load('lastnspikes_pulse.txt'); 


TIMEBIN = 1;
SMOOTHING = 1;

for STUDIEDCELL = clustercells
%No Time-binning, just smoothing
tt = reshape(ll(STUDIEDCELL, :), TIMEBIN, 350/TIMEBIN, 50); z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2); hold off;
if SMOOTHING 
    z = smooth(z, 10);
end
allspikes = [allspikes z];

%plot(1:TIMEBIN:350, smooth(z, 10), 'b', 'Linewidth', 2) 
tt = reshape(llNINL(STUDIEDCELL, :), TIMEBIN, 350/TIMEBIN, 50); z = mean(squeeze(1000*sum(tt,1)/TIMEBIN), 2); hold on;
%plot(1:TIMEBIN:350, smooth(z, 10), 'r', 'Linewidth', 2); axis([0  150 0 200]);
if SMOOTHING
    z = smooth(z, 10);
end
allspikesNINL = [allspikesNINL z];
end

avgNINL = mean(allspikesNINL);
% hold off; plot(mean(allspikes ./ (1.0 + repmat(avgNINL, 350, 1)), 2), 'r') ; hold on; plot(mean(allspikesNINL ./ (1.0 + repmat(avgNINL, 350, 1)), 2), 'b'); axis([0  150 0 10]); % Not much difference...
hold off; plot(1:TIMEBIN:350, mean(allspikes, 2), 'r', 'Linewidth', 2) ; hold on; plot(1:TIMEBIN:350, mean(allspikesNINL, 2), 'b', 'Linewidth', 2); axis([0  150 0 200]);

