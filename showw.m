function showw(fname)
    if (nargin <1)
        fname = 'wff.txt';
    end
    wff = load(fname);
    BORDER=2;
%    NBCELLS_S1 = size(wff,1);
    NBCELLS_S1 = 100;
    RFSIZE = sqrt((size(wff,2)) / 2);
    sqs = ceil(sqrt(NBCELLS_S1));
    NBLI = sqs; NBCO = sqs;


    allrf = -100*ones((RFSIZE + BORDER) * NBLI, (RFSIZE + BORDER) * NBCO);


    %c1c = 3; [tmp, idx] = sort(w_C1(:,c1c), 'descend');
    idx = 1:NBCELLS_S1; %NBCELLS_S1:-1:1;
    %idx=unique([nums(o) 1:144], 'stable'); % To be used after drawing the clustering matrix *between cells* (not between presentations)
    n = 1;
    for li = 1:NBLI
        for co = 1:NBCO
            if n <= numel(idx)
                thisw = reshape(wff(idx(n), :), RFSIZE, 2*RFSIZE);
                thisrf = thisw(:,1:end/2) - thisw(:,end/2+1:end);
                %thisrf(:) = zscore(thisrf(:));
                %thisrf(:) = thisrf(:) - mean(thisrf(:));
                %thisrf = thisrf / max(abs(thisrf(:)));
                allrf(((li-1) * (RFSIZE + BORDER) + 1):((li-1) * (RFSIZE + BORDER) + 1)+RFSIZE-1, ((co-1) * (RFSIZE + BORDER) + 1):((co-1) * (RFSIZE + BORDER) + 1)+RFSIZE-1) = thisrf;
            end
            n = n + 1;
        end
    end

    allrf(allrf < -99) = .6;%allmin;
    allrf = allrf - min(allrf(:)); allrf = allrf ./ max(allrf(:));
    allrf = imresize(allrf,2);

    %figure;
    %imagesc(allrf);
    imshow(allrf);
    colormap(gray);
    axis off;
end
