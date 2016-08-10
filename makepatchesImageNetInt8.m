RFSIZE = 17;
CROPMARGIN = 13;
CROPSIZE = RFSIZE + 2*CROPMARGIN;
NBFRAMES = 110000; %110000;

MAXRATE = 126; % Should be < 127 for Char conversion

RAND = 0;


%DoG = fspecial('gaussian', 9, .5) - fspecial('gaussian', 9, 1.0);
DoG = fspecial('gaussian', 21, 1) - fspecial('gaussian', 21, 2.0);
patches = zeros(RFSIZE*RFSIZE* NBFRAMES, 1);


% =====>>>> Zscoring or dividing by sum(abs) create patches with strong "dots"!!!!

rng(0);

FNAMES = dir('ImageNet/images'); FNAMES=FNAMES(3:end); NBFILES = numel(FNAMES);

warning('');

disp('Reading images....');
images={};
for nn=1:NBFILES
    if FNAMES(nn).isdir
        continue;
    end
    im = imread(['ImageNet/images/' FNAMES(nn).name]);
    if length(lastwarn) == 0 % No warning?
        images{nn} = mean(im,3);
    else
        disp('Skipping image...');
        warning('');
    end
end
NBFILES = numel(images);
disp('Read all images!');

RF_SD = 6;
z=fspecial('gaussian', RFSIZE, RF_SD); z=z./max(z(:)); GAUSSIANMASK = z;

sump = zeros(RFSIZE);

patchesim={};
patchesdata = zeros(NBFRAMES, RFSIZE * RFSIZE);
pos = 1;

numframe=1;

while numframe <= NBFRAMES
    if mod(numframe, 1000) == 1
        disp(numframe);
    end
    numim = randi(NBFILES); 
    im = images{numim};
    XSIZE = size(im, 2);  YSIZE = size(im, 1);
    if (XSIZE <= CROPSIZE+2) || (YSIZE <= CROPSIZE+2) % In case image is too small...
        continue;
    end
    x = randi(XSIZE -1 - CROPSIZE); y = randi(YSIZE - 1 - CROPSIZE);
    mypatch_large = im(y:y+CROPSIZE-1, x:x+CROPSIZE-1);

    
    %mypatch_large = imrotate(mypatch_large, -18, 'bilinear',  'crop');
    
    zangle = floor(rand() * 360);
    mypatch_large = imrotate(mypatch_large, zangle, 'bilinear',  'crop');
    
    mypatch_large = imfilter(mypatch_large, DoG);
    
    mypatch = mypatch_large(CROPMARGIN + 1:end-CROPMARGIN, CROPMARGIN + 1:end-CROPMARGIN); 


    % The following may not be useful given the rotation (except maybe for the transpose?) but why not..
    if (rand() > .5)
        mypatch = mypatch';
    end
    if (rand() > .5)
        mypatch = flipud(mypatch);
    end
    if (rand() > .5)
        mypatch = fliplr(mypatch); 
    end
    %mypatch(:) = zscore(mypatch(:));


    if (RAND)
        mypatch = mypatch(randperm(numel(mypatch)));
    end


    mypatch(:) = mypatch(:)-mean(mypatch(:));
    

    if mean(mypatch(:) ~= mode(mypatch(:))) < .1
        disp(['Skipping image due to insufficient variation - img num ' num2str(numim) ', frame num ' num2str(numframe)]);
    else
        mypatch(:) = MAXRATE * mypatch(:)./(1e-12 + max(abs(mypatch(:))));
        patches(pos:pos + RFSIZE * RFSIZE -1, 1) = mypatch(:);
        patchesim{numframe}  =  mypatch;
        patchesdata(numframe, :) = mypatch(:)';
        pos = pos + RFSIZE * RFSIZE;
        numframe = numframe + 1;
    end
    
    %mypatch = mypatch .* GAUSSIANMASK;
    
    
    %mypatch(:) = MAXRATE * mypatch(:)./(1e-12 + max(abs(mypatch(:))));
    
    
end

%patches = patches-mean(patches); patches = MAXRATE * patches./max(abs(patches));


patches = int8(patches);

%FNAME = 'patchesCenteredUnscaledImageNetONOFF_rotated';
FNAME = ['patchesCenteredScaledBySumTo' num2str(MAXRATE) 'ImageNetONOFFRotatedNewInt8'];
if (RAND)
    FNAME = [FNAME '_random'];
end
save([FNAME '.mat'], 'patches');

%fid=fopen(['patchesCenteredScaledBySumTo' num2str(MAXRATE) 'EachFrameRotated_Gauss' num2str(RF_SD) '.dat'], 'w');
fid=fopen([FNAME '.bin.dat'], 'w');
fwrite(fid, patches, 'int8');
fclose(fid);

if (1 == 0)
START = RFSIZE * RFSIZE * 20;
pp = reshape(patches(1+START:START+ RFSIZE*RFSIZE*30), RFSIZE, RFSIZE*30);
imagesc(pp); colormap(gray); axis equal;
end
