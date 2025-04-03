

function locFr = locGLRT(im, nParticles)
    filtIm = imgaussfilt(double(im),1);
    chi2 = 100; 
    FWHM_pix = 3;%anything between 5 and 10 seem to do the job
    delta = 5;
    [I_hat] = Localization.GLRTfiltering2(filtIm,delta,FWHM_pix,chi2);
    I_hat(I_hat<0)= 0;

    Loc2 = imregionalmax(I_hat);


    %delete localization in background
    [N,bins] = histcounts(I_hat(Loc2),20);
    % figure
    % bar(bins(1:end-1),N)
    threshold = bins(find(N==0,1,'first'));
    % 
    % 
     idx = find(Loc2);
     thresh = I_hat(Loc2)<threshold;
    % 
     Loc2(idx(thresh)) = 0;

    %  figure
    % imagesc(imfuse(I_hat,Loc2));
    %get localization
    loca = regionprops(Loc2,'Centroid');

    locFr = reshape([loca.Centroid],2,size(loca,1))';
    locFr = locFr+delta; %shift the localization to convert from conv image to normal image
    
    if size(locFr,1)>nParticles
       disp('test') 
        
    end
    
    locFr = [locFr(:,2), locFr(:,1)];
    
    
%     figure 
%     imagesc(im);
%     hold on
%     scatter(locFr(:,1),locFr(:,2))



end 