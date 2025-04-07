function locFr = locGLRT(im, nParticles)
        filtIm = imgaussfilt(double(im),1);
        filtIm = double(im);
        chi2 = 100; 
        FWHM_pix = 7;%anything between 5 and 10 seem to do the job
        delta = 6;
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

%       figure(i*(j-1)+i)
%       imagesc(imfuse(I_hat,Loc2));
%    
    
    %get localization
    loca = regionprops(Loc2,'Centroid');
    
    locFr = reshape([loca.Centroid],2,size(loca,1))';
    locFr = locFr+delta; %shift the localization to convert from conv image to normal image
    
    if size(locFr,1) < nParticles
       newIm = double(im);
        for i = 1:size(locFr,1)
            sigma_pix = FWHM_pix / sqrt(8*log(2));
            %   calculating a gaussian PSF located in the middle of the ROI
            psf_model.name = 'gaussian';
            psf_model.sigma_x = sigma_pix;
            psf_model.sigma_y = psf_model.sigma_x;
            [Xgrid, Ygrid] = meshgrid(-delta:1:delta,-delta:1:delta);
            [ G ] = emitter_psf( Xgrid, Ygrid, 0, 0, psf_model);
            
            G = G./max(G(:)).* newIm(locFr(i,2),locFr(i,1));
            
            newIm(locFr(i,2)-delta:locFr(i,2)+delta,locFr(i,1)-delta:locFr(i,1)+delta) = newIm(locFr(i,2)-delta:locFr(i,2)+delta,locFr(i,1)-delta:locFr(i,1)+delta) -G;
            
            
            
        end
        
        newLoc = goldProj.locGLRT(newIm,nParticles-size(locFr,1));
        
        locFr = [locFr; newLoc(:,2),newLoc(:,1)];
               
    elseif size(locFr,1) >nParticles
        disp('too many particles detected, removing the excess');
        %Repeat has many time as too many particles
        corr2do = size(locFr,1);
        for i = 1:corr2do-nParticles
            %get matrix distance to find points that are too close
            distM = pdist(locFr);
            distM = squareform(distM);
            % find minimum distance
            [minD] = min(nonzeros(distM));
            % get the index
            idx = find(distM==minD);
            
            %find the points
            [a,~] = ind2sub(size(distM),idx);
            
            loc2Test = locFr(a,:);
            
            %get the pixel in the image to find which is brightest
            intVal(1,1) = im(loc2Test(1,2),loc2Test(1,1));
            intVal(2,1) = im(loc2Test(2,2),loc2Test(2,1));
            %take the index of the minimum for deletion
            [~,id2delete] = min(intVal);
            
            idx2Delete = ismember(locFr,loc2Test(id2delete,:));
            idx2Delete = idx2Delete(:,1) & idx2Delete(:,2);
            locFr(idx2Delete,:) = [];
        end       
        
    end
    
    locFr = [locFr(:,2), locFr(:,1)];
    
    
%     figure 
%     imagesc(im);
%     hold on
%     scatter(locFr(:,1),locFr(:,2))

function [ psf ] = emitter_psf( Xgrid, Ygrid, xpos, ypos, model )
%EMITTER_PSF calculates the PSF using the desired model

if strcmp(model.name, 'gaussian')
    sigma_x = model.sigma_x;
    sigma_y = model.sigma_y;
    

    X1 = ((Xgrid - xpos).^2) / (2*(sigma_x.^2));
    X2 = ((Ygrid - ypos).^2) / (2*(sigma_y.^2));
    A  = 1/(2*pi*sigma_x*sigma_y);  
    psf = A * exp(-(X1+X2));

end
end

end 