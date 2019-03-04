clearvars;
dst = double(imread('Everest.jpg'));
src = double(imread('Marti_Snow.png')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
mask_src=logical(imread('Marti_Snow_Mask.png'));
mask_dst=logical(imread('Everest_Mask.png'));

for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    sourceGrad_i = sol_DiFwd( src(:,:,nC), param.hi );
    sourceGrad_j = sol_DjFwd( src(:,:,nC), param.hj );
    
    dstGrad_i = sol_DiFwd( dst(:,:,nC), param.hi );
    dstGrad_j = sol_DjFwd( dst(:,:,nC), param.hj );


    drivingGrad_i = zeros(size(dst(:,:,1)));
    drivingGrad_i(mask_dst(:)) = sourceGrad_i(mask_src(:));
    

    drivingGrad_j = zeros(size(dst(:,:,1)));
    drivingGrad_j(mask_dst(:)) = sourceGrad_j(mask_src(:));
    
    
    drivingDestGrad_i = zeros(size(dst(:,:,1)));
    drivingDestGrad_j = zeros(size(dst(:,:,1)));
    drivingDestGrad_i(mask_dst(:)) = dstGrad_i(mask_dst(:));
    drivingDestGrad_j(mask_dst(:)) = dstGrad_j(mask_dst(:));

    
    drivingGrad_i (abs(drivingDestGrad_i)>abs(drivingGrad_i)) = drivingDestGrad_i(abs(drivingDestGrad_i)>abs(drivingGrad_i));
    drivingGrad_j (abs(drivingDestGrad_j)>abs(drivingGrad_j)) = drivingDestGrad_j(abs(drivingDestGrad_j)>abs(drivingGrad_j));
    
    driving_on_dst = +sol_DiBwd( drivingGrad_i(:,:), param.hi ) + sol_DjBwd( drivingGrad_j(:,:), param.hj );
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end



imshow(dst1/256)