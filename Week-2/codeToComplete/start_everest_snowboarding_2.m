clearvars;
dst = double(imread('Everest.jpg'));
src = double(imread('Marti_Snow.png')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
mask_src=logical(imread('Marti_Snow_Mask_2.png'));
mask_dst=logical(imread('Everest_Mask_2.png'));

for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    drivingGrad_i = sol_DiFwd( src(:,:,nC), param.hi ) - sol_DiBwd( src(:,:,nC), param.hi );
    drivingGrad_j = sol_DjFwd( src(:,:,nC), param.hi ) - sol_DjBwd( src(:,:,nC), param.hi );

    driving_on_src = drivingGrad_i + drivingGrad_j;
    
    driving_on_dst = zeros(size(dst(:,:,1)));   
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst(:,:,nC), mask_dst,  param);
end



imshow(dst1/256)