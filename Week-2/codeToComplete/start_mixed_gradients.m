clearvars;
dst = double(imread('lena.png'));
src = double(imread('girl.png')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%masks to exchange: Eyes
mask_src=logical(imread('mask_src_eyes.png'));
mask_dst=logical(imread('mask_dst_eyes.png'));

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

%Mouth
%masks to exchange: Mouth
mask_src=logical(imread('mask_src_mouth.png'));
mask_dst=logical(imread('mask_dst_mouth.png'));
for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    drivingGrad_i = sol_DiFwd( src(:,:,nC), param.hi ) - sol_DiBwd( src(:,:,nC), param.hi );
    drivingGrad_j = sol_DjFwd( src(:,:,nC), param.hi ) - sol_DjBwd( src(:,:,nC), param.hi );

    driving_on_src = drivingGrad_i + drivingGrad_j;
    
    drivingGrad_i_dst = sol_DiFwd( dst(:,:,nC), param.hi ) - sol_DiBwd( dst(:,:,nC), param.hi );
    drivingGrad_j_dst = sol_DjFwd( dst(:,:,nC), param.hi ) - sol_DjBwd( dst(:,:,nC), param.hi );
    
    driving_on_dst0 = drivingGrad_i_dst + drivingGrad_j_dst;
 
    driving_on_dst = zeros(size(src(:,:,1))); 
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    driving_on_dst(abs(driving_on_dst0)>abs(driving_on_dst)) = driving_on_dst0(abs(driving_on_dst0)>abs(driving_on_dst));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = G7_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst,  param);
end

imshow(dst1/256)