clear all;
close all;
clc;

% im_name='3_12_s.bmp';
im_name='7_9_s.bmp';
% im_name='2_1_s.bmp';


% TODO: Update library path
% Add  library paths
basedir='C:\Users\FV4RSR0\Desktop\M2\exercisis\low_level_segmentation\low_level_segmentation\UGM_2011\UGM';
addpath(genpath(basedir));

convert2LAB = true;


%Set model parameters
%cluster color
K=3; % Number of color clusters (=number of states of hidden variables)

%Pair-wise parameters
smooth_term=[0.0 2]; % Potts Model

%Load images
im = double(imread(im_name));


NumFils = size(im,1);
NumCols = size(im,2);

%Convert to LAB colors space
% TODO: Uncomment if you want to work in the LAB space
%

if convert2LAB
    im = RGB2Lab(im);
end


%Preparing data for GMM fiting
%
% TODO: define the unary energy term: data_term
% nodePot = P( color at pixel 'x' | Cluster color 'c' )  

im=double(im);
x=reshape(im,[size(im,1)*size(im,2), size(im,3)]);
gmm_color = gmdistribution.fit(x,K);
mu_color=gmm_color.mu;

% Estimate Unary potentials
data_term=gmm_color.posterior(x);
[~,c] = max(data_term,[],2);

% Standardize Features
%Xstd = UGM_standardizeCols(reshape(im,[1 1 NumFils*NumCols*3]),1);
%Xstd = UGM_standardizeCols(reshape(data_term,[NumFils*NumCols,K]),1);

nodePot = zeros(NumFils*NumCols,K);
for i = 1:K
    %nodePot(:,1) = exp(-1-2.5*Xstd(:,i));
    %nodePot(:,i) = exp(-1-2.5*data_term(:,i));
    nodePot(:,i) = data_term(:,i);
end


%Building 4-grid
%Build UGM Model for 4-connected segmentation
disp('create UGM model');

% Create UGM data
[edgePot,edgeStruct] = CreateGridUGMModel(NumFils, NumCols, K ,smooth_term, data_term);


if ~isempty(edgePot)

    % color clustering
    [~,c] = max(reshape(data_term,[NumFils*NumCols K]),[],2);
    im_c= reshape(mu_color(c,:),size(im));
    
    % Call different UGM inference algorithms
    display('Loopy Belief Propagation'); tic;
    [nodeBelLBP,edgeBelLBP,logZLBP] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);toc;
    [im_aux,idx] = max(nodeBelLBP,[],2);
    im_lbp= reshape(mu_color(idx,:),size(im));
    
    % Max-sum
    display('Max-sum'); tic;
    decodeLBP = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
    im_bp= reshape(mu_color(decodeLBP,:),size(im));
    toc;
    
    
    % TODO: apply other inference algorithms and compare their performance
    %
    % - Graph Cut
    % - Linear Programing Relaxation
    
    figure
    if convert2LAB
        subplot(2,2,1),imshow(Lab2RGB(im));xlabel('Original');
        subplot(2,2,2),imshow(Lab2RGB(im_c),[]);xlabel('Clustering without GM');
        subplot(2,2,3),imshow(Lab2RGB(im_bp),[]);xlabel('Max-Sum');
        subplot(2,2,4),imshow(Lab2RGB(im_lbp),[]);xlabel('Loopy Belief Propagation');
    else
        subplot(2,2,1),imshow(uint8(im));xlabel('Original');
        subplot(2,2,2),imshow(uint8(im_c),[]);xlabel('Clustering without GM');
        subplot(2,2,3),imshow(uint8(im_bp),[]);xlabel('Max-Sum');
        subplot(2,2,4),imshow(uint8(im_lbp),[]);xlabel('Loopy Belief Propagation');
    end
else
   
    error('You have to implement the CreateGridUGMModel.m function');

end