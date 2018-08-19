clear;
clc;
addpath my_derain_code
addpath ToolBox

optt=[0.1  0.3 1 3 10 30 100 300 1000];




padsize=5;
%%
%%%---Load Rainy Video---%%%
load('waterfall_rainy_case4.mat');
load('waterfall_clean.mat');

%%%---Load Rainy Video---%%%
[O_clean, ~] =  rgb2gray_hsv(B_clean(:,:,:,1:100));
[tmpRain, O_hsv] = rgb2gray_hsv(Rainy);%rgb2hsv

for op=50
    %%% change it for your video
    opts.tol = 5*0.001;
    opts.beta = 50;
    opts.alpha1 = 100;
    opts.alpha2 = 10;
    opts.alpha3 = 100;
    opts.alpha4 = 100;
    for it=250
        opts.maxit=it;
        opts.tol=5*0.001;
        %--- Algorithm ---%
        tic
        O_Rainy=biger(tmpRain,padsize);
        [B_1 iter]=rain_removal_RGSYTL_ver_TV_2_1(O_Rainy,opts);
        B_1=smaller(B_1,padsize);
        toc
        PSNR= psnr(B_1,O_clean);
        SSIM=ssim(B_1,O_clean);
        [PSNR SSIM]
    end
end
B_1=gray2color_hsv(O_hsv,B_1);
implay(B_1);






