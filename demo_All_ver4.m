%%% ---demo--- %%%
%%%, date: 18/10/2016

%clear all;close all;clc;
%%
%%%--- Load Video Name ---%%%
clear;
clc;
addpath data
addpath my_derain_code
addpath ToolBox
addpath synthetic

kkk=1;frames=150;big=1;


optt=[0.1  0.3 1 3 10 30 100 300 1000];
totnum=length(optt);



padsize=5;
for video_num=1:1
    %%
    %%%---Load Rainy Video---%%%
    videofilename=['waterfall_rainy_case4.mat'];
    load(videofilename);
    load('waterfall_clean.mat');

    %%%---Load Rainy Video---%%%
    [O_clean, O_hsv] =  rgb2gray_hsv(B_clean(:,:,:,1:100));
    [O_Rainy,~]=rgb2gray_hsv(Rainy);%rgb2hsv
    tmpRain=O_Rainy;
    for op=50
        opts.tol=5*0.001;
        n=1;
        name=['.\opt\' num2str(video_num) '.txt'];
        fid=fopen(name,'r');
        opts.beta = 50;
        for j=1:4
            eval(['opts.alpha' num2str(j) '=fscanf(fid,''%f'',1);']);
        end
        fclose(fid);
        bestpsnr=0;
        bestopt=opts;
        bestopt.beta=50;
        lastopt=bestopt;
        bestssim=0;
        for it=250
            opts.maxit=it;
            opts.tol=5*0.001;
            %--- Algorithm ---%
            tic
            O_Rainy=biger(tmpRain,padsize);
            [B_1 iter]=rain_removal_RGSYTL_ver_TV_2_1(O_Rainy,opts);
            B_1=smaller(B_1,padsize);
            PSNR1= psnr(B_1,O_clean);
            SSIM1=ssim(B_1,O_clean);
            [PSNR1 SSIM1]
            if PSNR1>bestpsnr
                bestpsnr=PSNR1;
                bestopt=opts;
                bestssim=SSIM1;
                bestopt.maxit=iter;
                BB=B_1;
            end
            toc
        end
              
    end
end







