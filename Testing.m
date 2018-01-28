clc;clear;close all;warning off all;
 
image_folder = 'BPN Tes';
filenames = dir(fullfile(image_folder, '*.png'));
total_images = numel(filenames);

for n = 1:total_images
    full_name = fullfile(image_folder, filenames(n).name);
    
    %preprocessing
    img = imread(full_name);
    img_gray=rgb2gray(img);
    
    cc = medfilt2(img_gray);
   
    %Tumor Segmentation
    T = 155;
    bw = im2bw(cc,T/255);

    %fs = get(0,'ScreenSize');
    %figure('Position',[0 0 fs(3)/2 fs(4)])

    SE = strel('disk',2);
    bw1 = imerode(bw,SE);

    %SE = strel('disk',0);
    bw2 = imdilate(bw1,SE);

    SE2 = strel('disk',4);
    bw3 = imerode(bw2,SE2);

    SE3 = strel('disk',4);

    %SE2 = strel('disk',1);
    bw4 = imdilate(bw3,SE3);

    bw6 = bwareaopen(bw4,350,8);
    %SE4 = strel('disk',2);
    %bw6 = imerode(bw5,SE4);

    bw6a = imfill(bw6, 'holes');

    SE5 = strel('disk',3);
    bw7 = imdilate(bw6a,SE5);

    if bw7 == 0
        brain1 = 255*uint8(bw7);
        brain_glcm = double(brain1);
        brain_glcm(brain_glcm==0) = NaN;
    else
        cc2 = cc;
        cc2(~bw7)=0;

        cc_resize = imresize(cc2, [256 256]);
        brain1 = cc_resize(cc_resize>0);
        brain_glcm = double(cc_resize);
        brain_glcm(brain_glcm==0) = NaN;
    end
    
    [N,M,L] = size(brain1);

    His2 =imhist(brain1)/(N*M); %Histogram Ternormalisasi

    Mean(n) = 0;
    for zi = 0:255
        Mean(n) = Mean(n) + (zi*His2(zi+1)); %Mean
    end

    Std2(n) = 0;
    for zi = 0:255
        Std2(n) = Std2(n) + ((zi-Mean(n)).^2*His2(zi+1));
    end
    Std2(n) = sqrt(Std2(n)); %Standard Deviation

    Ent(n) = 0;
    for zi = 0:255
        if His2(zi+1)>0
            Ent(n) = Ent(n) - (His2(zi+1)*log2(His2(zi+1))); %Entropy
        end
    end
    
    %Gray Level Coocurence Matrix (GLCM) Texture Feature Extraction
    jarak = 1; %distance beetwen pixel
    warning('off','Images:graycomatrix:scaledImageContainsNan');
    GLCM = graycomatrix(brain_glcm,'NumLevels',8, 'GrayLimits',[], 'Offset',[0 jarak; -jarak jarak; -jarak 0; -jarak -jarak]);
    stats = graycoprops(GLCM,{'contrast','homogeneity'});
    warning('on','Images:graycomatrix:scaledImageContainsNan');

    contrast = stats.Contrast;
    contrast0(n) = contrast(1);
    contrast45(n) = contrast(2);
    contrast90(n) = contrast(3);
    contrast135(n) = contrast(4);
    homogeneity = stats.Homogeneity;
    homogeneity0(n) = homogeneity(1);
    homogeneity45(n) = homogeneity(2);
    homogeneity90(n) = homogeneity(3);
    homogeneity135(n) = homogeneity(4);
end

input = [Mean;Std2;Ent;contrast0;contrast45;contrast90;contrast135;homogeneity0;homogeneity45;homogeneity90;homogeneity135]
t1 = [1;0;0;0]; %Glioma
t2 = [0;1;0;0]; %Meningioma
t3 = [0;0;1;0]; %Metastatic
t4 = [0;0;0;1]; %Normal
target = [t1 t1 t1 t1 t1 t1 t2 t2 t2 t3 t3 t3 t4 t4 t4]
target_ind = vec2ind(target);

a = 0;
b = 255;
ra = 0.9;
rb = 0.1;
pa = (((ra-rb) * (input - a)) / (b - a)) + rb

load jst

output = sim(net,pa)
testIndices = vec2ind(output) %predict

plotconfusion(target,output) %testing accuracy
output_test = testIndices.'
cp = classperf(target_ind, output_test);
get(cp)
cp.CountingMatrix
cp.CorrectRate * 100
cp.ErrorRate * 100
