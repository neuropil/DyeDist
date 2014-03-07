clear

cd('G:\DropBox\Dropbox\JS_ImageAnalysis\GUI_settings_beta')

image = imread('Cfos_CGRP_test_1.tif');
blueChannel = image(:,:,3);
blueChannel = imresize(blueChannel,0.1);
figure; imshow(blueChannel);




%%

i = 3

numY = numel(yCor);

minSorty = sort(yCor);

minNow = minSorty(1:i);

yVal = ismember(yCor, minNow)