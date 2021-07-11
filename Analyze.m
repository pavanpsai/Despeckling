clear all;clc
% Get Image File from the user
[FileName,PathName] = uigetfile(...
                            {'*.jpg;*.tif;*.png;*.gif','All Image Files';...
                            '*.*','All Files'},...
                            'Select Images','MultiSelect','off');

% Constructing FileName and FilePath for reading selected image
I = strcat(PathName,FileName);
RGB = imread(I);        % Read Selected Image
figure
imshow(RGB);title('Original Image');

% RGB = imread('football.jpg');
OI = preprocess(RGB);   % Preprocess Seleted Image
figure
imshow(OI);title('Preprocessed Image');
% Get variance of noise from user
v = input('Enter variance of speckle noise = ');
NI = AddSpecNoise(OI,v);
figure
imshow(NI);title('Noisy Image');
% Applying Savitzky-Golay Filter on Noisy Image
B = sgolayfilt(NI,3,41,[],2);
% Applying Median Filter on Noisy Image
C = medfilt2(NI,[7 7]);
% Get level of wavelet decomposition from user
L = input('Enter level of wavelet decomposition = ');
% Compute Non-Decimated Two Dimensional Wavelet Transform
AI = ndwt2(OI,L,'db1');
BI = ndwt2(B,L,'db1');
CI = ndwt2(C,L,'db1');
% Applying Bilateral Filtering 
[threshtemp MSEtemp PSNRtemp] = bft(NI,AI,BI,CI,L,2,'try');
% Selecting best threshold value from previous BFT ouput which gives
% maximum PSNR as selecting for minimum MSE degrades the visual quality of
% image.
thresh = threshtemp(PSNRtemp==max(max(PSNRtemp)));
thresh = max(max(thresh));
% Applying Bilateral Filter Algorithm for computing best result
[thresh MSE PSNR DI] = bft(NI,AI,BI,CI,L,2,'execute',thresh);
figure
% Visualize Image
subplot(2,3,1);imshow(OI);title('Original Image');
subplot(2,3,2);imshow(NI);title('Speckled Image');
subplot(2,3,3);imshow(B);title('Savitzky-Golay Filetered Image');
subplot(2,3,4);imshow(C);title('Median Filetered Image');
subplot(2,3,5);imshow(DI);title(' Bilater filter De-Speckled Image');
xlabel(['PSNR = ',num2str(PSNR),' dB','  ','MSE = ',num2str(MSE)]);

figure;
subplot(1,2,1);imshow(NI);title('Speckled Image');
subplot(1,2,2);imshow(DI);title('Bilater filter De-Speckled Image');
xlabel(['PSNR = ',num2str(PSNR),' dB','  ','MSE = ',num2str(MSE)]);