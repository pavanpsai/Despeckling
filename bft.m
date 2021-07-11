function [thresh MSE PSNR EI] = bft(OI,A,B,C,level,n,type,thresh)
% This function bft implements Brute Force Thresholding Algorithm as
% specified in the paper entitled as Despeckling of SAR Image Using 
% Adaptive and Mean Filters, by Syed Musharaf Ali, Muhammad Younus Javed, 
% and Naveed Sarfraz Khattak
%
% INPUT:
%           OI = Original Image
%           A = Nondecimated 2-D wavelet transform of OI
%           B = Nondecimated 2-D wavelet transform of Filtered OI
%           C = Nondecimated 2-D wavelet transform of Filtered OI
%           level = Decomposition level
%           n = Number of spacing between minimum and maximum threshold
%              (Larger n results time consuming)
% OUTPUT:
%           thresh = matrix of dimenshion (3*level,n)containing threshold 
%                    value
%           MSE = Mean Square Error
%           PSNR = Peak Singal to Noise Ratio
%
% 
% Implemented by ASHISH MESHRAM
% meetashish85@gmail.com http://www.facebook.com/ashishmeet

% Checking Input Arguments
if nargout==4 && nargin<8
    error('Not enough input arguments');
end
if nargin<7||isempty(type)
    type = 'try';
end
if nargin<6,error('Not enough input arguments');end
if nargin<5,error('Not enough input arguments');end
if nargin<4,error('Not enough input arguments');end
if nargin<3,error('Not enough input arguments');end
if nargin<2,error('Not enough input arguments');end
if nargin<1,error('Not enough input arguments');end



% Implementation Starts Here

switch lower(type)
    case 'try'  % Case for finding Best Threshold
                
        if level == 1   % First level decomposition
            
            ALH1 = A.dec{2,1};BLH1 = B.dec{2,1};CLH1 = C.dec{2,1}; % Vertical Detail
            AHL1 = A.dec{3,1};BHL1 = B.dec{3,1};CHL1 = C.dec{3,1}; % Horizonatal Detail
            AHH1 = A.dec{4,1};BHH1 = B.dec{4,1};CHH1 = C.dec{4,1}; % Diagonal Detail
   
            % Find minimum and maximum value of coefficients of each detail sub band
            Aminmax = [min(min(ALH1)) max(max(ALH1));...
               min(min(AHL1)) max(max(AHL1));...
               min(min(AHH1)) max(max(AHH1))];
    
            % Preallocating
            thresh = zeros(3,n);
            MSE = zeros(3,n);
            PSNR = zeros(3,n);
    
            % Vertical Detail
            coeff = linspace(Aminmax(1,1),Aminmax(1,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                ALH1(ALH1<coeff(k)) = CLH1(ALH1<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                ALH1 = ddm(ALH1,BLH1,coeff(k),'LH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                ALH1 = dsf(ALH1,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH1;
                A.dec{3,1} = AHL1;
                A.dec{4,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(1,k) PSNR(1,k)] = MetricsMeasurement(OI,EI);
                thresh(1,k) = coeff(k);
            end
    
            % Horizontal Detail
            coeff = linspace(Aminmax(2,1),Aminmax(2,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHL1(AHL1<coeff(k)) = CHL1(AHL1<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHL1 = ddm(AHL1,BHL1,coeff(k),'HL');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHL1 = dsf(AHL1,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH1;
                A.dec{3,1} = AHL1;
                A.dec{4,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(2,k) PSNR(2,k)] = MetricsMeasurement(OI,EI);
                thresh(2,k) = coeff(k);
            end
    
            % Diagonal Detail
            coeff = linspace(Aminmax(3,1),Aminmax(3,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHH1(AHH1<coeff(k)) = CHH1(AHH1<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHH1 = ddm(AHH1,BHH1,coeff(k),'HH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHH1 = dsf(AHH1,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH1;
                A.dec{3,1} = AHL1;
                A.dec{4,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(3,k) PSNR(3,k)] = MetricsMeasurement(OI,EI);
                thresh(3,k) = coeff(k);
            end
            
        elseif level == 2   % Second level decomposition
            ALH2 = A.dec{2,1};BLH2 = B.dec{2,1};CLH2 = C.dec{2,1}; % Vertical Detail
            AHL2 = A.dec{3,1};BHL2 = B.dec{3,1};CHL2 = C.dec{3,1}; % Horizonatal Detail
            AHH2 = A.dec{4,1};BHH2 = B.dec{4,1};CHH2 = C.dec{4,1}; % Diagonal Detail
            ALH1 = A.dec{5,1};BLH1 = B.dec{5,1};CLH1 = C.dec{5,1}; % Vertical Detail
            AHL1 = A.dec{6,1};BHL1 = B.dec{6,1};CHL1 = C.dec{6,1}; % Horizonatal Detail
            AHH1 = A.dec{7,1};BHH1 = B.dec{7,1};CHH1 = C.dec{7,1}; % Diagonal Detail
    
            % Find minimum and maximum value of coefficients of each detail sub band
            Aminmax = [min(min(ALH2)) max(max(ALH2));...
               min(min(AHL2)) max(max(AHL2));...
               min(min(AHH2)) max(max(AHH2));...
               min(min(ALH1)) max(max(ALH1));...
               min(min(AHL1)) max(max(AHL1));...
               min(min(AHH1)) max(max(AHH1))];
    
            % Preallocating
            thresh = zeros(6,n);
            MSE = zeros(6,n);
            PSNR = zeros(6,n);
    
            % Vertical Detail
            coeff = linspace(Aminmax(1,1),Aminmax(1,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                ALH2(ALH2<coeff(k)) = CLH2(ALH2<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                ALH2 = ddm(ALH2,BLH2,coeff(k),'LH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                ALH2 = dsf(ALH2,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH2;
                A.dec{3,1} = AHL2;
                A.dec{4,1} = AHH2;
                A.dec{5,1} = ALH1;
                A.dec{6,1} = AHL1;
                A.dec{7,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(1,k) PSNR(1,k)] = MetricsMeasurement(OI,EI);
                thresh(1,k) = coeff(k);
            end
    
            % Horizontal Detail
            coeff = linspace(Aminmax(2,1),Aminmax(2,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHL2(AHL2<coeff(k)) = CHL2(AHL2<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHL2 = ddm(AHL2,BHL2,coeff(k),'HL');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHL2 = dsf(AHL2,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH2;
                A.dec{3,1} = AHL2;
                A.dec{4,1} = AHH2;
                A.dec{5,1} = ALH1;
                A.dec{6,1} = AHL1;
                A.dec{7,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(2,k) PSNR(2,k)] = MetricsMeasurement(OI,EI);
                thresh(2,k) = coeff(k);
            end
    
            % Diagonal Detail
            coeff = linspace(Aminmax(3,1),Aminmax(3,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHH2(AHH2<coeff(k)) = CHH2(AHH2<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHH2 = ddm(AHH2,BHH2,coeff(k),'HH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHH2 = dsf(AHH2,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH2;
                A.dec{3,1} = AHL2;
                A.dec{4,1} = AHH2;
                A.dec{5,1} = ALH1;
                A.dec{6,1} = AHL1;
                A.dec{7,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(3,k) PSNR(3,k)] = MetricsMeasurement(OI,EI);
                thresh(3,k) = coeff(k);
            end
    
            % Vertical Detail
            coeff = linspace(Aminmax(4,1),Aminmax(4,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                ALH1(ALH1<coeff(k)) = CLH1(ALH1<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                ALH1 = ddm(ALH1,BLH1,coeff(k),'LH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                ALH1 = dsf(ALH1,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH2;
                A.dec{3,1} = AHL2;
                A.dec{4,1} = AHH2;
                A.dec{5,1} = ALH1;
                A.dec{6,1} = AHL1;
                A.dec{7,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(4,k) PSNR(4,k)] = MetricsMeasurement(OI,EI);
                thresh(4,k) = coeff(k);
            end
    
            % Horizontal Detail
            coeff = linspace(Aminmax(5,1),Aminmax(5,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHL1(AHL1<coeff(k)) = CHL1(AHL1<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHL1 = ddm(AHL1,BHL1,coeff(k),'HL');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHL1 = dsf(AHL1,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH2;
                A.dec{3,1} = AHL2;
                A.dec{4,1} = AHH2;
                A.dec{5,1} = ALH1;
                A.dec{6,1} = AHL1;
                A.dec{7,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(5,k) PSNR(5,k)] = MetricsMeasurement(OI,EI);
                thresh(5,k) = coeff(k);
            end
    
            % Diagonal Detail
            coeff = linspace(Aminmax(6,1),Aminmax(6,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHH1(AHH1<coeff(k)) = CHH1(AHH1<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHH1 = ddm(AHH1,BHH1,coeff(k),'HH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHH1 = dsf(AHH1,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH2;
                A.dec{3,1} = AHL2;
                A.dec{4,1} = AHH2;
                A.dec{5,1} = ALH1;
                A.dec{6,1} = AHL1;
                A.dec{7,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(6,k) PSNR(6,k)] = MetricsMeasurement(OI,EI);
                thresh(6,k) = coeff(k);
            end
    
        elseif level == 3   % Third level decomposition
    
            ALH3 = A.dec{2,1};BLH3 = B.dec{2,1};CLH3 = C.dec{2,1};    % Vertical Detail
            AHL3 = A.dec{3,1};BHL3 = B.dec{3,1};CHL3 = C.dec{3,1};    % Horizonatal Detail
            AHH3 = A.dec{4,1};BHH3 = B.dec{4,1};CHH3 = C.dec{4,1};    % Diagonal Detail
            ALH2 = A.dec{5,1};BLH2 = B.dec{5,1};CLH2 = C.dec{5,1};    % Vertical Detail
            AHL2 = A.dec{6,1};BHL2 = B.dec{6,1};CHL2 = C.dec{6,1};    % Horizonatal Detail
            AHH2 = A.dec{7,1};BHH2 = B.dec{7,1};CHH2 = C.dec{7,1};    % Diagonal Detail
            ALH1 = A.dec{8,1};BLH1 = B.dec{8,1};CLH1 = C.dec{8,1};    % Vertical Detail
            AHL1 = A.dec{9,1};BHL1 = B.dec{9,1};CHL1 = C.dec{9,1};    % Horizonatal Detail
            AHH1 = A.dec{10,1};BHH1 = B.dec{10,1};CHH1 = C.dec{10,1}; % Diagonal Detail
    
            % Find minimum and maximum value of coefficients of each detail sub band
            Aminmax = [min(min(ALH3)) max(max(ALH3));...
               min(min(AHL3)) max(max(AHL3));...
               min(min(AHH3)) max(max(AHH3));...
               min(min(ALH2)) max(max(ALH2));...
               min(min(AHL2)) max(max(AHL2));...
               min(min(AHH2)) max(max(AHH2));...
               min(min(ALH1)) max(max(ALH1));...
               min(min(AHL1)) max(max(AHL1));...
               min(min(AHH1)) max(max(AHH1))];
    
            % Preallocating
            thresh = zeros(9,n);
            MSE = zeros(9,n);
            PSNR = zeros(9,n);
            % Vertical Detail
            coeff = linspace(Aminmax(1,1),Aminmax(1,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                ALH3(ALH3<coeff(k)) = CLH3(ALH3<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                ALH3 = ddm(ALH3,BLH3,coeff(k),'LH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                ALH3 = dsf(ALH3,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH3;
                A.dec{3,1} = AHL3;
                A.dec{4,1} = AHH3;
                A.dec{5,1} = ALH2;
                A.dec{6,1} = AHL2;
                A.dec{7,1} = AHH2;
                A.dec{8,1} = ALH1;
                A.dec{9,1} = AHL1;
                A.dec{10,1} = AHH1;
        
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(1,k) PSNR(1,k)] = MetricsMeasurement(OI,EI);
                thresh(1,k) = coeff(k);
            end
    
            % Horizontal Detail
            coeff = linspace(Aminmax(2,1),Aminmax(2,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHL3(AHL3<coeff(k)) = CHL3(AHL3<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHL3 = ddm(AHL3,BHL3,coeff(k),'HL');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHL3 = dsf(AHL3,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH3;
                A.dec{3,1} = AHL3;
                A.dec{4,1} = AHH3;
                A.dec{5,1} = ALH2;
                A.dec{6,1} = AHL2;
                A.dec{7,1} = AHH2;
                A.dec{8,1} = ALH1;
                A.dec{9,1} = AHL1;
                A.dec{10,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(2,k) PSNR(2,k)] = MetricsMeasurement(OI,EI);
                thresh(2,k) = coeff(k);
            end
    
            % Diagonal Detail
            coeff = linspace(Aminmax(3,1),Aminmax(3,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHH3(AHH3<coeff(k)) = CHH3(AHH3<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHH3 = ddm(AHH3,BHH3,coeff(k),'HH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHH3 = dsf(AHH3,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH3;
                A.dec{3,1} = AHL3;
                A.dec{4,1} = AHH3;
                A.dec{5,1} = ALH2;
                A.dec{6,1} = AHL2;
                A.dec{7,1} = AHH2;
                A.dec{8,1} = ALH1;
                A.dec{9,1} = AHL1;
                A.dec{10,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(3,k) PSNR(3,k)] = MetricsMeasurement(OI,EI);
                thresh(3,k) = coeff(k);
            end
    
            % Vertical Detail
            coeff = linspace(Aminmax(4,1),Aminmax(4,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                ALH2(ALH2<coeff(k)) = CLH2(ALH2<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                ALH2 = ddm(ALH2,BLH2,coeff(k),'LH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                ALH2 = dsf(ALH2,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH3;
                A.dec{3,1} = AHL3;
                A.dec{4,1} = AHH3;
                A.dec{5,1} = ALH2;
                A.dec{6,1} = AHL2;
                A.dec{7,1} = AHH2;
                A.dec{8,1} = ALH1;
                A.dec{9,1} = AHL1;
                A.dec{10,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(4,k) PSNR(4,k)] = MetricsMeasurement(OI,EI);
                thresh(4,k) = coeff(k);
            end
    
            % Horizontal Detail
            coeff = linspace(Aminmax(5,1),Aminmax(5,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHL2(AHL2<coeff(k)) = CHL2(AHL2<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHL2 = ddm(AHL2,BHL2,coeff(k),'HL');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHL2 = dsf(AHL2,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH3;
                A.dec{3,1} = AHL3;
                A.dec{4,1} = AHH3;
                A.dec{5,1} = ALH2;
                A.dec{6,1} = AHL2;
                A.dec{7,1} = AHH2;
                A.dec{8,1} = ALH1;
                A.dec{9,1} = AHL1;
                A.dec{10,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(5,k) PSNR(5,k)] = MetricsMeasurement(OI,EI);
                thresh(5,k) = coeff(k);
            end
    
            % Diagonal Detail
            coeff = linspace(Aminmax(6,1),Aminmax(6,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHH2(AHH2<coeff(k)) = CHH2(AHH2<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHH2 = ddm(AHH2,BHH2,coeff(k),'HH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHH2 = dsf(AHH2,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH3;
                A.dec{3,1} = AHL3;
                A.dec{4,1} = AHH3;
                A.dec{5,1} = ALH2;
                A.dec{6,1} = AHL2;
                A.dec{7,1} = AHH2;
                A.dec{8,1} = ALH1;
                A.dec{9,1} = AHL1;
                A.dec{10,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(6,k) PSNR(6,k)] = MetricsMeasurement(OI,EI);
                thresh(6,k) = coeff(k);
            end
    
            % Vertical Detail
            coeff = linspace(Aminmax(7,1),Aminmax(7,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                ALH1(ALH1<coeff(k)) = CLH1(ALH1<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                ALH1 = ddm(ALH1,BLH1,coeff(k),'LH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                ALH1 = dsf(ALH1,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH3;
                A.dec{3,1} = AHL3;
                A.dec{4,1} = AHH3;
                A.dec{5,1} = ALH2;
                A.dec{6,1} = AHL2;
                A.dec{7,1} = AHH2;
                A.dec{8,1} = ALH1;
                A.dec{9,1} = AHL1;
                A.dec{10,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(7,k) PSNR(7,k)] = MetricsMeasurement(OI,EI);
                thresh(7,k) = coeff(k);
            end
    
            % Horizontal Detail
            coeff = linspace(Aminmax(8,1),Aminmax(8,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHL1(AHL1<coeff(k)) = CHL1(AHL1<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Dependent Mask');
                AHL1 = ddm(AHL1,BHL1,coeff(k),'HL');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHL1 = dsf(AHL1,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH3;
                A.dec{3,1} = AHL3;
                A.dec{4,1} = AHH3;
                A.dec{5,1} = ALH2;
                A.dec{6,1} = AHL2;
                A.dec{7,1} = AHH2;
                A.dec{8,1} = ALH1;
                A.dec{9,1} = AHL1;
                A.dec{10,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(8,k) PSNR(8,k)] = MetricsMeasurement(OI,EI);
                thresh(8,k) = coeff(k);
            end
    
            % Diagonal Detail
            coeff = linspace(Aminmax(9,1),Aminmax(9,2),n);
            for k = 1:length(coeff)
                % Replacing Low Level pixel of Image A with Image B
                disp('Replacing Low Level pixel of Image A with Image B');
                AHH1(AHH1<coeff(k)) = CHH1(AHH1<coeff(k));
                % Applying Direction Dependent Mask
                disp('Applying Direction Mask');
                AHH1 = ddm(AHH1,BHH1,coeff(k),'HH');
                % Applying Direction Smoothing Filter
                disp('Applying Direction Smoothing Filter');
                AHH1 = dsf(AHH1,3);
                % Replacing decomposition coefficients with new coefficients
                disp('Replacing decomposition coefficients with new coefficients');
                A.dec{2,1} = ALH3;
                A.dec{3,1} = AHL3;
                A.dec{4,1} = AHH3;
                A.dec{5,1} = ALH2;
                A.dec{6,1} = AHL2;
                A.dec{7,1} = AHH2;
                A.dec{8,1} = ALH1;
                A.dec{9,1} = AHL1;
                A.dec{10,1} = AHH1;
                % Inverse nondecimated 2-D wavelet transform
                disp('Inverse nondecimated 2-D wavelet transform');
                EI = indwt2(A,'aa',1);
                % Computing Results
                disp('Computing Results');
                [MSE(9,k) PSNR(9,k)] = MetricsMeasurement(OI,EI);
                thresh(9,k) = coeff(k);
            end
    
        else
            error(['Unhandeled level = ', num2str(level)]);
        end

    case 'execute'  % Case for executing BFT by the best threshold

        if level == 1
            
            ALH1 = A.dec{2,1};BLH1 = B.dec{2,1};CLH1 = C.dec{2,1}; % Vertical Detail
            AHL1 = A.dec{3,1};BHL1 = B.dec{3,1};CHL1 = C.dec{3,1}; % Horizonatal Detail
            AHH1 = A.dec{4,1};BHH1 = B.dec{4,1};CHH1 = C.dec{4,1}; % Diagonal Detail

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image B');
            ALH1(ALH1<thresh) = CLH1(ALH1<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            ALH1 = ddm(ALH1,BLH1,thresh,'LH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            ALH1 = dsf(ALH1,3);


            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image B');
            AHL1(AHL1<thresh) = CHL1(AHL1<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHL1 = ddm(AHL1,BHL1,thresh,'HL');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHL1 = dsf(AHL1,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image B');
            AHH1(AHH1<thresh) = CHH1(AHH1<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHH1 = ddm(AHH1,BHH1,thresh,'HH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHH1 = dsf(AHH1,3);
            
            
            % Replacing decomposition coefficients with new coefficients
            disp('Replacing decomposition coefficients with new coefficients');
            A.dec{2,1} = ALH1;
            A.dec{3,1} = AHL1;
            A.dec{4,1} = AHH1;
            % Inverse nondecimated 2-D wavelet transform
            disp('Inverse nondecimated 2-D wavelet transform');
            EI = indwt2(A,'aa',1);
            % Computing Results
            disp('Computing Results');
            [MSE PSNR] = MetricsMeasurement(OI,EI);
            
        elseif level == 2
            
            ALH2 = A.dec{2,1};BLH2 = B.dec{2,1};CLH2 = C.dec{2,1}; % Vertical Detail
            AHL2 = A.dec{3,1};BHL2 = B.dec{3,1};CHL2 = C.dec{3,1}; % Horizonatal Detail
            AHH2 = A.dec{4,1};BHH2 = B.dec{4,1};CHH2 = C.dec{4,1}; % Diagonal Detail
            ALH1 = A.dec{5,1};BLH1 = B.dec{5,1};CLH1 = C.dec{5,1}; % Vertical Detail
            AHL1 = A.dec{6,1};BHL1 = B.dec{6,1};CHL1 = C.dec{6,1}; % Horizonatal Detail
            AHH1 = A.dec{7,1};BHH1 = B.dec{7,1};CHH1 = C.dec{7,1}; % Diagonal Detail

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image B');
            ALH2(ALH2<thresh) = CLH2(ALH2<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            ALH2 = ddm(ALH2,BLH2,thresh,'LH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            ALH2 = dsf(ALH2,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image B');
            AHL2(AHL2<thresh) = CHL2(AHL2<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHL2 = ddm(AHL2,BHL2,thresh,'HL');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHL2 = dsf(AHL2,3);
            
            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image B');
            AHH2(AHH2<thresh) = CHH2(AHH2<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHH2 = ddm(AHH2,BHH2,thresh,'HH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHH2 = dsf(AHH2,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image B');
            ALH1(ALH1<thresh) = CLH1(ALH1<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            ALH1 = ddm(ALH1,BLH1,thresh,'LH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            ALH1 = dsf(ALH1,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image B');
            AHL1(AHL1<thresh) = CHL1(AHL1<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHL1 = ddm(AHL1,BHL1,thresh,'HL');
            %  Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHL1 = dsf(AHL1,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image B');
            AHH1(AHH1<thresh) = CHH1(AHH1<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHH1 = ddm(AHH1,BHH1,thresh,'HH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHH1 = dsf(AHH1,3);
            
            
            % Replacing decomposition coefficients with new coefficients
            disp('Replacing decomposition coefficients with new coefficients');
            
            A.dec{2,1} = ALH2;
            A.dec{3,1} = AHL2;
            A.dec{4,1} = AHH2;
            A.dec{5,1} = ALH1;
            A.dec{6,1} = AHL1;
            A.dec{7,1} = AHH1;
            % Inverse nondecimated 2-D wavelet transform
            disp('Inverse nondecimated 2-D wavelet transform');
            EI = indwt2(A,'aa',1);
            % Computing Results
            disp('Computing Results');
            [MSE PSNR] = MetricsMeasurement(OI,EI);
            
        elseif level == 3
            
            ALH3 = A.dec{2,1};BLH3 = B.dec{2,1};CLH3 = C.dec{2,1};    % Vertical Detail
            AHL3 = A.dec{3,1};BHL3 = B.dec{3,1};CHL3 = C.dec{3,1};    % Horizonatal Detail
            AHH3 = A.dec{4,1};BHH3 = B.dec{4,1};CHH3 = C.dec{4,1};    % Diagonal Detail
            ALH2 = A.dec{5,1};BLH2 = B.dec{5,1};CLH2 = C.dec{5,1};    % Vertical Detail
            AHL2 = A.dec{6,1};BHL2 = B.dec{6,1};CHL2 = C.dec{6,1};    % Horizonatal Detail
            AHH2 = A.dec{7,1};BHH2 = B.dec{7,1};CHH2 = C.dec{7,1};    % Diagonal Detail
            ALH1 = A.dec{8,1};BLH1 = B.dec{8,1};CLH1 = C.dec{8,1};    % Vertical Detail
            AHL1 = A.dec{9,1};BHL1 = B.dec{9,1};CHL1 = C.dec{9,1};    % Horizonatal Detail
            AHH1 = A.dec{10,1};BHH1 = B.dec{10,1};CHH1 = C.dec{10,1}; % Diagonal Detail

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image C');
            ALH3(ALH3<thresh) = CLH3(ALH3<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            ALH3 = ddm(ALH3,BLH3,thresh,'LH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            ALH3 = dsf(ALH3,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image C');
            AHL3(AHL3<thresh) = CHL3(AHL3<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHL3 = ddm(AHL3,BHL3,thresh,'HL');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHL3 = dsf(AHL3,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image C');
            AHH3(AHH3<thresh) = CHH3(AHH3<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHH3 = ddm(AHH3,BHH3,thresh,'HH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHH3 = dsf(AHH3,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image C');
            ALH2(ALH2<thresh) = CLH2(ALH2<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            ALH2 = ddm(ALH2,BLH2,thresh,'LH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            ALH2 = dsf(ALH2,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image C');
            AHL2(AHL2<thresh) = CHL2(AHL2<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHL2 = ddm(AHL2,BHL2,thresh,'HL');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHL2 = dsf(AHL2,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image C');
            AHH2(AHH2<thresh) = CHH2(AHH2<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHH2 = ddm(AHH2,BHH2,thresh,'HH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHH2 = dsf(AHH2,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image C');
            ALH1(ALH1<thresh) = CLH1(ALH1<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            ALH1 = ddm(ALH1,BLH1,thresh,'LH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            ALH1 = dsf(ALH1,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image C');
            AHL1(AHL1<thresh) = CHL1(AHL1<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Dependent Mask');
            AHL1 = ddm(AHL1,BHL1,thresh,'HL');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHL1 = dsf(AHL1,3);

            % Replacing Low Level pixel of Image A with Image B
            disp('Replacing Low Level pixel of Image A with Image C');
            AHH1(AHH1<thresh) = CHH1(AHH1<thresh);
            % Applying Direction Dependent Mask
            disp('Applying Direction Mask');
            AHH1 = ddm(AHH1,BHH1,thresh,'HH');
            % Applying Direction Smoothing Filter
            disp('Applying Direction Smoothing Filter');
            AHH1 = dsf(AHH1,3);
            
            % Replacing decomposition coefficients with new coefficients
            disp('Replacing decomposition coefficients with new coefficients');
            
            A.dec{2,1} = ALH3;
            A.dec{3,1} = AHL3;
            A.dec{4,1} = AHH3;
            A.dec{5,1} = ALH2;
            A.dec{6,1} = AHL2;
            A.dec{7,1} = AHH2;
            A.dec{8,1} = ALH1;
            A.dec{9,1} = AHL1;
            A.dec{10,1} = AHH1;
            % Inverse nondecimated 2-D wavelet transform
            disp('Inverse nondecimated 2-D wavelet transform');
            EI = indwt2(A,'aa',1);
            % Computing Results
            disp('Computing Results');
            [MSE PSNR] = MetricsMeasurement(OI,EI);
            
        else
            error(['Unhandeled level = ', num2str(level)]);
        end
        
    otherwise
        error(['Unhandeled type = ', type]);
end

