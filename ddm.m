function A = ddm(A,B,thresh,detail)
% This function ddm implements Direction Dependent Mask as mentioned in
% paper entitled as Despeckling of SAR Image Using Adaptive and Mean 
% Filters, by Syed Musharaf Ali, Muhammad Younus Javed, and Naveed Sarfraz
% Khattak
%
% INPUT:
%       A = Two Dimensional Matrix
%       B = Two Dimensional Matrix
%       thersh = Threshold value
% OUTPUT:
%       A = Processed Two Dimensonal Matrix
% 
% Algorithm Details:
%
% Pixel Level Classification
% if (pixel<threshold) 
%   ----> Low Level Pixel
% else
%   ----> High Level Pixel
%
% If the dark portion of mask or window contains at least one high level 
% pixel, then keep the original value of pixel intact else replace the 
% value of pixel of A by the corresponding pixel value of B.
% For this an edge detection mask is used for each detail coefficients
% For visualization run the following code
% 
% maskVr = [1 0 1;1 1 1;1 0 1];
% maskHr = [1 1 1;0 1 0;1 1 1];
% maskdiag = [0 1 0;1 1 1;0 1 0];
% figure(1);
% subplot(1,3,1);imshow(maskVr,'InitialMagnification','fit');title('Vertical detail edge detection mask');
% subplot(1,3,2);imshow(maskHr,'InitialMagnification','fit');title('Horizontal detail edge detection mask');
% subplot(1,3,3);imshow(maskdiag,'InitialMagnification','fit');title('Diagonal detail edge detection mask');
%
% USAGE EXAMPLES:
%
% RGB = imread('football.jpg');
% OI = rgb2gray(RGB);
% NI = imnoise(OI,'speckle');
% % As the maximum value of pixel on OI is 225, so taking that value as
% % thresholding to apply Direction Dependent Mask
% C = ddm(NI,OI,255,'HH');  % ddm on Detail Coefficient
% subplot(1,2,1);imshow(NI);title('Speckled Image');
% subplot(1,2,2);imshow(C);title('Directional Dependent Mask result');

% Implemented by ASHISH MESHRAM
% meetashish85@gmail.com http://www.facebook.com/ashishmeet


% Checking Input Arguments
if nargin<4,error('Not enough input argument');end
if nargin<3,error('Not enough input argument');end
if nargin<2,error('Not enough input argument');end
if nargin<1
    error('Not enough input argument');
elseif size(A)~=size(B)
    error('Dimension of A and B should be same');
end

% Implementation starts here
[ra ca] = size(A);              % Get Dimension of matrix A
PA = zeros(ra + 2,ca + 2);      % Initializing Padding Matrix
PA(2:end - 1,2:end - 1) = A;    % Padding Image

switch detail
    case 'LH'   % Vertical Detail Coefficients
        % Scanning Original Image by traversing mask on it
        for p = 1:ra
            for q = 1:ca
                mask = PA(p:2 + p,q:2 + q);
                % Vertical detail edge detection mask
                % Execute below two lines of code for visualization
                % maskVr = [1 0 1;1 1 1;1 0 1];
                % Execute below two lines of code for visualization
                % imshow(maskVr,'InitalMagnification','fit');
                if mask(4)<thresh || mask(6)<thresh
                    A(p,q) = B(p,q);
                end
            end
        end

    case 'HL'   % Horizontal Detail Coefficients
        % Scanning Original Image by traversing mask on it
        for p = 1:ra
            for q = 1:ca
                mask = PA(p:2 + p,q:2 + q);
                % Horizontal detail edge detection mask
                % Execute below two lines of code for visualization
                % maskHr = [1 1 1;0 1 0;1 1 1];
                % imshow(maskHr,'InitalMagnification','fit');
                if mask(2)<thresh || mask(8)<thresh
                    A(p,q) = B(p,q);
                end

            end
        end
        
    case 'HH'   % Diagonal Detail Coefficients
        % Scanning Original Image by traversing mask on it
        for p = 1:ra
            for q = 1:ca
                mask = PA(p:2 + p,q:2 + q);
                % Diagonal detail edge detection mask
                % Execute below two lines of code for visualization
                % maskdiag = [0 1 0;1 1 1;0 1 0];
                % imshow(maskdiag,'InitalMagnification','fit');
                if mask(1)<thresh || mask(3)<thresh ||...
                   mask(7)<thresh || mask(9)<thresh
                    A(p,q) = B(p,q);
                end
            end
        end
        
    otherwise
        error(['Unhandeled detail = ', detail]);
        
end