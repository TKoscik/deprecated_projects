function outputImg = imscramble(inputimg)

imgGray = mat2gray(double(inputimg));

% read and rescale (0-1) image

imgSize = size(imgGray);

% generate random phase structure
randPhase = angle(fft2(rand(imgSize(1), imgSize(2))));


for layer = 1:imgSize(3)
    imgFourier(:,:,layer) = fft2(imgGray(:,:,layer));       
%Fast-Fourier transform
    Amp(:,:,layer) = abs(imgFourier(:,:,layer));       
%amplitude spectrum
    Phase(:,:,layer) = angle(imgFourier(:,:,layer));   
%phase spectrum
    Phase(:,:,layer) = Phase(:,:,layer) + randPhase;
%add random phase to original phase
    imgScrambled(:,:,layer) = ifft2(Amp(:,:,layer).*exp(sqrt(-1)*(Phase(:,:,layer))));   
%combine Amp and Phase then perform inverse Fourier
end

outputImg = real(imgScrambled); %get rid of imaginery part in image (due to rounding error)

