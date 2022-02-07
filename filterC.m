%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  C-weighting Filter                  %
%              with MATLAB Implementation              %
%                                                      %
% Author: Ph.D. Eng. Hristo Zhivomirov        12/17/18 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xC = filterC(x, fs)

% function: xC = filterC(x, fs)
% x - original signal in the time domain
% fs - sampling frequency, Hz
% xA - filtered signal in the time domain

% Note: The C-weighting filter's coefficients 
% are acccording to IEC 61672-1:2002 standard 

% determine the signal length
xlen = length(x);

% number of unique points
NumUniquePts = ceil((xlen+1)/2);

% FFT
X = fft(x);

% fft is symmetric, throw away the second half
X = X(1:NumUniquePts);

% frequency vector with NumUniquePts points
f = (0:NumUniquePts-1)*fs/xlen;

% C-weighting filter coefficients
c1 = 12194.217^2;
c2 = 20.598997^2;

% evaluate the C-weighting filter in the frequency domain
f = f.^2;
num = c1*f;
den = (f+c2) .* (f+c1);
C = 1.0072*num./den;

% filtering in the frequency domain
XC = X(:).*C(:);

% reconstruct the whole spectrum
if rem(xlen, 2)                     % odd xlen excludes the Nyquist point
    XC = [XC; conj(XC(end:-1:2))];
else                                % even xlen includes the Nyquist point
    XC = [XC; conj(XC(end-1:-1:2))];
end

% IFFT
xC = real(ifft(XC));
    
% represent the filtered signal in the form of the original one
xC = reshape(xC, size(x));

end