%**************************************************************************
% Fuzzy Impulse Noise Detection and Reduction Method
%
%  The paper of the FIDRM is proposed in: 
%
%  Stefan Schulte, Mike Nachtegael, Valérie De Witte, 
%  Dietrich Van Der Weken and  Etienne E. Kerre:
%  A Fuzzy Impulse Noise Detection and Reduction Method.
%  IEEE Transactions on Image Processing 15(5), 2006, 1153-1162
%  
% Stefan Schulte (stefan.schulte@Ugent.be):
% Last modified: 03/01/06
%
% Inputs:  A = the noisy input image
%          O = the original noise-free image
% Outputs:  Fs = the filtered image 
%**************************************************************************
function Fs = FIDRM(A,O)
   disp(sprintf(''));
   disp(sprintf('FIXED VALUE IMPULSE NOISE REDUCTION'));
   disp(sprintf(''));
   
   [M,N,DIM] = size(A);
   A = double(A);
   if (DIM ~= 1)
       A1 = A(:,:,1);
       clear A;
       A = A1;
       clear A1;
   end
   
   [M2,N2,DIM2] = size(O);
   if (DIM2 ~= 1)
       O1 = O(:,:,1);
       clear O;
       O = O1;
       clear O1;
   end
   clear O1;

   
   F = double(A);

   
   fixed = Detection(A);
   
   if (fixed(1,1) ~= -1)
      Windowsize = 1;
      F = Denoise(A,fixed,Windowsize);

      Fs(:,:,1) = F;   
      Fs(:,:,2) = F;   
      Fs(:,:,3) = F;   
   else
      Fs(:,:,1) = A;   
      Fs(:,:,2) = A;   
      Fs(:,:,3) = A;   
   end
   
   disp(sprintf('======================='));
   disp(sprintf('Output MSE:  %10.5f',MSEgray(F,O,2)));
   disp(sprintf('Output PSNR: %10.5f',(log(255^2/MSEgray(F,O,2))/log(10)*10)));
   disp(sprintf('======================='));

   disp(sprintf(''));
   disp(sprintf('RANDOM VALUE IMPULSE NOISE REDUCTION'));
   disp(sprintf(''));
   Fs = FIDRMb(F,O);
   disp(sprintf('======================='));
   disp(sprintf('Output MSE:  %10.5f',MSEgray(Fs(:,:,1),O,2)));
   disp(sprintf('Output PSNR: %10.5f',(log(255^2/MSEgray(Fs(:,:,1),O,2))/log(10)*10)));
   disp(sprintf('======================='));

   