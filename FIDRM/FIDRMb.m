
%**************************************************************************
% FRINRM: fuzzy random imulse noise reduction method
%
%
% Stefan Schulte (stefan.schulte@Ugent.be):
% Last modified: 15/01/06
%**************************************************************************


function Fs = FIDRMb(A,O)
  [M,N,DIM] = size(A);
   a = 100.0; b = 190.0;
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

   Windowsize = 1;
   F = double(A);
   F2 = double(A);

   begin = 0;
   loop = 1;
   while (begin == 0)
      for k = 0:4 
         F = double(F2);
         mem(M,N) = 0;
         
         mem= Detectionb(F,a/(2^k),b/(2^k));
         F2 = Denoise2(F,mem,Windowsize);
       %  MSEgray(F2,O(:,:,1),2)

         if (MSEgray(F,O,3) < MSEgray(F2,O,3)) | (loop == 40)
            begin = 1;
            break;
         end
      end
%      disp(sprintf('==> MSE %2d: %10.5f',loop,MSEgray(F2,O,3)));
      loop = loop +1;
   end
   
%   mem= Detectionb(A,a,b);
%   F = Denoise2(A,mem,Windowsize);
%   MSEgray(F,O(:,:,1),2)
   
%   mem2= Detectionb(F,a/2,b/2); % 30 en 60
%   F2 = Denoise2(F,mem2,Windowsize);
%   MSEgray(F2,O(:,:,1),2)
   
%   mem3= Detectionb(F2,a/4,b/4);
%   F3 = Denoise2(F2,mem3,Windowsize);
%   MSEgray(F3,O(:,:,1),2)
   
%   mem4= Detectionb(F3,a/8,b/8);
%   F4 = Denoise2(F3,mem4,Windowsize);
%   MSEgray(F4,O(:,:,1),2)
   
   
   Fs(:,:,1) = F;   
   Fs(:,:,2) = F;   
   Fs(:,:,3) = F;   
%   disp(sprintf('======================='));
%   disp(sprintf('Output MSE:  %10.5f',MSEgray(F,O(:,:,1),2)));
%   disp(sprintf('Output PSNR: %10.5f',(log(255^2/MSEgray(F,O(:,:,1),2))/log(10)*10)));
%   disp(sprintf('======================='));