function result = MSEgray(A,B,rand);
[mA nA,d] = size(A);
[mB nB,d] = size(B);
if (mA ~= mB | nA ~= nB)
   error('Arrays A and B must have the same size');
end;
if (d ~= 1)
   error('Only grayscale images!');
end;
A = double(A);
B = double(B);
result = 0;
% we laten de rand buiten beschouwing
for i = 1+rand:mA-rand
   for j = 1+rand:nA-rand
      result = result + (A(i,j,1) - B(i,j,1))^2;
   end;
end;
result = result / ((mA-2*rand)*(nA-2*rand));