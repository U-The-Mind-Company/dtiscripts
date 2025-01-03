function res = ifft3c(x)
fctr = size(x,1)*size(x,2)*size(x,3);
for n=1:size(x,4)
res(:,:,:,n) = sqrt(fctr)*fftshift(fftshift(fftshift(ifftn(ifftshift(ifftshift(ifftshift(x(:,:,:,n),1),2),3)),1),2),3);
end
