clear
clc, %close all
nt=1900;
nx=320;
nz=320;
% nt=10000;
% nx=2301;
% nz=751;
% nt=1000;
% nx=737;
% nz=240;

npml=32;
Block_SizeX=16;
Block_SizeY=16;
nnx=2*npml+ceil(nx/Block_SizeX)*Block_SizeX
nnz=2*npml+ceil(nz/Block_SizeY)*Block_SizeY

% figure(1),clf
% fid=fopen('wav.dat','rb');
% for it=1:50:nt
%     fseek(fid, (it-1)*nnx*nnz*4,'bof');
%     x=fread(fid,[nnz nnx],'float');
%     xx=x(npml+1:nnz-npml,npml+1:nnx-npml);
%     imagesc(xx),colorbar
% %     colormap(seismic(3));colorbar
%     drawnow;
% %      pause;
% end
% fclose(fid);


figure(2),clf
fid=fopen('wav.dat','rb');     
x=fread(fid,[nnz nnx],'float');
    xx=x(npml+1:nnz-npml,npml+1:nnx-npml);
    imagesc(xx),
    %pcolor(xx),shading interp,axis ij,
    colormap(gray),colorbar
    
fclose(fid);




% ng=320;
% figure(3),clf
% fid=fopen('wav.dat','rb');
% x=fread(fid,[nt ng],'float');
% imagesc(x);
% fclose(fid);