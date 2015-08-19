% Example for matlab interface on Madagascar

% Add m8r api to matlab search path
path(path, '$RSFROOT/lib'); % replace $RSFROOT with your m8r find path

dims = rsf_dim('lena.rsf');

lena = zeros(dims');
rsf_read(lena, 'lena.rsf');

figure, imagesc(lena'), colormap('gray'), axis image, axis off

% Call RSF program
nlena  = m8r('sfnoise seed=2015 var=1400', lena);

figure, imagesc(nlena'), colormap('gray'), axis image, axis off

dnlena = wiener2(nlena,[5,5]);

figure, imagesc(dnlena'), colormap('gray'), axis image, axis off

% Creat RSF head from a exsit RSF file
rsf_create('nlena.rsf','lena.rsf');

% write data to RSF file
rsf_write(nlena, 'nlena.rsf');
        
dims2 = [512; 1024];
rsf_create('twolena.rsf', dims2)
rsf_write(nlena', 'twolena.rsf','same')
rsf_write(dnlena', 'twolena.rsf', 'same')





