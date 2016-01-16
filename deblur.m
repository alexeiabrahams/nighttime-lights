clear
clc
tic
%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETERS:
input_image_pathname = 'E:\nighttime_lights\F152000_topcoded_addis_ababa.tif';%avg_vis image
input_image_perc_pathname = 'E:\nighttime_lights\F152000_topcoded_addis_ababa.tif';%frequency image
output_image_pathname = 'E:\nighttime_lights\F152000_topcoded_addis_ababa_DEBLURRED.tif';%name of the file to be produced by this script
pixel_width = .915;%user should measure the width of a pixel in the image (kilometers)
pixel_length = .921;%user should measure the length of a pixel in the image (kilometers)
uniform_cutoff=4;%user should set this based on inspection of image -- what seems to be the typical 'background' pixel value?
satellite='F15';%user should set this
year='2000';%user should set this
divide_input_image_perc_by = 1.0;%default value: 1.0. this is to make sure your frequency image is scaled between 0 and 100.
sigma_x2 = 1.65^2;%default value: 1.65^2
sigma_y2 = 1.65^2;%default value: 1.65^2
radiance_calibrated=0; %set to 1 if you are using radiance-calibrated imagery, 0 otherwise.
%%NOTHING BELOW THIS LINE NEEDS TO BE CHANGED BY THE USER :-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[input_image, input_georeference] = geotiffread(input_image_pathname);
input_image_perc = imread(input_image_perc_pathname);
input_image_perc = double(input_image_perc)/divide_input_image_perc_by;
num_rows = size(input_image,1);
num_cols = size(input_image',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wu_et_al_parameters = {'1', 'F10', '1992', '0.8959', '1.0310', '0.9492';
    '2', 'F10', '1993', '0.6821', '1.1181', '0.8731';
    '3', 'F10', '1994', '0.9127', '1.0640', '0.9112';
    '4', 'F12', '1994', '0.4225', '1.3025', '0.8559';
    '5', 'F12', '1995', '0.3413', '1.3604', '0.9275';
    '6', 'F12', '1996', '0.9247', '1.0576', '0.9541';
    '7', 'F12', '1997', '0.3912', '1.3182', '0.9042';
    '8', 'F12', '1998', '0.9734', '1.0312', '0.9125';
    '9', 'F12', '1999', '1.2743', '0.9539', '0.8846';
    '10', 'F14', '1997', '1.3041', '0.9986', '0.8945';
    '11', 'F14', '1998', '0.9824', '1.1070', '0.9589';
    '12', 'F14', '1999', '1.0347', '1.0904', '0.9479';
    '13', 'F14', '2000', '0.9885', '1.0702', '0.9047';
    '14', 'F14', '2001', '0.9282', '1.0928', '0.9706';
    '15', 'F14', '2002', '0.9748', '1.0857', '0.9752';
    '16', 'F14', '2003', '0.9144', '1.1062', '0.9156';
    '17', 'F15', '2000', '0.8028', '1.0855', '0.9242';
    '18', 'F15', '2001', '0.8678', '1.0646', '0.8700';
    '19', 'F15', '2002', '0.7706', '1.0920', '0.8854';
    '20', 'F15', '2003', '0.9852', '1.1141', '0.9544';
    '21', 'F15', '2004', '0.8640', '1.1671', '0.9352';
    '22', 'F15', '2005', '0.5918', '1.2894', '0.9322';
    '23', 'F15', '2006', '0.9926', '1.1226', '0.9145';
    '24', 'F15', '2007', '1.1823', '1.0850', '0.9041';
    '25', 'F16', '2004', '0.7638', '1.1507', '0.9123';
    '26', 'F16', '2005', '0.6984', '1.2292', '0.8620';
    '27', 'F16', '2006', '0.9028', '1.1306', '0.9412';
    '28', 'F16', '2007', '0.8864', '1.1112', '0.9576';
    '29', 'F16', '2008', '0.9971', '1.0977', '0.9653';
    '30', 'F16', '2009', '1.4637', '0.9858', '0.8735';
    '31', 'F18', '2010', '0.8114', '1.0849', '0.9542'}

noaa_radcal_intercalibrate = {'1', 'F12', '1997', '0.915', '4.336';
    '2', 'F12', '1999', '0.78', '1.423';
    '3', 'F15', '2000', '0.71', '3.658';
    '4', 'F15', '2003', '0.797', '3.736';
    '5', 'F14', '2004', '0.761', '1.062';
    '6', 'F16', '2006', '1.000', '0.0';
    '7', 'F16', '2010', '1.195', '2.196';
    '8', 'F16', '2011', '1.246', '-1.987'}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust standard deviations to appropriate
%pixel size:
sigma_x2 = sigma_x2/(pixel_width^2); %converts sigma_x2 from units of km to units of bins
sigma_y2 = sigma_y2/(pixel_length^2); %converts sigma_y2 from units of km to units of bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTERCALIBRATE IMAGERY USING WU ET AL 2012 FOR NON-RADIANCE-CALIBRATED (TOPCODED) IMAGERY, OR HSU ET AL (2015) FOR RADIANCE-CALIBRATED IMAGERY:
yes_calibrate=1;
if yes_calibrate==1
if radiance_calibrated==0
    for i = 1:31
        if strcmp(wu_et_al_parameters(i, 2), satellite)==1 & strcmp(wu_et_al_parameters(i, 3), year)==1
            a = wu_et_al_parameters(i, 4);
            b = wu_et_al_parameters(i, 5);
            a = str2num(a{1});
            b = str2num(b{1});
            input_image = (double(input_image + 1.0).^b)*a-1.0;
            uniform_cutoff = (double(uniform_cutoff + 1.0)^b)*a-1.0;
        end
    end
end
if radiance_calibrated==1
    for i=1:8
        if strcmp(noaa_radcal_intercalibrate(i, 2), satellite)==1 & strcmp(noaa_radcal_intercalibrate(i, 3), year)==1
            b = noaa_radcal_intercalibrate(i, 4);
            a = noaa_radcal_intercalibrate(i, 5);
            input_image = a + b*double(input_image);
            uniform_cutoff = a + b*uniform_cutoff;
        end
    end
end
end
%subtract background noise:
input_image = imadd(input_image,-uniform_cutoff);
%set negative values to zero:
input_image(input_image<0.0) = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Generate PSF
mu = [0 0];
Sigma2 = [sigma_x2  0; 0 sigma_y2];
x1 = -10:1:10; x2 = -10:1:10;
[X1,X2] = meshgrid(x1,x2);
Gaussian_PSF = mvnpdf([X1(:) X2(:)],mu,Sigma2);
Gaussian_PSF = reshape(Gaussian_PSF,length(x2),length(x1));
%DISPLAY PSF
surf(x1,x2,Gaussian_PSF);
xlabel('x1'); ylabel('x2'); zlabel('bivariate normal surface');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEBLUR IMAGE:
signal_var = var(Gaussian_PSF(:));
noise_var = .001*signal_var;
deblurred_image = deconvwnr(input_image, Gaussian_PSF, noise_var/signal_var);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IDENTIFY LOCAL MAXIMA IN THE FREQUENCY SURFACE EXCEEDING THRESHOLD:
tolerance=5.0;
threshold=20.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_rows = size(input_image_perc,1);
num_cols = size(input_image_perc',1);
local_maxima = zeros(num_rows, num_cols);
for row = 1:num_rows
for col = 1:num_cols
current_val = input_image_perc(row, col);
if current_val >=threshold
neighborhood = [];
nghbr_row = row-1;
nghbr_col = col-1;
if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
neighborhood = [neighborhood input_image_perc(nghbr_row, nghbr_col)];
end
nghbr_row = row-1;
nghbr_col = col;
if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
neighborhood = [neighborhood input_image_perc(nghbr_row, nghbr_col)];
end
nghbr_row = row-1;
nghbr_col = col+1;
if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
neighborhood = [neighborhood input_image_perc(nghbr_row, nghbr_col)];
end
nghbr_row = row;
nghbr_col = col-1;
if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
neighborhood = [neighborhood input_image_perc(nghbr_row, nghbr_col)];
end
nghbr_row = row;
nghbr_col = col+1;
if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
neighborhood = [neighborhood input_image_perc(nghbr_row, nghbr_col)];
end
nghbr_row = row+1;
nghbr_col = col-1;
if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
neighborhood = [neighborhood input_image_perc(nghbr_row, nghbr_col)];
end
nghbr_row = row+1;
nghbr_col = col;
if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
neighborhood = [neighborhood input_image_perc(nghbr_row, nghbr_col)];
end
nghbr_row = row+1;
nghbr_col = col+1;
if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
neighborhood = [neighborhood input_image_perc(nghbr_row, nghbr_col)];
end

if current_val>=max(neighborhood)-tolerance
local_maxima(row, col)=1;
end
end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%turn off pixels that are not local maxima
deblurred_image = double(deblurred_image) .* double(local_maxima);
deblurred_image(deblurred_image<0.0)=0.0;
%renormalize so that total light in deblurred_image equals total light in
%input_image
deblurred_image = immultiply(deblurred_image, sum(sum(input_image))/sum(sum(deblurred_image)));
%export deblurred_image to hard disk:
geotiffwrite(output_image_pathname,deblurred_image,input_georeference);
imagesc(deblurred_image);
toc
