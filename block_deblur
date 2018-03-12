function result = block_deblur(block, perc_block, pixel_length, pixel_width, use_default_noise_sig_ratio)
    freq_filter = 1;
    threshold=20.0;
    default_noise_sig_ratio = 0.011;
    
    blurred_image = block;
    blurred_image_perc = perc_block;
    
    %set negative values (if there are any) to zero:
    blurred_image(blurred_image<0.0) = 0.0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ESTIMATE NOISE-SIGNAL RATIO: (assumes image is >10 rows and >10 columns)
    num_rows = size(blurred_image,1);
    num_cols = size(blurred_image',1);

    double_blurred_image = im2double(blurred_image);
    double_interior_blurred_image = double_blurred_image(5:num_rows-5, 5:num_cols-5);
    var_interior = var(double_interior_blurred_image(:));

    %calculate variance of periperhy:
    top_five_rows = double_blurred_image(1:5, :);
    bottom_five_rows = double_blurred_image(num_rows-4:num_rows, :);
    five_left_columns = double_blurred_image(6:num_rows-5, 1:5);
    five_right_columns = double_blurred_image(6:num_rows-5, num_cols-4:num_cols);
    vector_of_exterior_pixels = [top_five_rows(:); bottom_five_rows(:); five_left_columns(:); five_right_columns(:)];
    var_periphery = var(vector_of_exterior_pixels);

    noise_sig_ratio = var_periphery/var_interior;
    if use_default_noise_sig_ratio==1
        noise_sig_ratio = default_noise_sig_ratio;
    end
    %signal_var = var(Gaussian_PSF(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigmas = [1.0:.1:4.0];
    normalizers = [];
    for sigma = sigmas
    %disp('now deblurring using this sigma:')
    %sigma
    sigma_x2 = sigma^2;%default value: 2.5^2
    sigma_y2 = sigma^2;%default value: 2.5^2
    sigma_x2 = sigma_x2/(pixel_width^2); %converts sigma_x2 from units of km to units of bins
    sigma_y2 = sigma_y2/(pixel_length^2); %converts sigma_y2 from units of km to units of bins

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
    %signal_var = var(Gaussian_PSF(:));
    deblurred_image = deconvwnr(blurred_image, Gaussian_PSF, noise_sig_ratio);
    if sum(sum(deblurred_image(deblurred_image<0)))<0
        disp('neg vals after deconvwnr')
        sum(sum(deblurred_image(deblurred_image<0)))
        pause
    end
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_rows = size(blurred_image_perc,1);
    num_cols = size(blurred_image_perc',1);
    local_maxima = zeros(num_rows, num_cols);
    for row = 1:num_rows
    for col = 1:num_cols
    current_val = blurred_image_perc(row, col);
    %if current_val >=threshold
    neighborhood = [];
    nghbr_row = row-1;
    nghbr_col = col-1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row-1;
    nghbr_col = col;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row-1;
    nghbr_col = col+1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row;
    nghbr_col = col-1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row;
    nghbr_col = col+1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row+1;
    nghbr_col = col-1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row+1;
    nghbr_col = col;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row+1;
    nghbr_col = col+1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end

    if current_val>=max(neighborhood)-tolerance
    local_maxima(row, col)=1;
    end
    %end
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
    total_light_before_freq_filter = sum(sum(deblurred_image));
    %turn off pixels that are not local maxima
    if freq_filter==1
    deblurred_image = double(deblurred_image) .* double(local_maxima);
    end
    deblurred_image(deblurred_image<0.0)=0.0;
    %normalize so that total light in deblurred_image equals total light in
    if sum(sum(deblurred_image)) ~= 0
    normalizer = total_light_before_freq_filter/sum(sum(deblurred_image));
    normalizers = [normalizers normalizer];
    deblurred_image = immultiply(deblurred_image, normalizer);
    else
    %deblurred_image contain zero light
    normalizer = 1.0;
    normalizers = [normalizers normalizer];
    end
    
    %APPLY FREQUENCY THRESHOLD: (DEFAULT: sources lit <20% of cloud-free nights should be switched off)
    lit_often = blurred_image_perc;
    lit_often(blurred_image_perc<threshold) = 0;
    lit_often(blurred_image_perc>=threshold) = 1;
    deblurred_image = double(deblurred_image) .* double(lit_often);
    end

    %calculate best sigma for this noise-signal ratio:
    index_of_best_sigma = find(normalizers==min(normalizers));
    index_of_best_sigma = index_of_best_sigma(1);
    best_sigma = sigmas(index_of_best_sigma);
    
    %store best sigma and its normalizer for outer loop:
    disp('results from looping thru sigmas:')
    noise_sig_ratio
    best_sigma
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %NOW THAT WE KNOW THE BEST SIGMA, WE DEBLUR ONCE USING THAT SIGMA, AND SAVE
    %THE RESULT TO DISK

    sigma = best_sigma;
    
    sigma_x2 = sigma^2;%default value: 2.5^2
    sigma_y2 = sigma^2;%default value: 2.5^2
    sigma_x2 = sigma_x2/(pixel_width^2); %converts sigma_x2 from units of km to units of bins
    sigma_y2 = sigma_y2/(pixel_length^2); %converts sigma_y2 from units of km to units of bins

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
    %signal_var = var(Gaussian_PSF(:));
    deblurred_image = deconvwnr(blurred_image, Gaussian_PSF, noise_sig_ratio);
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
    num_rows = size(blurred_image_perc,1);
    num_cols = size(blurred_image_perc',1);
    local_maxima = zeros(num_rows, num_cols);
    for row = 1:num_rows
    for col = 1:num_cols
    current_val = blurred_image_perc(row, col);
    %if current_val >=threshold
    neighborhood = [];
    nghbr_row = row-1;
    nghbr_col = col-1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row-1;
    nghbr_col = col;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row-1;
    nghbr_col = col+1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row;
    nghbr_col = col-1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row;
    nghbr_col = col+1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row+1;
    nghbr_col = col-1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row+1;
    nghbr_col = col;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end
    nghbr_row = row+1;
    nghbr_col = col+1;
    if nghbr_row>=1 & nghbr_col>=1 & nghbr_row<=num_rows & nghbr_col<=num_cols
    neighborhood = [neighborhood blurred_image_perc(nghbr_row, nghbr_col)];
    end

    if current_val>=max(neighborhood)-tolerance
    local_maxima(row, col)=1;
    end
    %end
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

    total_light_before_freq_filter = sum(sum(deblurred_image));
    %turn off pixels that are not local maxima
    if freq_filter==1
    deblurred_image = double(deblurred_image) .* double(local_maxima);
    end
    deblurred_image(deblurred_image<0.0)=0.0;
    %normalize so that total light in deblurred_image equals total light in
    normalizer = total_light_before_freq_filter/sum(sum(deblurred_image));
    deblurred_image = immultiply(deblurred_image, normalizer);

    %APPLY FREQUENCY THRESHOLD: (DEFAULT: sources lit <20% of cloud-free nights should be switched off)
    lit_often = blurred_image_perc;
    lit_often(blurred_image_perc<threshold) = 0;
    lit_often(blurred_image_perc>=threshold) = 1;
    deblurred_image = double(deblurred_image) .* double(lit_often);
    
    noise_sig_ratio
    best_sigma
    result = deblurred_image;
