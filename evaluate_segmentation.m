% �����ָ���
function [psnr_value, ssim_value, uiqi_value, hpsi_value] =evaluate_segmentation(original_image, segmented_image)
    % ���� PSNR
    psnr_value = psnr(segmented_image, original_image);
    fprintf('PSNR: %.4f dB\n', psnr_value);
    
    % ���� SSIM
    ssim_value = ssim(segmented_image, original_image);
    fprintf('SSIM: %.4f\n', ssim_value);
    
    % ���� UIQI
    uiqi_value = universal_image_quality_index(original_image, segmented_image);
    fprintf('UIQI: %.4f\n', uiqi_value);
    
    % ���� HPSI
    hpsi_value = haar_wavelet_similarity(original_image, segmented_image);
    fprintf('HPSI: %.4f\n', hpsi_value);
end

% UIQI ����
function uiqi_value = universal_image_quality_index(original_image, segmented_image)
    original_image = double(original_image);
    segmented_image = double(segmented_image);
    
    mean_original = mean(original_image(:));
    mean_segmented = mean(segmented_image(:));
    
    var_original = var(original_image(:));
    var_segmented = var(segmented_image(:));
    
    cov_value = cov(original_image(:), segmented_image(:));
    cov_value = cov_value(1, 2);
    
    uiqi_value = (4 * cov_value * mean_original * mean_segmented) / ...
                 ((var_original + var_segmented) * (mean_original^2 + mean_segmented^2));
end

% HPSI ����
function hpsi_value = haar_wavelet_similarity(original_image, segmented_image)
    % ʹ�� Haar С���任����������
    [c_original, ~] = wavedec2(original_image, 1, 'haar');
    [c_segmented, ~] = wavedec2(segmented_image, 1, 'haar');
    
    % ����������
    hpsi_value = corr2(c_original, c_segmented);
end