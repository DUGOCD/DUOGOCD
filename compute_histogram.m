% 计算图像的直方图
function [P, average] = compute_histogram(image)
    [h, w] = size(image);
    MN = h * w;
    P = zeros(1, 256);
    
    for i = 1:h
        for j = 1:w
            pixel_value = image(i, j) + 1; % MATLAB 索引从 1 开始
            P(pixel_value) = P(pixel_value) + 1;
        end
    end
    
    P = P / MN; % 归一化
    average = sum((0:255) .* P); % 计算平均值
end