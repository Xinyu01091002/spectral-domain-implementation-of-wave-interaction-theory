function [kx, ky, KX, KY, dx, dy] = compute_kxky(x_vec, y_vec)
% compute_kxky  由空间坐标向量生成对应的 kx, ky（FFT/谱分析用）
%
% 输入:
%   x_vec  - x 方向的一维坐标向量（必须等间距）
%   y_vec  - y 方向的一维坐标向量（必须等间距）
%
% 输出:
%   kx, ky - 分别为 x 和 y 方向的 1D 波数向量（已用 fftshift 中心化，单位 rad/m）
%   KX, KY - 通过 meshgrid(kx,ky) 生成的 2D 波数网格（大小 ny x nx）
%   dx, dy - x_vec, y_vec 的格距
%
% 说明:
%  - 频谱/FFT 分析时常用此类中心化 k 网格配合 fftshift/ifftshift
%  - 本函数假定输入向量均从小到大且等距

% 计算尺寸与格距
nx = numel(x_vec);
ny = numel(y_vec);
if nx < 2 || ny < 2
    error('x_vec 和 y_vec 必须至少包含两个点');
end

dx = x_vec(2) - x_vec(1);
dy = y_vec(2) - y_vec(1);

% 采样频率（空间频率）基准
fx = (-floor(nx/2):ceil(nx/2)-1) / (nx * dx);
fy = (-floor(ny/2):ceil(ny/2)-1) / (ny * dy);

% 转换为空间波数 (rad/m)
kx = 2 * pi * fx;
ky = 2 * pi * fy;

% 生成 2D 网格（注意 meshgrid 的顺序：meshgrid(x,y) -> size = [length(y), length(x)])
[KX, KY] = meshgrid(kx, ky);

end
