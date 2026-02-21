% function xA = analytic2D(x)
% % ANALYTIC2D  计算二维实信号的解析信号
% %
% %   xA = analytic2D(x)
% %
% % 输入：
% %   x   — M×N 实值矩阵（例如一个二维波面或图像）
% %
% % 输出：
% %   xA  — M×N 复值矩阵，即 x 的二维解析信号
% %
% % 说明：
% %   二维解析信号的定义：在频域 PSD（二维傅里叶变换）上，保留“半边”频率分量，
% %   令对应谱值加倍，丢弃另一半（负频率）；同时 DC 分量保持不变。具体做法：
% %     1. 计算 X = fft2(x)，得到未中心化频谱。
% %     2. 做中心化： X0 = fftshift(X)。此时频谱中心位于 (M_c, N_c)。
% %     3. 构造掩码 H（大小 M×N），对于频率坐标 (ky,kx)：
% %          若 (kx > 0) 或者 (kx == 0 且 ky > 0)，则 H = 2；否则 H = 0。
% %        对 DC 点 (kx=0, ky=0) 以 H = 1 处理。
% %     4. X_analytic_shifted = H .* X0；然后做逆中心化： X_analytic = ifftshift(X_analytic_shifted)。
% %     5. xA = ifft2(X_analytic)，直接得到复值解析信号矩阵。
% %
% %   返回的 xA 是复数，实部等于原始信号 x，虚部是其二维 Hilbert 变换。
% %
% % 示例：
% %   img = peaks(256);      % 256×256 实矩阵
% %   imgA = analytic2D(img);
% %   % now imag(imgA) is the 2D Hilbert transform of img.
% %
% % 作者：ChatGPT
% % 日期：2025-06-05
% 
%     [M, N] = size(x);
% 
%     % 1. 原始二维 FFT
%     X = fft2(x);                 % M×N 复谱（未中心化）
% 
%     % 2. 中心化
%     X0 = fftshift(X);            % 中心化谱：零频在 (M_c, N_c)
% 
%     % 坐标偏移量，中心索引
%     Mc = floor(M/2) + 1;         % 第 M/2+1 行是 ky = 0
%     Nc = floor(N/2) + 1;         % 第 N/2+1 列是 kx = 0
% 
%     % 3. 构造掩码 H
%     %    对每个频率位置 (i,j)，其 ky = i - Mc，kx = j - Nc
%     %    保留条件： kx > 0  或者  (kx == 0 且 ky > 0)
%     [J, I] = meshgrid(1:N, 1:M); % I(i,j)=i，J(i,j)=j
%     kx = J - Nc;                 % 大小 M×N
%     ky = I - Mc;                 % 大小 M×N
% 
%     H = zeros(M, N);
%     % 条件 kx > 0  或  (kx == 0 且 ky > 0)
%     mask = (kx > 0) | ( (kx == 0) & (ky > 0) );
%     H(mask) = 2;
% 
%     % DC 点 (kx=0, ky=0) 单独设为 1
%     H(Mc, Nc) = 1;
% 
%     % 4. 应用掩码并逆中心化
%     X0_analytic = H .* X0;                       % 只保留半边频率并加倍
%     X_analytic  = ifftshift(X0_analytic);        % 恢复到未中心化频谱
% 
%     % 5. 逆傅里叶变换得到解析信号
%     xA = ifft2(X_analytic);                      % M×N 复数矩阵
% 
% end
function xA = analytic2D(x)
    S_temp=fft2(x);
    % S_temp(1:end/2+1,:)=0;
S_temp(:,1:end/2+1)=0;
    xA=ifft2(S_temp)*2;
end