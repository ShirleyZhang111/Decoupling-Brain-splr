function data = surrogate_data(inCfg)
% Generate surrogate data

% input configuration
data_path = inCfg.data_path;
file_name = inCfg.file_name;
save_path = inCfg.save_path;
if isfield(inCfg,'T');T = inCfg.T;else;T = [];end
if ~contains(file_name,'.mat');file_name = [file_name '.mat'];end

% start processing
load([data_path file_name],'TC');% N*T
TC = TC(:,6:end-5);
% if isempty(T);T = 5:size(TC,2)-4;end
data = fft(TC,[],2);
absData = abs(data);
phaseData = angle(data);
phaseDataRand = phaseData;
for i_row = 1:size(data,1)
    phaseDataRand(i_row,1:600) = rand(600,1).*(2*pi)-pi; % random phase values -pi~pi
    for iInd = 1:600
        phaseDataRand(i_row,600+iInd) = -phaseDataRand(i_row,600-iInd+1);
    end
end
data = absData.*exp(1i*phaseData+1i*phaseDataRand);
data = real(ifft(data,[],2));

save([save_path file_name(4:end-4) '_surrogate.mat'], 'data');