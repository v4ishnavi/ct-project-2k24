clc;clear;close all;
%A/D Conversion
% [audio1, fs] = audioread('project.wav');

%This is the A/D converter, completely copied honestly

% audio = audio1(:,1);
% audio_normalized = int16(audio * 32767); 
% audio_binary = dec2bin(typecast(audio_normalized(:), 'uint16'), 16); 
% binary_vector = audio_binary(:)';

%------------------------------------------------------------------------------
% l = length(binary_vector);
k = 8; 
m = 4; %mapping to 4 bits 

%This section of the code encodes the binary vector using 4-PAM 
% symbols_vector = fourpammapA(binary_vector);

symbols_vector = [1 1 1 3 3 -3 -1 1 1 1 3 3 -3 -1 1 1 1 3 3 -3 -1];
%-------------------------------------------------------------------------------
%LINE CODING 
symbols_vector_up = zeros(1,m*k);
for j = 1:k
    symbols_vector_up(j + (j-1)*(m-1)) = symbols_vector(j);
end %you could have just used upsample, but okie
% Rectangular pulse
p1 = [ones(1,k/2)];

%raised cosine pulse 
a = 0.5;
[p2, ~] = raised_cosine(a, 4, k/2);

line_vector_rect = conv(symbols_vector_up,p1,'same');
line_vector_rcos = conv(symbols_vector_up, p2, 'same');
figure;
plot(1:m*k,line_vector_rect);
hold on
plot(1:m*k, line_vector_rcos);
hold on
stem(1:m*k, symbols_vector_up);
hold off
%-------------------------------------------------------------------------------------
%MODULATION PART: here there is no doppler effect so dsb-sc
%rectangular encoded: 
fc = 1e6;
len_ip = length(symbols_vector_up);
fs = 44100;%rough: for now. 
% ts = 1/fs;
upsample_factor = fs;

t = linspace(1, len_ip, len_ip * upsample_factor);
lv_rect_upsampled = interp1(1:len_ip, line_vector_rect, t, 'linear');
lv_rcos_upsampled = interp1(1:len_ip, line_vector_rcos, t, 'linear');
% len_up = length(lv_rect_upsampled);
% t = linspace(0, (len_up-1)/fs,len_up);

carrier = cos(fc*2*pi*t);
% disp(length(t)); disp(length(t));
modulated_rect_vec = lv_rect_upsampled .* carrier;
modulated_rcos_vec = lv_rcos_upsampled .* carrier;

figure;
subplot(4,1,1)
plot(modulated_rect_vec);
subplot(4,1,2)
plot(lv_rect_upsampled);
subplot(4,1,3)
plot(modulated_rcos_vec);
subplot(4,1,4)
plot(lv_rcos_upsampled);
%---------------------------------------------------------------------------------------
%This section is just to check stuff/ ignore unless curious
% modulation failures 
% line_vector_rect = repelem(line_vector_rect, fs);
% line_vector_rcos = repelem(line_vector_rcos, fs); %i think we need to upsample this using linear upsampling. 
% disp(length(line_vector_rcos)); disp(length(line_vector_rect)); disp(length(carrier));
% 
% binary_vectorA = fourpamunmapA(symbols_vector);
% 
% binary_matrix = reshape(binary_vectorA, [], 16);
% audio_integers = bin2dec(binary_matrix); 
% audio_reconstructed = typecast(uint16(audio_integers), 'int16'); 
% audio_reconstructed_normalized = double(audio_reconstructed) / 32767; 
% 
% audiowrite('reconstructed_project.wav', audio_reconstructed_normalized, fs);
% % 
% sound(audio1, fs); 
% pause(length(audio)/fs + 1);
% sound(audio_reconstructed_normalized, fs);