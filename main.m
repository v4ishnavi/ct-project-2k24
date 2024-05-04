clc;clear;close all;

[audio1, fs] = audioread('project.wav');

%This is the A/D converter, completely copied honestly
audio = audio1(:,1);
audio_normalized = int16(audio * 32767); 
audio_binary = dec2bin(typecast(audio_normalized(:), 'uint16'), 16); 
binary_vector = audio_binary(:)';


l = length(binary_vector);
k = 8;
m = 4;

%This section of the code encodes the binary vector using 4-PAM 
% symbols_vector = fourpammapA(binary_vector);
symbols_vector = [1 1 1 3 3 -3 -1 1 1 1 3 3 -3 -1 1 1 1 3 3 -3 -1];

symbols_vector_up = zeros(1,m*k);
for j = 1:k
    symbols_vector_up(j + (j-1)*(m-1)) = symbols_vector(j);
end

%and that's it hehe

% figure;
% stem(1:l/2,symbols_vector);

% Rectangular pulse
p = [ones(1,k/2)];
% figure;
% stem(1:k*2,p);

line_vector = conv(symbols_vector_up,p,'same');

% figure;
% plot(1:32,line_vector);

%This section is just to check stuff/ ignore unless curious
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