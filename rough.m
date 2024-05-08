clc;clear;close all;
%A/D Conversion
[audio1, fs] = audioread('project.wav');

% This is the A/D converter, completely copied honestly

audio2 = audio1(:,1);
% disp(length(audio));

bl = 100000; %block length
% n = ceil(length(audio2)/bl);
n = 5;
%number of fullsize block lengths
% bl = 100;
% n = 1;
% op_rect_ml = zeros(1,16*n*bl); %16 for 16 binary bits?
% op_rcos_ml = zeros(1, 16*n*bl);
% op_rect_m_ = zeros(1, 16*n*bl);
% op_rcos_m_ = zeros(1, 16*n*bl);

% op_rect_ml = zeros(1, 16*length(audio2)); %16 for 16 binary bits?
% op_rcos_ml = zeros(1, 16*length(audio2));
% op_rect_m_ = zeros(1, 16*length(audio2));
% op_rcos_m_ = zeros(1, 16*length(audio2));

op_rect_ml = '';
op_rcos_ml = '';
op_rect_m_ = '';
op_rcos_m_ = '';

for i = 1:n
    
    if i == n
        infocus = 1+(i-1)*bl:length(audio2);
        outfocus = 1+(i-1)*16*bl:16*length(audio2);
    else
        infocus = 1+(i-1)*bl:bl+(i-1)*bl;
        outfocus = 1+(i-1)*16*bl:16*bl+(i-1)*16*bl;         
    end
    % disp(infocus(1) + ":" + infocus(length(infocus)));
    % disp(outfocus(1) + ":" + outfocus(length(outfocus)));
    
    length(outfocus)

    audio = audio2(infocus);
    
    audio_normalized = int16(audio * 32767); 
    audio_binary = dec2bin(typecast(audio_normalized(:), 'uint16'), 16); 
    binary_vector = audio_binary(:)';
    % disp(binary_vector);
    
    %------------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             encoding                                %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    l = length(binary_vector);
    k = 8; 
    m = 4; %mapping to 4 bits 
    
    %This section of the code encodes the binary vector using 4-PAM 
    symbols_vector = fourpammapA(binary_vector);
    % disp(length(symbols_vector))
    % symbols_vector = [1 1 1 3 3 -3 -1 1 1 1 3 3 -3 -1 1 1 1 3 3 -3 -1];
    
    %-------------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             line encoding                           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % disp(l);
    symbols_vector_up = zeros(1,m*l/2);
    for j = 1:l/2
        symbols_vector_up(j + (j-1)*(m-1)) = symbols_vector(j);
    end %you could have just used upsample, but okie
    % Rectangular pulse
    p1 = [ones(1,k/2)];
    
    %raised cosine pulse 
    a = 0.5;
    [p2, ~] = raised_cosine(a, 4, k/2);
    
    line_vector_rect = conv(symbols_vector_up,p1,'same');
    line_vector_rcos = conv(symbols_vector_up, p2, 'same');
    % figure;
    % plot(line_vector_rect);
    % hold on
    % % plot(line_vector_rcos);
    % % hold on
    % plot(1:m*l/2, symbols_vector_up);
    % hold off
    %-------------------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             modulation                              %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %rectangular encoded: 
     fc = 1e6;
    len_ip = length(symbols_vector_up); 
    %Tb = Ts*logM
    t = (0:len_ip-1)/(0.5*fs);
    carrier = cos(fc*2*pi*t);
    % disp(length(t)); disp(length(t));
    % disp(length(line_vector_rect)); disp(length(carrier))
    modulated_rect_vec = line_vector_rect .* carrier;
    modulated_rcos_vec = line_vector_rcos .* carrier;
    
    % figure;
    % plot(t,modulated_rect_vec);
    % plot(t,modulated_rcos_vec);
    
    % figure; 
    % subplot(2,1,1)
    % stem(1:length(fft(line_vector_rect)),abs(fft(line_vector_rect)))
    % subplot(2,1,2)
    % stem(1:length(fft(modulated_rect_vec)), abs(fft(modulated_rect_vec)))
    
    %_---------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                                 channel                             %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    snr_val = 10000;
    % Channel memoryless: 
    % assume channel noise has SNR = 3
    channel_op_rcos = awgn(modulated_rcos_vec, snr_val);
    channel_op_rect = awgn(modulated_rect_vec, snr_val); 
    % figure; 
    % subplot(2,1,1)
    % plot(channel_op_rcos);
    % subplot(2,1,2);
    % plot(channel_op_rect);
    % plot(channel_op_rcos);
    
    
    %channel memory: 
    % Parameters
    a = 0.5; 
    Tb = 2/fs;
    b = 5;
    t_h = 0:Tb:(10*Tb);  % Covering 0 to 2Tb for two delta functions
    h = a * (t_h == 0) + (1 - a) * (t_h == b*Tb);  % Impulse response
    
    % ch_m_o_rcos = conv(modulated_rcos_vec, h, 'same'); 
    % chmo_rcos = awgn(ch_m_o_rcos, snr_val);
    
    ch_m_o_rect = conv(modulated_rect_vec, h, 'same');
    chmo_rect = awgn(ch_m_o_rect, snr_val);
    
    ch_m_o_rcos = conv(modulated_rcos_vec, h, 'same');
    chmo_rcos = awgn(ch_m_o_rcos, snr_val);
    
    % figure; 
    % subplot(3,1,1);
    % stem(modulated_rcos_vec);
    % subplot(3,1,2);
    % stem(ch_m_o_rcos);
    % subplot(3,1,3);
    % stem(chmo_rcos);
    
    
    % figure; 
    % subplot(3,1,1);
    % plot(modulated_rect_vec);
    % subplot(3,1,2);
    % plot(ch_m_o_rect);
    % subplot(3,1,3);
    % plot(chmo_rect);
    
    %---------------------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             demodulation                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % demod_rect_vec = lowpass(channel_op_rect.*carrier, fc, 2*fc);
    demod_rect_vec = lowpass(channel_op_rect.*carrier,fc,2*fc);
    demod_rcos_vec = lowpass(channel_op_rcos.*carrier, fc, 2*fc); 
    % figure;
    % plot(t, demod_rect_vec);
    % plot(t, demod_rcos_vec);
    
    demod_rect_vec_mem = lowpass(chmo_rect.*carrier, fc, 2*fc);
    demod_rcos_vec_mem = lowpass(chmo_rcos.*carrier, fc, 2*fc);
    
    
    %----------------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             line decoding                           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    line_decoded_rect_vec = conv(demod_rect_vec, p1, 'same');
    % line_decoded_rect_vec = 3*line_decoded_rect_vec/max(line_decoded_rect_vec);
    
    line_decoded_rcos_vec = conv(demod_rcos_vec, p2, 'same');
    % line_decoded_rcos_vec = 3*line_decoded_rcos_vec/max(line_decoded_rcos_vec);
    % figure;
    % plot(line_decoded_rect_vec)
    % plot(line_decoded_rcos_vec);
    
    line_decoded_rect_vec_mem = conv(demod_rect_vec_mem,p1,'same');
    line_decoded_rcos_vec_mem = conv(demod_rect_vec_mem,p2,'same');
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                            4PAM decoding                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_rect = (line_decoded_rect_vec);
    decoded_downsample_rect = downsample(decoded_rect, m);

    decoded_rcos = (line_decoded_rcos_vec);
    decoded_downsample_rcos = downsample(decoded_rcos, m);
    
    decoded_rect_mem = (line_decoded_rect_vec_mem);
    decoded_downsample_rect_mem = downsample(decoded_rect_mem,m);
    
    decoded_rcos_mem = (line_decoded_rcos_vec_mem);
    decoded_downsample_rcos_mem = downsample(decoded_rcos_mem,m);
    nc = 3;
    X = zeros(length(decoded_downsample_rect),nc);
    X(:, 1) = decoded_downsample_rect;
    if i == 1
        stem(decoded_downsample_rect);
    end
    disp(decoded_downsample_rect(1:100))
    threshold = iterative_threshold(X,nc);
    disp(threshold)
    % final_output_rect = fourpamunmapA(decoded_downsample_rect,0);
    % final_output_rcos = fourpamunmapA(decoded_downsample_rcos,0);
    % 
    % final_output_rect_mem = fourpamunmapA(decoded_downsample_rect_mem,1);
    % final_output_rcos_mem = fourpamunmapA(decoded_downsample_rcos_mem,1);
end