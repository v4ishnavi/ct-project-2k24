clc;clear;close all;
%A/D Conversion
[audio1, fs] = audioread('project.wav');

% This is the A/D converter, completely copied honestly

audio2 = audio1(:,2);
% disp(length(audio));

% <<<<<<< HEAD:main.asv
bl = 1000000; %block length
n = ceil(length(audio2)/bl);
% n = 5;
=======
bl = 100000; %block length
% n = ceil(length(audio2)/bl);
n = 5;
% >>>>>>> 1ec08f7c55e6c40f65b260c5d16b2394de004c91:rough.m
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
<<<<<<< HEAD:main.asv
    
    % figure;
    % stem(real(decoded_downsample_rcos_mem),imag(decoded_downsample_rcos_mem));
    
    final_output_rect = fourpamunmapA(decoded_downsample_rect,0);
    final_output_rcos = fourpamunmapA(decoded_downsample_rcos,0);
    
    final_output_rect_mem = fourpamunmapA(decoded_downsample_rect_mem,1);
    final_output_rcos_mem = fourpamunmapA(decoded_downsample_rcos_mem,1);
    
    count_rect = 0; 
    count_rcos = 0; 
    count_rect_mem = 0;
    count_rcos_mem = 0;
    % disp(length(binary_vector)); disp(length(final_output_rcos))
    
    % final_output_rect = binary_vector;
    
    for p = 1: length(binary_vector)
        % if binary_vector(p) ~= final_output_rect(p)
        %     final_output_rect(p) = binary_vector(p);
        % end
        if binary_vector(p) == final_output_rect(p)
            count_rect = count_rect + 1;
        end
        if binary_vector(p) == final_output_rcos(p)
            count_rcos = count_rcos + 1;
        end
        if binary_vector(p) == final_output_rect_mem(p)
            count_rect_mem = count_rect_mem + 1;
        end
        if binary_vector(p) == final_output_rcos_mem(p)
            count_rcos_mem = count_rcos_mem + 1;
        end
    end
    disp(count_rcos_mem/length(binary_vector))
    disp(count_rect_mem/length(binary_vector))
    disp(count_rcos/length(binary_vector));
    disp(count_rect/length(binary_vector));

    % op_rect_ml(outfocus) = final_output_rect;
    % op_rcos_ml(outfocus) = final_output_rcos;
    % op_rect_m_(outfocus) = final_output_rect_mem;
    % op_rcos_m_(outfocus) = final_output_rcos_mem;

    op_rect_ml(i,1:length(binary_vector)) = final_output_rect;
    op_rcos_ml(i,1:length(binary_vector)) = final_output_rcos;
    op_rect_m_(i,1:length(binary_vector)) = final_output_rect_mem;
    op_rcos_m_(i,1:length(binary_vector)) = final_output_rcos_mem;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            D/A conversion                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(op_rcos_m_);

% binary_matrix = num2str(reshape(op_rect_ml, [], 16));
% audio_integers = bin2dec(binary_matrix); 
% audio_reconstructed = typecast(uint16(audio_integers), 'int16'); 
% audio_reconstructed_normalized = double(audio_reconstructed) / 32767; 
% 
% audiowrite('audio_rectml.wav', audio_reconstructed_normalized, fs);

% sound(double(typecast(uint16(bin2dec(num2str(reshape(final_output_rect, [], 16))), 'int16'))/32767, fs))

% for l = 1:n
%     op_rect_ml = 
% end

% op_rect_ml = reshape(op_rect_ml',[1,n*length(op_rect_ml)]);
sound(double(typecast(uint16(bin2dec(num2str(reshape(op_rcos_m_, [], 16)))), 'int16'))/32767, fs)

%------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
%This section is just to check stuff/ ignore unless curious
% plot(modulated_rect_vec);
% hold off;
% subplot(2,1,1)
% plot(modulated_rcos_vec);
% subplot(2,1,2)
% hold on;
% plot(line_vector_rect);
% % ts = 1/fs;
% upsample_factor = fs;
% % 
% t = 1/m:1/m:length(symbols_vector);
% t = linspace(1, len_ip, len_ip * upsample_factor);
% lv_rect_upsampled = interp1(1:len_ip, line_vector_rect, t, 'linear');
% lv_rcos_upsampled = interp1(1:len_ip, line_vector_rcos, t, 'linear');
% len_up = length(lv_rect_upsampled);
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
=======
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
>>>>>>> 1ec08f7c55e6c40f65b260c5d16b2394de004c91:rough.m



    % % cutoff = iterative_threshold(outvector, 4);
    % % cutoff = sort(cutoff); 
    % % disp(cutoff)
    % 
    % for k = 0:2:l*2-1
    % 
    %     val = zeros(1, 4);
    %     for i = 1:4
    %     val(i) = abs(outvector(k/2 + 1) - cutoff(i));
    %     end
    %     if min(val) == val(4)
    %         invector(k+1) = '1';
    %         invector(k+2) = '1';
    %     elseif min(val) == val(3)
    %         invector(k+1) = '1';
    %         invector(k+2) = '0';
    %     elseif min(val) == val(2)
    %         invector(k+1) = '0';
    %         invector(k+2) = '1';
    %     elseif min(val) == val(1)
    %         invector(k+1) = '0';
    %         invector(k+2) = '0';
    %     end
    % 
    %     % if outvector(k/2 + 1) >= cutoff(1)
    %     %     invector(k+1) = '1';
    %     %     invector(k+2) = '1';
    %     % elseif outvector(k/2 + 1) >= cutoff(2)
    %     %     invector(k+1) = '1';
    %     %     invector(k+2) = '0';
    %     % elseif outvector(k/2 + 1) >= cutoff(3)
    %     %     invector(k+1) = '0';
    %     %     invector(k+2) = '1';
    %     % else
    %     %     invector(k+1) = '0';
    %     %     invector(k+2) = '0';
    %     % end
    %     % if 0 >= outvector(k/2 + 1) && outvector(k/2 + 1) > -2
    %     %     invector(k+1) = '0';
    %     %     invector(k+2) = '0';
    %     % elseif outvector(k/2 + 1) > 0 && outvector(k/2 +1) <= 2
    %     %     invector(k+1) = '0';
    %     %     invector(k+2) = '1';
    %     % elseif outvector(k/2 + 1) <= -2 
    %     %     invector(k+1) = '1';
    %     %     invector(k+2) = '0';
    %     % else
    %     %     invector(k+1) = '1';
    %     %     invector(k+2) = '1';
    %     % end