

%%%%Code-1- Signal Clipping Technique

clc
clear all
close all  
 
 
M = 4;                        %   QPSK signal 
data_points = 128;            %   have 128 data points
blk_size = 8;                 %   size of each ofdm block
cp_len = ceil(0.1*blk_size);  %   length of cyclic prefix
ifft_points = blk_size;       %   128 points for the FFT/IFFT
fft_points = blk_size;
 
%%%% Transmitter Part %%%%%%
data_source = randsrc(1, data_points, 0:M-1); %%%%Generating 1 x 128 vector of random data points
figure(1)
stem(data_source); 
title('Transmitted Data')
xlabel('Data Points');
ylabel('transmitted data phase representation')
grid on; 
 
%%%%%%%%%%%% QPSK modulation%%%%%%%%%
qpsk_mod_data = pskmod(data_source, M);
 
%%%%%%%%%%%% IFFT on each block %%%%%%
 
num_cols=length(qpsk_mod_data)/blk_size;%   First:Find out the number of colums that will exist after reshaping
data_matrix = reshape(qpsk_mod_data, blk_size, num_cols);
 
 
cp_start = blk_size-cp_len;             %   Second: Create empty matix to put the IFFT'd data
cp_end = blk_size;
 
 
for i=1:num_cols                        %   Third: Operate columnwise & do CP
    ifft_matrix(:,i) = ifft((data_matrix(:,i)),ifft_points);
   
    for j=1:cp_len                      %   Compute and append Cyclic Prefix
       actual_cp(j,i) = ifft_matrix(j+cp_start,i);
    end
    ifft_data(:,i) = vertcat(actual_cp(:,i),ifft_matrix(:,i)); %   Append the CP to the existing block to create the actual OFDM block
end
 
%%%%%%Converting  to serial stream for transmission%%%%%%%%%%%5
[rows_ifft_data cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;
 
%%%%%%%Actual OFDM signal to be Transmitted %%%%%%%%%%%%%%%%%
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);
figure(3)
plot(real(ofdm_signal));
title('OFDM Signal');
xlabel('Time');
ylabel('Amplitude');
grid on;
 
%%%%%%%%%%%%%% CLIPPING TECHNIQUE %%%%%%%%%%%%
 
avg=0.4;
clipped=ofdm_signal;
for i=1:length(clipped)
    if clipped(i) > avg
        clipped(i) = avg;
    end
    if clipped(i) < -avg
        clipped(i) = -avg;
    end
end
figure(4)
plot(real(clipped)); 
title('clipped Signal');
xlabel('Time'); 
ylabel('Amplitude');
grid on;
 
%%%%%%%%%%% High Pass Amplifier %%%%%%%%%%%%%%%%%%%%%%%%%
 
noise = randn(1,len_ofdm_data) +  sqrt(-1)*randn(1,len_ofdm_data);%%%%%%Generating random complex noise%%%%%%
 
%%%%%%%%%% Transmitted OFDM signal after passing through HPA
 
%%%%%%Without clipping%%%%%%%%%%%%
for i=1:length(ofdm_signal)
    if ofdm_signal(i) > avg
        ofdm_signal(i) = ofdm_signal(i)+noise(i);
    end
    if ofdm_signal(i) < -avg
        ofdm_signal(i) = ofdm_signal(i)+noise(i);
    end
end
figure(5)
plot(real(ofdm_signal)); 
title('OFDM Signal after HPA');
xlabel('Time'); 
ylabel('Amplitude');
grid on;
 
%%%%%%With clipping%%%%%%%%%%%%
avg=0.4;
for i=1:length(clipped)
    if clipped(i) > avg
        clipped(i) = clipped(i)+noise(i);
    end
    if clipped(i) < -avg
        clipped(i) = clipped(i)+noise(i);
    end
end
figure(6)
plot(real(clipped));
title('clipped Signal after HPA');
xlabel('Time');
ylabel('Amplitude');
grid on;
 
%%%%%%%%%%%% Channel %%%%%%%%%%%%%%%%%%%%
 
channel = randn(1,blk_size) + sqrt(-1)*randn(1,blk_size);  %   Creating a complex multipath channel
 
%%%%%%%%%%%% Receiver Without Clipped Signal %%%%%%%%%%%%%%%%%%
 
channel_signal = filter(channel, 1, ofdm_signal);%%%%Passing the ofdm signal through the channel
 
awgn_noise = awgn(zeros(1,length(channel_signal)),0);%%%%%% Adding AWGN Noise
 
recvd_signal = awgn_noise+channel_signal;   %%%%%%%%% Adding  whole noise to signal...
 
recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);%%%%%%Converting Data back to "parallel" form to perform FFT
 
recvd_signal_matrix(1:cp_len,:)=[]          %%%%%%%%%%%% Removing Cylic Prefix
 
 
for i=1:cols_ifft_data,
    
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),fft_points);%%%%%% Performing FFT
end
 
recvd_serial_data = reshape(fft_data_matrix, 1,(blk_size*num_cols));%%%%%%%% Converting to Serial stream
 
qpsk_demod_data = pskdemod(recvd_serial_data,M);% %%%%%%%% Demodulating the data
 
figure(7)
stem(qpsk_demod_data,'rx');
title('Received Data "X"')
xlabel('Data Points');
ylabel('received data phase representation');
grid on;
   
%%%%%%%%%%% Reciever With Clipped Signal %%%%%%%%%%%%%
 
channel_signal = filter(channel, 1, clipped);  %  Pass the ofdm signal through the channel
 
awgn_noise = awgn(zeros(1,length(channel_signal)),0);%%%%%%%% Adding AWGN Noise
 
recvd_signal = awgn_noise+channel_signal;%%%%%%%%%%Adding Whole noise to signal
 
recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);%%%%%%%%%Converting Data back to "parallel" form to perform FFT
 
recvd_signal_matrix(1:cp_len,:)=[];%%%%%%%%Removing CP
 
 
for i=1:cols_ifft_data,
    
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),fft_points);%%%%%%%Performing FFT
end
 
recvd_serial_data = reshape(fft_data_matrix, 1,(blk_size*num_cols));%%%%%%%%Converting again to serial stream
 
qpsk_demod_data = pskdemod(recvd_serial_data,M);%%%%%%%%%%%%Demodulating the data
figure(8)
stem(qpsk_demod_data,'rx');
title('Received Data clipped "X"')  
xlabel('Data Points');
ylabel('received data phase representation');
grid on;
