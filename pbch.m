enb.CyclicPrefix = 'Normal';
enb.PHICHDuration = 'Normal';
enb.Ng = 'Sixth';
enb.NDLRB = 25;
enb.CellRefP = 1
enb.DuplexMode = 'FDD';
enb.NCellID = 0;
enb.NSubframe = 0;
enb.CFI = 1;
enb.NFrame = 828;

mib_bit=zeros(1,24);
if enb.NDLRB == 6
    mib_bit(1,1:3) = [0 0 0];  
elseif enb.NDLRB == 15
    mib_bit(1,1:3) = [0 0 1];
elseif enb.NDLRB == 25
    mib_bit(1,1:3) = [0 1 0];
elseif enb.NDLRB == 50
    mib_bit(1,1:3) = [0 1 1]; 
elseif  enb.NDLRB == 75
    mib_bit(1,1:3) = [1 0 0];
elseif enb.NDLRB == 100
    mib_bit(1,1:3) = [1 0 1];
end

if enb.PHICHDuration == 'Normal'
    mib_bit(1,4) = (0);
    
elseif enb.PHICHDuration == 'Extended'
    mib_bit(1,4) = [1];
end
if  enb.Ng == 'Sixth'
    mib_bit(1,5:6) = [0 0];
   
elseif enb.Ng == 'Half'
    mib_bit(1,5:6) = [0 1];
    
elseif enb.Ng == 'One'
    mib_bit(1,5:6) = [1 0];
    
elseif enb.Ng == 'Two'
    mib_bit(1,5:6) = [1 1];   
end
A = dec2bin(enb.NFrame)
result = reshape(str2num(A(:)), size(A)) %char' in8'e dönüştürmek için 
result_NFrame = int8(result) %char' in8'e dönüştürmek için 
mib_bit(7:16) = result_NFrame 

%mib_bit_crc = lteCRCEncode(mib_bit,'16',0);%crc bitini ekleme
mib_bit = [0,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0]

  crc_poly = [1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1]
  crc_len  = 16;
  
    crc_rem   = zeros(1,crc_len+1);
    tmp_array = [mib_bit, zeros(1,crc_len)];

    % Loop to calculate CRC bits
    for(n=1:length(mib_bit)+crc_len)
        for(m=1:crc_len)
            crc_rem(m) = crc_rem(m+1);
        end
        crc_rem(crc_len+1) = tmp_array(n);

        if(crc_rem(1) ~= 0)
            for(m=1:crc_len+1)
                crc_rem(m) = mod(crc_rem(m)+crc_poly(m), 2);
            end
        end
    end

    for(n=1:crc_len)
        crc_bit(n) = crc_rem(n+1);
        
    end
    mib_bit_crc = [mib_bit crc_bit]; %mib biti ile crc bitini toplar ve 40 bitlik total bit oluşur

Mib_Convolutional = lteConvolutionalEncode(mib_bit_crc)%40*3
Mib_RateMatching = repmat(Mib_Convolutional,16,1)%120*16

random_data = randi([0 1],1920,1); %randi fonksiyonu ile random data ürettim
scrambler_data = xor(random_data,Mib_RateMatching); %random data ile rateMatchin xor yapılıp scrambler data oluşturulur
descrambler_data = xor(random_data,scrambler_data) %random data ile geri elde edilir
isequal(descrambler_data,Mib_RateMatching) %sağlama yapılır 

n = 1920
index=1;
for i=1:2:n-1;
    if(scrambler_data(i) == 0 && scrambler_data(i+1) == 0)
        qpsk_signal(index) = 0.707+0.707j;
        
    elseif(scrambler_data(i) == 0 && scrambler_data(i+1) == 1)   
        qpsk_signal(index) = 0.707-0.707j;

    elseif(scrambler_data(i) == 1 && scrambler_data(i+1) == 0)
       qpsk_signal(index) = -0.707+0.707j;

    elseif(scrambler_data(i) == 1 && scrambler_data(i+1) == 1)
       qpsk_signal(index) = -0.707-0.707j;

    end
    index = index + 1;
end

N = length(qpsk_signal)
SNR_dB_vec = 0:5:30;

for snr_idx = 1:length(SNR_dB_vec)
SNR_dB = SNR_dB_vec(snr_idx);% Sembol Enerjisini Hesapla
Eavg = sum(abs(qpsk_signal) .^ 2)/N;% SNR'yi (dB cinsinden) SNR'ye (Doğrusal olarak) dönüştür
SNR_lin = 10 .^ (SNR_dB/10);% AWGN'nin Sigmasını (Standart Sapma) hesaplayın
awgnSigma = sqrt(Eavg/(2*SNR_lin));% Normal Dağılım ile bir gürültü dizisi oluşturun ve bunu sigma ile yeniden ölçeklendirin
noise = awgnSigma*(randn(1,N)+j*randn(1,N));% Gürültüyü orijinal sinyale ekleyin
Noise_corrupted_signal = qpsk_signal + noise;

B = [];
Received_Symbol = [];
for i = 1:size(Noise_corrupted_signal,2)
    Received_Signal_Real = real(Noise_corrupted_signal(1,i));
    Received_Signal_Imaginary = imag(Noise_corrupted_signal(1,i));
    if (Received_Signal_Real >= 0 && Received_Signal_Imaginary >= 0)
        Received_Symbol = [0 0];
    elseif (Received_Signal_Real > 0 && Received_Signal_Imaginary < 0)
        Received_Symbol = [0 1];
    elseif (Received_Signal_Real < 0 && Received_Signal_Imaginary < 0)
        Received_Symbol = [1 1];
    else
        Received_Symbol = [1 0];
    end
    B = [B Received_Symbol];
    All_Received_Symbols = transpose(B);
end
%BER CALCULATİON
bits_flipped = 0;
idx = 1;
m = 1920;
for ik = 1:m
    total_bits = m;
    if (scrambler_data(ik) ~= All_Received_Symbols(ik))
        bits_flipped = bits_flipped + 1;
    end
        idx = idx + 1;
end
Bit_Error_Rate(snr_idx) = bits_flipped/total_bits;
C = xor(scrambler_data,All_Received_Symbols);

[num,ratio(snr_idx)] = biterr(All_Received_Symbols,scrambler_data)
end
semilogy(SNR_dB_vec, Bit_Error_Rate),grid,xlabel('SNR dB'),ylabel('BER')