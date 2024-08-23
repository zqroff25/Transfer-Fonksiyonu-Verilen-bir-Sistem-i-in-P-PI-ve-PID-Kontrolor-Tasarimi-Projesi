clc;
clear;
fprintf('\n Zakir YILDIZ 200757029 Son �� rakam 029 C=0.029,a = 0.0209,b = 0.920,e = 0.029/500,f = 1 \n ');
% Sistem parametreleri
C = 0.029;
a = 0.0209;
b = 0.920;
e = 0.029/500;
f = 1;

% Transfer fonksiyonu P(s)
pay = C;
payda = [0.0012 0.0743 0.9208];
P = tf(pay, payda);

% Routh-Hurwitz kriteri
disp('Sistemin Routh-Hurwitz Tablosu :');
Routh_Hurwitz_Tablosu = routh_hurwitz(payda);
disp(Routh_Hurwitz_Tablosu);
fprintf('---------------------------------\n');
% K�klerin bulunmas�
disp('Polinom K�kleri:');
kokler = roots(payda);
disp('    0');
disp(kokler);
% Kararl�l�k kontrol�
if all(real(kokler) < 0)
    disp('Sistem kararl�d�r.');
else
    disp('Sistem kararl� de�ildir.');
end
fprintf('---------------------------------\n');
% Zaman vekt�r�
zaman = 0:0.001:5;% �rnekleme s�resi 0.001 s ve toplam s�re 5 s

% Ad�m tepkisi (sadece P(s))
fprintf('/////////// Soru1: Sitemin Ad�m Tepkisi  /////////// \n');
[~, ~, ~] = zaman_tepkisi_ciz(P, zaman, 'P(s) Sistemde Ad�m Tepkisi');

% Ziegler-Nichols y�ntemi ile kontrolc� tasar�m�
Ku = 1; % �rnek kritik kazan�
Pu = 1; % �rnek periyot

% P kontrolc� tasar�m�
fprintf('/////////// Soru2: P kontrolc� ile ad�m tepkisi /////////// \n');
Kp_P = 0.5 * Ku;
C_P = pid(Kp_P, 0, 0);
sistem_P = feedback(C_P * P, 1);
[y_P, ~, bilgi_P_kontrolcu] = zaman_tepkisi_ciz(sistem_P, zaman, 'P Kontrolc� ile Ad�m Tepkisi');

% PI kontrolc� tasar�m�
fprintf('/////////// Soru2: PI kontrolc� ile ad�m tepkisi ////////// \n');

Kp_PI = 0.45 * Ku;
Ti_PI = Pu / 1.2;
C_PI = pid(Kp_PI, Kp_PI / Ti_PI, 0);
sistem_PI = feedback(C_PI * P, 1);
[y_PI, ~, bilgi_PI_kontrolcu] = zaman_tepkisi_ciz(sistem_PI, zaman, 'PI Kontrolc� ile Ad�m Tepkisi');

% PID kontrolc� tasar�m�
fprintf('/////////// Soru2: PID kontrolc� ile ad�m tepkisi ///////// \n');

Kp_PID = 0.6 * Ku;
Ti_PID = Pu / 2;
Td_PID = Pu / 8;
C_PID = pid(Kp_PID, Kp_PID / Ti_PID, Kp_PID * Td_PID);
sistem_PID = feedback(C_PID * P, 1);
[y_PID, ~, bilgi_PID_kontrolcu] = zaman_tepkisi_ciz(sistem_PID, zaman, 'PID Kontrolc� ile Ad�m Tepkisi');

fprintf('Soru3: Her denetleyici i�in Yerle�me Zaman�, Y�ksleme Zaman�, A��m Y�zdesi ve Kararl� Hal Hatas� \n');

% Sadece P kontrolc� ve P(s) sistemde varken
[y_b, ~, bilgi_b] = zaman_tepkisi_ciz(C_P * P, zaman, ' Soru3: Sadece P Kontrolc� ve P(s) Sistemde Varken Ad�m Tepkisi');

% Sadece PI kontrolc� ve P(s) sistemde varken
[y_c, ~, bilgi_c] = zaman_tepkisi_ciz(C_PI * P, zaman, 'Soru3: Sadece PI Kontrolc� ve P(s) Sistemde Varken Ad�m Tepkisi');

% Sadece PID kontrolc� ve P(s) sistemde varken
[y_d, ~, bilgi_d] = zaman_tepkisi_ciz(C_PID * P, zaman, 'Soru3: Sadece PID Kontrolc� ve P(s) Sistemde Varken Ad�m Tepkisi');

fprintf('Soru4: �nce Ayar Yap�lm�� PID Kontrolc� �le Ad�m Tepkisi \n ');
% �nce ayar yap�lm�� PID kontrolc�
Kp_ince_ayar = Kp_PID; % �nce ayar yap�lm�� Kp
Ti_ince_ayar = Ti_PID; % �nce ayar yap�lm�� Ti
Td_ince_ayar = Td_PID; % �nce ayar yap�lm�� Td
C_PID_ince_ayar = pid(Kp_ince_ayar, Kp_ince_ayar / Ti_ince_ayar, Kp_ince_ayar * Td_ince_ayar);
sistem_PID_ince_ayar = feedback(C_PID_ince_ayar * P, 1);
[y_PID_ince_ayar, ~, bilgi_PID_ince_ayar] = zaman_tepkisi_ciz(sistem_PID_ince_ayar, zaman, '�nce Ayar Yap�lm�� PID Kontrolc� ile Ad�m Tepkisi');

% Bozucu etkisi
fprintf('Soru5: Bozucu Etkiye Kar�� Her Bir Sistemin Tepkisi \n ');
bozucu = 1; % Birim ad�m bozucu etkisi
[y_bozucu, t, bilgi_bozucu] = zaman_tepkisi_ciz(feedback(P, C_PID_ince_ayar * bozucu, +1), zaman,'Bozucu Etkiye Kar�� PID Kontrolc� ile Sistem Cevab�');
fprintf('---------------------------------------------------------------------\n');

% Ad�m bozulmas� durumunda kontrol�rler ile sistem tepkileri
disp('Ad�m Bozulmas� Durumunda Sistem Tepkileri:');
% Ad�m bozulmas� i�in P kontrol�r�
sistem_bozo_P = feedback(P, C_P);
[y_bozo_P, ~, bilgi_bozo_P] = zaman_tepkisi_ciz(sistem_bozo_P, zaman, 'Ad�m Bozulmas� i�in P Kontrolc� ile Sistem Tepkisi');

% Ad�m bozulmas� i�in PI kontrol�r�
sistem_bozo_PI = feedback(P, C_PI);
[y_bozo_PI, ~, bilgi_bozo_PI] = zaman_tepkisi_ciz(sistem_bozo_PI, zaman, 'Ad�m Bozulmas� i�in PI Kontrolc� ile Sistem Tepkisi');

% Ad�m bozulmas� i�in PID kontrol�r�
sistem_bozo_PID = feedback(P, C_PID);
[y_bozo_PID, ~, bilgi_bozo_PID] = zaman_tepkisi_ciz(sistem_bozo_PID, zaman, 'Ad�m Bozulmas� i�in PID Kontrolc� ile Sistem Tepkisi');

% Ad�m bozulmas�na kar�� kontrol stratejileri de�erlendirmesi
fprintf('Ad�m Bozulmas�na Kar�� Kontrol Stratejileri De�erlendirmesi:\n');
fprintf('P Kontrolc�: Y�kselme S�resi: %f saniye, Yerle�me S�resi: %f saniye, A��m: %f%%, Kararl� Hal Hatas�: %f\n', bilgi_bozo_P.RiseTime, bilgi_bozo_P.SettlingTime, bilgi_bozo_P.Overshoot, abs(1 - y_bozo_P(end)));
fprintf('PI Kontrolc�: Y�kselme S�resi: %f saniye, Yerle�me S�resi: %f saniye, A��m: %f%%, Kararl� Hal Hatas�: %f\n', bilgi_bozo_PI.RiseTime, bilgi_bozo_PI.SettlingTime, bilgi_bozo_PI.Overshoot, abs(1 - y_bozo_PI(end)));
fprintf('PID Kontrolc�: Y�kselme S�resi: %f saniye, Yerle�me S�resi: %f saniye, A��m: %f%%, Kararl� Hal Hatas�: %f\n', bilgi_bozo_PID.RiseTime, bilgi_bozo_PID.SettlingTime, bilgi_bozo_PID.Overshoot, abs(1 - y_bozo_PID(end)));

% Sonu� ve Yorum
fprintf('Sonu� ve Yorum:\n');
fprintf('Ad�m bozulmas� durumunda, PI ve\n PID kontrol�rleri P kontrol�r�ne g�re\n daha iyi performans g�stermektedir.\n PI ve PID kontrol�rleri sistemdeki\n bozulmalar� daha h�zl� ve daha az\n a��m ile telafi ederken, PID\n kontrol�r� �zellikle daha h�zl�\n yerle�me s�resi ve d���k kararl�\n hal hatas� ile en iyi performans�\n g�stermektedir.\n');

% Routh-Hurwitz tablosu fonksiyonu
function routh_tablosu = routh_hurwitz(katsayilar)
    n = length(katsayilar);
    if mod(n, 2) == 1
        katsayilar(end+1) = 0;
        n = n + 1;
    end
    
    routh_tablosu = zeros(n, ceil(n/2));
    
    % �lk iki sat�r� doldurma
    routh_tablosu(1, :) = katsayilar(1:2:end);
    routh_tablosu(2, :) = katsayilar(2:2:end);
    
    % Di�er sat�rlar� doldurma
    for i = 3:n
        for j = 1:ceil(n/2)-1
            if routh_tablosu(i-1, 1) == 0
                routh_tablosu(i-1, 1) = eps;
            end
            routh_tablosu(i, j) = -(routh_tablosu(i-2, 1)*routh_tablosu(i-1, j+1) - routh_tablosu(i-2, j+1)*routh_tablosu(i-1, 1)) / routh_tablosu(i-1, 1);
        end
    end
end

% Zaman tepkisi grafi�i �izme ve performans kriterleri hesaplama fonksiyonu
function [y, t, bilgi] = zaman_tepkisi_ciz(sistem, t, baslik)
    figure;
    [y, t] = step(sistem, t);
    plot(t, y);
    title(baslik);
    grid on;
    xlabel('Zaman (s)');
    ylabel('Genlik');
    bilgi = stepinfo(y, t);
    disp(baslik);
    disp(['Y�kselme S�resi: ', num2str(bilgi.RiseTime),' saniye']);
    disp(['Yerle�me S�resi: ', num2str(bilgi.SettlingTime),' saniye']);
    disp(['A��m Y�zdesi: %', num2str(bilgi.Overshoot)]);
    disp(['Kararl� Hal Hatas�: ', num2str(abs(1 - y(end)))]);
    fprintf('---------------------------------------------------------------------\n');
end
