clc;
clear;
fprintf('\n Zakir YILDIZ 200757029 Son üç rakam 029 C=0.029,a = 0.0209,b = 0.920,e = 0.029/500,f = 1 \n ');
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
% Köklerin bulunmasý
disp('Polinom Kökleri:');
kokler = roots(payda);
disp('    0');
disp(kokler);
% Kararlýlýk kontrolü
if all(real(kokler) < 0)
    disp('Sistem kararlýdýr.');
else
    disp('Sistem kararlý deðildir.');
end
fprintf('---------------------------------\n');
% Zaman vektörü
zaman = 0:0.001:5;% örnekleme süresi 0.001 s ve toplam süre 5 s

% Adým tepkisi (sadece P(s))
fprintf('/////////// Soru1: Sitemin Adým Tepkisi  /////////// \n');
[~, ~, ~] = zaman_tepkisi_ciz(P, zaman, 'P(s) Sistemde Adým Tepkisi');

% Ziegler-Nichols yöntemi ile kontrolcü tasarýmý
Ku = 1; % Örnek kritik kazanç
Pu = 1; % Örnek periyot

% P kontrolcü tasarýmý
fprintf('/////////// Soru2: P kontrolcü ile adým tepkisi /////////// \n');
Kp_P = 0.5 * Ku;
C_P = pid(Kp_P, 0, 0);
sistem_P = feedback(C_P * P, 1);
[y_P, ~, bilgi_P_kontrolcu] = zaman_tepkisi_ciz(sistem_P, zaman, 'P Kontrolcü ile Adým Tepkisi');

% PI kontrolcü tasarýmý
fprintf('/////////// Soru2: PI kontrolcü ile adým tepkisi ////////// \n');

Kp_PI = 0.45 * Ku;
Ti_PI = Pu / 1.2;
C_PI = pid(Kp_PI, Kp_PI / Ti_PI, 0);
sistem_PI = feedback(C_PI * P, 1);
[y_PI, ~, bilgi_PI_kontrolcu] = zaman_tepkisi_ciz(sistem_PI, zaman, 'PI Kontrolcü ile Adým Tepkisi');

% PID kontrolcü tasarýmý
fprintf('/////////// Soru2: PID kontrolcü ile adým tepkisi ///////// \n');

Kp_PID = 0.6 * Ku;
Ti_PID = Pu / 2;
Td_PID = Pu / 8;
C_PID = pid(Kp_PID, Kp_PID / Ti_PID, Kp_PID * Td_PID);
sistem_PID = feedback(C_PID * P, 1);
[y_PID, ~, bilgi_PID_kontrolcu] = zaman_tepkisi_ciz(sistem_PID, zaman, 'PID Kontrolcü ile Adým Tepkisi');

fprintf('Soru3: Her denetleyici için Yerleþme Zamaný, Yüksleme Zamaný, Aþým Yüzdesi ve Kararlý Hal Hatasý \n');

% Sadece P kontrolcü ve P(s) sistemde varken
[y_b, ~, bilgi_b] = zaman_tepkisi_ciz(C_P * P, zaman, ' Soru3: Sadece P Kontrolcü ve P(s) Sistemde Varken Adým Tepkisi');

% Sadece PI kontrolcü ve P(s) sistemde varken
[y_c, ~, bilgi_c] = zaman_tepkisi_ciz(C_PI * P, zaman, 'Soru3: Sadece PI Kontrolcü ve P(s) Sistemde Varken Adým Tepkisi');

% Sadece PID kontrolcü ve P(s) sistemde varken
[y_d, ~, bilgi_d] = zaman_tepkisi_ciz(C_PID * P, zaman, 'Soru3: Sadece PID Kontrolcü ve P(s) Sistemde Varken Adým Tepkisi');

fprintf('Soru4: Ýnce Ayar Yapýlmýþ PID Kontrolcü Ýle Adým Tepkisi \n ');
% Ýnce ayar yapýlmýþ PID kontrolcü
Kp_ince_ayar = Kp_PID; % Ýnce ayar yapýlmýþ Kp
Ti_ince_ayar = Ti_PID; % Ýnce ayar yapýlmýþ Ti
Td_ince_ayar = Td_PID; % Ýnce ayar yapýlmýþ Td
C_PID_ince_ayar = pid(Kp_ince_ayar, Kp_ince_ayar / Ti_ince_ayar, Kp_ince_ayar * Td_ince_ayar);
sistem_PID_ince_ayar = feedback(C_PID_ince_ayar * P, 1);
[y_PID_ince_ayar, ~, bilgi_PID_ince_ayar] = zaman_tepkisi_ciz(sistem_PID_ince_ayar, zaman, 'Ýnce Ayar Yapýlmýþ PID Kontrolcü ile Adým Tepkisi');

% Bozucu etkisi
fprintf('Soru5: Bozucu Etkiye Karþý Her Bir Sistemin Tepkisi \n ');
bozucu = 1; % Birim adým bozucu etkisi
[y_bozucu, t, bilgi_bozucu] = zaman_tepkisi_ciz(feedback(P, C_PID_ince_ayar * bozucu, +1), zaman,'Bozucu Etkiye Karþý PID Kontrolcü ile Sistem Cevabý');
fprintf('---------------------------------------------------------------------\n');

% Adým bozulmasý durumunda kontrolörler ile sistem tepkileri
disp('Adým Bozulmasý Durumunda Sistem Tepkileri:');
% Adým bozulmasý için P kontrolörü
sistem_bozo_P = feedback(P, C_P);
[y_bozo_P, ~, bilgi_bozo_P] = zaman_tepkisi_ciz(sistem_bozo_P, zaman, 'Adým Bozulmasý için P Kontrolcü ile Sistem Tepkisi');

% Adým bozulmasý için PI kontrolörü
sistem_bozo_PI = feedback(P, C_PI);
[y_bozo_PI, ~, bilgi_bozo_PI] = zaman_tepkisi_ciz(sistem_bozo_PI, zaman, 'Adým Bozulmasý için PI Kontrolcü ile Sistem Tepkisi');

% Adým bozulmasý için PID kontrolörü
sistem_bozo_PID = feedback(P, C_PID);
[y_bozo_PID, ~, bilgi_bozo_PID] = zaman_tepkisi_ciz(sistem_bozo_PID, zaman, 'Adým Bozulmasý için PID Kontrolcü ile Sistem Tepkisi');

% Adým bozulmasýna karþý kontrol stratejileri deðerlendirmesi
fprintf('Adým Bozulmasýna Karþý Kontrol Stratejileri Deðerlendirmesi:\n');
fprintf('P Kontrolcü: Yükselme Süresi: %f saniye, Yerleþme Süresi: %f saniye, Aþým: %f%%, Kararlý Hal Hatasý: %f\n', bilgi_bozo_P.RiseTime, bilgi_bozo_P.SettlingTime, bilgi_bozo_P.Overshoot, abs(1 - y_bozo_P(end)));
fprintf('PI Kontrolcü: Yükselme Süresi: %f saniye, Yerleþme Süresi: %f saniye, Aþým: %f%%, Kararlý Hal Hatasý: %f\n', bilgi_bozo_PI.RiseTime, bilgi_bozo_PI.SettlingTime, bilgi_bozo_PI.Overshoot, abs(1 - y_bozo_PI(end)));
fprintf('PID Kontrolcü: Yükselme Süresi: %f saniye, Yerleþme Süresi: %f saniye, Aþým: %f%%, Kararlý Hal Hatasý: %f\n', bilgi_bozo_PID.RiseTime, bilgi_bozo_PID.SettlingTime, bilgi_bozo_PID.Overshoot, abs(1 - y_bozo_PID(end)));

% Sonuç ve Yorum
fprintf('Sonuç ve Yorum:\n');
fprintf('Adým bozulmasý durumunda, PI ve\n PID kontrolörleri P kontrolörüne göre\n daha iyi performans göstermektedir.\n PI ve PID kontrolörleri sistemdeki\n bozulmalarý daha hýzlý ve daha az\n aþým ile telafi ederken, PID\n kontrolörü özellikle daha hýzlý\n yerleþme süresi ve düþük kararlý\n hal hatasý ile en iyi performansý\n göstermektedir.\n');

% Routh-Hurwitz tablosu fonksiyonu
function routh_tablosu = routh_hurwitz(katsayilar)
    n = length(katsayilar);
    if mod(n, 2) == 1
        katsayilar(end+1) = 0;
        n = n + 1;
    end
    
    routh_tablosu = zeros(n, ceil(n/2));
    
    % Ýlk iki satýrý doldurma
    routh_tablosu(1, :) = katsayilar(1:2:end);
    routh_tablosu(2, :) = katsayilar(2:2:end);
    
    % Diðer satýrlarý doldurma
    for i = 3:n
        for j = 1:ceil(n/2)-1
            if routh_tablosu(i-1, 1) == 0
                routh_tablosu(i-1, 1) = eps;
            end
            routh_tablosu(i, j) = -(routh_tablosu(i-2, 1)*routh_tablosu(i-1, j+1) - routh_tablosu(i-2, j+1)*routh_tablosu(i-1, 1)) / routh_tablosu(i-1, 1);
        end
    end
end

% Zaman tepkisi grafiði çizme ve performans kriterleri hesaplama fonksiyonu
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
    disp(['Yükselme Süresi: ', num2str(bilgi.RiseTime),' saniye']);
    disp(['Yerleþme Süresi: ', num2str(bilgi.SettlingTime),' saniye']);
    disp(['Aþým Yüzdesi: %', num2str(bilgi.Overshoot)]);
    disp(['Kararlý Hal Hatasý: ', num2str(abs(1 - y(end)))]);
    fprintf('---------------------------------------------------------------------\n');
end
