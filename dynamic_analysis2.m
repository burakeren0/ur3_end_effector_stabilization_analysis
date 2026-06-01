clear; clc;

%% 1. Sembolik Değişkenlerin Tanımlanması
syms th1 th2 th3 th4 th5 th6 th1_dot th2_dot th3_dot th4_dot th5_dot th6_dot real
q_dot = [th1_dot; th2_dot; th3_dot; th4_dot; th5_dot; th6_dot];
q = [th1; th2; th3; th4; th5; th6];

%% 2. Robot Parametreleri (UR3)
% DH Parametreleri: [a, d, alpha, theta]
dh_params = [
    0,          0.1519,   pi/2,   th1;             
    -0.24365,   0,        0,      th2;      
    -0.21325,   0,        0,      th3;             
    0,          0.11235,  pi/2,   th4;      
    0,          0.08535, -pi/2,   th5;             
    0,          0.0819,   0,      th6              
];

% Link Kütleleri
m = [1.98; 3.4445; 1.437; 0.871; 0.805; 0.261];

% Kütle Merkezleri (Local Frame)
Oi_list = [
    0,     -0.02,  0;
    0.13,   0,     0.1157;
    0.05,   0,     0.0238;
    0,      0,     0.01;
    0,      0,     0.01;
    0,      0,    -0.02
];

% Atalet (Inertia) Matrisleri (Local Frame - YAML'dan)
I_local = cell(6,1);
I_local{1} = diag([0.008093, 0.008093, 0.005625]);
I_local{2} = diag([0.02172, 0.02172, 0.009618]);
I_local{3} = diag([0.006546, 0.006546, 0.003543]);
I_local{4} = diag([0.001610, 0.001610, 0.002250]);
I_local{5} = diag([0.001572, 0.001572, 0.002250]);
I_local{6} = diag([0.000136, 0.000136, 0.000179]);

%% 3. İleri Kinematik ve Koordinat Dönüşümleri
T0_i = eye(4);

Org_all = sym(zeros(3,7));   % Orijinler (Eklem dönüş noktaları)
CoM_all = sym(zeros(3,6));   % Kütle Merkezleri
Z_all   = sym(zeros(3,7));   % Dönme eksenleri
R_all   = cell(6,1);         % Base frame'e göre dönüş matrisleri

Org_all(:,1) = [0;0;0];      
Z_all(:,1)   = [0;0;1];      

fprintf('İleri Kinematik çözülüyor...\n');
for i = 1:6
    a   = dh_params(i, 1);
    d   = dh_params(i, 2);
    alp = dh_params(i, 3);
    th  = dh_params(i, 4);
    
    Ai = [cos(th), -sin(th)*cos(alp),  sin(th)*sin(alp), a*cos(th);
          sin(th),  cos(th)*cos(alp), -cos(th)*sin(alp), a*sin(th);
          0,        sin(alp),          cos(alp),         d;
          0,        0,                 0,                1];
    
    T0_i = T0_i * Ai; 
    
    % 1. Eklem Orijinini Kaydet (T matrisinin 4. sütunu)
    Org_all(:, i+1) = T0_i(1:3, 4);
    
    % 2. Kütle Merkezini (CoM) Base Frame'e Taşı
    Oi_hom = [Oi_list(i, :) 1]';
    Ob_i = T0_i * Oi_hom;
    CoM_all(:, i) = Ob_i(1:3);
    
    % 3. Z eksenini ve Rotasyon Matrisini Kaydet
    Z_all(:, i+1) = T0_i(1:3, 3);
    R_all{i} = T0_i(1:3, 1:3);
end

%% 4. Jacobian Matrislerinin Hesaplanması
fprintf('Jacobian Matrisleri oluşturuluyor...\n');
J_all = cell(6,1); 
for i = 1:6
    Jv = sym(zeros(3,6)); 
    Jw = sym(zeros(3,6)); 
    
    for j = 1:i
        % Doğru Formül: Yarıçap = Kütle Merkezi (CoM) - Orijin (Org)
        Jv(:,j) = cross(Z_all(:,j), CoM_all(:,i) - Org_all(:,j));
        Jw(:,j) = Z_all(:,j);
    end
    
    J_all{i} = simplify([Jv; Jw]); 
end

%% 5. D(q) - Eylemsizlik (Inertia) Matrisi
fprintf('D(q) Eylemsizlik matrisi hesaplanıyor...\n');
D = sym(zeros(6,6));
for i = 1:6
    Jv_i = J_all{i}(1:3,:);   % Lineer 
    Jw_i = J_all{i}(4:6,:);   % Açısal 
    
    % Local Inertia'yı Base Frame'e döndür: I_base = R * I_local * R^T
    I_base = R_all{i} * I_local{i} * (R_all{i}.');
    
    % Kinetik Enerji = Öteleme (m*Jv'*Jv) + Dönme (Jw'*I_base*Jw)
    D = D + m(i) * (Jv_i.' * Jv_i) + (Jw_i.' * I_base * Jw_i);
end
D = simplify(D);
D = vpa(D, 4); % Dev kesirleri 4 anlamlı rakamlı ondalıklara çevirir
disp('D(q) Matrisi (Sadeleştirilmiş):');
disp(D(1:2, 1:2)); % Ekrana sığması için sadece 2x2'lik kısmını örnek gösteriyoruz

%% 6. G(q) - Yerçekimi (Gravity) Matrisi
fprintf('\ng(q) Yerçekimi matrisi hesaplanıyor...\n');
g_acc = 9.81;              
G_vec = [0; 0; -g_acc];    
P = sym(0); 

for i = 1:6
    % Potansiyel enerji: m * g * h
    P = P + m(i) * (-G_vec.' * CoM_all(:, i));
end

g_q = sym(zeros(6, 1));
for k = 1:6
    g_q(k) = diff(P, q(k));
end
g_q = simplify(g_q);
g_q = vpa(g_q, 4); % Ondalık sayı formatına çevir
disp('g(q) Vektörü (Sadeleştirilmiş):');
disp(g_q);

%% 7. C(q, q_dot) - Coriolis ve Merkezkaç Matrisi
fprintf('\nC(q, q_dot) matrisi hesaplanıyor... (LÜTFEN BEKLEYİN, BİRAZ SÜREBİLİR)\n');
n = 6;
C = sym(zeros(n));
for k = 1:n
    for j = 1:n
        for i = 1:n
            % Christoffel Sembolleri
            C(k,j) = C(k,j) + 0.5 * ...
                ( diff(D(k,j), q(i)) + ...
                  diff(D(k,i), q(j)) - ...
                  diff(D(i,j), q(k)) ) * q_dot(i);
        end
    end
end
C = simplify(C);
C = collect(C, q_dot);
C = vpa(C, 4); % Karmaşık kesirleri 4 haneli ondalıklara çevir
fprintf('C(q, q_dot) hesaplaması başarıyla tamamlandı!\n');

fprintf('\n=== TÜM DİNAMİK MATRİSLER BAŞARIYLA ÇIKARILDI ===\n');
%% 8. PYTHON İÇİN OTOMATİK KOD ÜRETİMİ (Code Generation)
fprintf('\nPython Numba kodu (robot_dynamics.py) üretiliyor...\n');

% Dosyayı oluştur ve aç
fid = fopen('robot_dynamics.py', 'w');

% Python başlıkları ve Numba dekoratörü
fprintf(fid, 'import numpy as np\n');
fprintf(fid, 'from numba import jit\n\n');
fprintf(fid, '@jit(nopython=True, cache=True)\n');
fprintf(fid, 'def robot_dynamics(q, dq):\n');

% Girdi vektörlerini parçalama (th1, th2... olarak)
fprintf(fid, '    th1, th2, th3, th4, th5, th6 = q\n');
fprintf(fid, '    th1_dot, th2_dot, th3_dot, th4_dot, th5_dot, th6_dot = dq\n\n');

% D(q) Matrisi
fprintf(fid, '    D = np.zeros((6, 6))\n');
for i = 1:6
    for j = 1:6
        if D(i,j) ~= 0
            % Sembolik ifadeyi karaktere çevir
            expr = char(vpa(D(i,j), 4)); 

            % MATLAB sözdizimini Python (NumPy) sözdizimine uyarla
            expr = strrep(expr, 'cos', 'np.cos');
            expr = strrep(expr, 'sin', 'np.sin');
            expr = strrep(expr, '^', '**'); % Üslü ifadeler için

            fprintf(fid, '    D[%d, %d] = %s\n', i-1, j-1, expr);
        end
    end
end
fprintf(fid, '\n');

% C(q, dq) Matrisi
fprintf(fid, '    C = np.zeros((6, 6))\n');
for i = 1:6
    for j = 1:6
        if C(i,j) ~= 0
            expr = char(vpa(C(i,j), 4));
            expr = strrep(expr, 'cos', 'np.cos');
            expr = strrep(expr, 'sin', 'np.sin');
            expr = strrep(expr, '^', '**');
            fprintf(fid, '    C[%d, %d] = %s\n', i-1, j-1, expr);
        end
    end
end
fprintf(fid, '\n');

% g(q) Vektörü
fprintf(fid, '    g = np.zeros(6)\n');
for i = 1:6
    if g_q(i) ~= 0
        expr = char(vpa(g_q(i), 4));
        expr = strrep(expr, 'cos', 'np.cos');
        expr = strrep(expr, 'sin', 'np.sin');
        expr = strrep(expr, '^', '**');
        fprintf(fid, '    g[%d] = %s\n', i-1, expr);
    end
end

% Return satırı
fprintf(fid, '\n    return D, C, g\n');

fclose(fid);
fprintf('HARİKA! "robot_dynamics.py" dosyası çalışma klasörünüze kaydedildi.\n');