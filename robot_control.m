load('robot_dinamik_model.mat')
matlabFunction(D_simplifed, 'File', 'calc_D', 'Vars', {q});
matlabFunction(C, 'File', 'calc_C', 'Vars', {q, q_dot});
matlabFunction(g_q_vpa, 'File', 'calc_g', 'Vars', {q});

Kp_degerleri = [15, 12.5, 12.5, 12.0, 5, 2];
Kp = diag(Kp_degerleri);


Kv_degerleri = 2 * sqrt(Kp_degerleri);
Kv = diag(Kv_degerleri);

Ki_degerleri = [0.000, 0.000, 0.0000, 0.0000, 0.0000, 0.00000];
Ki = diag(Ki_degerleri);