clear; clc;
syms th1 th2 th3 th4 th5 th6 real

% UR3 DH Parameters
dh_params = [
    0,          0.1519,   pi/2,   th1;             
    -0.24365,   0,        0,      th2;      
    -0.21325,   0,        0,      th3;             
    0,          0.11235,  pi/2,   th4;              
    0,          0.08535, -pi/2,   th5;             
    0,          0.0819,   0,      th6              
];

% Forward Kinematics
T0_6 = eye(4); 
fprintf('Calculating symbolic matrix...\n');
for i = 1:6
    a = dh_params(i, 1);
    d = dh_params(i, 2);
    alp = dh_params(i, 3);
    th = dh_params(i, 4); 
    
    A_i = [cos(th), -sin(th)*cos(alp),  sin(th)*sin(alp), a*cos(th);
           sin(th),  cos(th)*cos(alp), -cos(th)*sin(alp), a*sin(th);
           0,        sin(alp),          cos(alp),         d;
           0,        0,                 0,                1];
       
    T0_6 = T0_6 * A_i;
end

T0_6 = simplify(T0_6);

% Substitute numerical joint angles (rad)
th_vals = [0, -pi/4, -pi/12, -pi/4, pi/2, 0]; 
T_safe = double(subs(T0_6, [th1, th2, th3, th4, th5, th6], th_vals));

disp('--------------------------------------------------');
disp('T0_6 Matrix:');
disp(T_safe);

% Extract Position
X = T_safe(1,4);
Y = T_safe(2,4);
Z = T_safe(3,4);

% Extract Orientation (Euler ZYX)
eul_angles = tform2eul(T_safe, 'ZYX'); 
Rx = eul_angles(3);
Ry = eul_angles(2);
Rz = eul_angles(1);

x_Reference = [X; Y; Z; Rx; Ry; Rz];

disp('--------------------------------------------------');
disp('Euler Reference Vector [X, Y, Z, Rx, Ry, Rz]:');
fprintf('X  : %10.4f m\n', x_Reference(1));
fprintf('Y  : %10.4f m\n', x_Reference(2));
fprintf('Z  : %10.4f m\n', x_Reference(3));
fprintf('Rx : %10.4f rad\n', x_Reference(4));
fprintf('Ry : %10.4f rad\n', x_Reference(5));
fprintf('Rz : %10.4f rad\n', x_Reference(6));

% Extract Orientation (Quaternion - ROS2 format [x,y,z,w])
quat_matlab = tform2quat(T_safe); 
quat_ros = [quat_matlab(2), quat_matlab(3), quat_matlab(4), quat_matlab(1)];

pose_Reference = [X; Y; Z; quat_ros(1); quat_ros(2); quat_ros(3); quat_ros(4)]';

disp('--------------------------------------------------');
disp('Quaternion Reference Vector [X, Y, Z, qx, qy, qz, qw]:');
fprintf('X  : %10.4f m\n', pose_Reference(1));
fprintf('Y  : %10.4f m\n', pose_Reference(2));
fprintf('Z  : %10.4f m\n', pose_Reference(3));
fprintf('qx : %10.4f\n', pose_Reference(4));
fprintf('qy : %10.4f\n', pose_Reference(5));
fprintf('qz : %10.4f\n', pose_Reference(6));
fprintf('qw : %10.4f\n', pose_Reference(7));