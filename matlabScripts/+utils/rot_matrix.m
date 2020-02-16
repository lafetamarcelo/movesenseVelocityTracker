function [rot] = rot_matrix(phi,theta,psi)
    %ROT_MATRIX Summary of this function goes here
    %   Detailed explanation goes here

    phi = phi * pi / 180;
    theta = theta * pi / 180;
    psi = psi * pi / 180;
    
    % Create the rotation matrices
    Rvi_v1 = [ cos(-psi),  sin(-psi),   0; ...
              -sin(-psi),  cos(-psi),   0; ...
                   0    ,       0   ,   1];
    
    Rv1_v2 = [ cos(-theta),  0, -sin(-theta); ...
                    0     ,  1,       0     ; ...
               sin(-theta),  0,  cos(-theta)];
    
    Rv2_vb = [ 1,      0    ,    0     ; ...
               0,  cos(-phi), sin(-phi); ...
               0, -sin(-phi), cos(-phi)];
    
    rot = Rvi_v1 * Rv1_v2 * Rv2_vb;
end

