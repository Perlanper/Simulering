function T = transform_2d_func(t, th)
% -------------------------------
% Compute a transformation matrix
% -------------------------------
%
% Input:
% ------
% t - 2d vector containing the translation
% th - scalar containing the rotation in radians
%
% Output:
% -------
% T - transformation matrix in 2d



% Form a transformation matrix as:
% | r11 r12 t(1) |
% | r21 r22 t(2) |
% | 0   0   1    |

% Create an empty matrix of size 3x3.
T = zeros(4,4);

% Add translation values
T(1,4) = t(1);
T(2,4) = t(2);
T(3,4) = t(3);

% Bottom right 1.
T(4,4) = 1;

% The rotation matrix. (around z)
R = [ cos(th)  -sin(th)   0 ;
      sin(th)   cos(th)   0 ;
           0        0   1 ];
 
% Copy the rotation matrix into T.
T(1:3,1:3) = R;
