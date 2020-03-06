S = System2D;


delta_t = 0.01;
ue = zeros(12,1);
Gravity = [0 -9.82 0 ]';

figure
M = calcM(S); % merge massmatrix from link 1 & 2 because it's needed in order to make the links move independently
for i=1:S.n
    ue(6*i-5:6*i-3) = Gravity*S.L(i).m; %external forces, we only care about gravity
end
        
for t = [0:delta_t:50]
    
   S = calcPos2D(S); 
   Jc1 = [ -sin(S.q(1))*S.L(1).CoM(1), 0; % Jacobian contains positions of CoM , i think. 
            cos(S.q(1))*S.L(1).CoM(1), 0; % Contains relations of joint velocities and end effector velocities
                                0 , 0;    % how does end effector(CoM) change when q change?
                                0 , 0;    
                                0 , 0;     
                                1 , 0];
                                
    Jc2 = [ -sin(S.q(1))*S.L(1).Length(1)-sin(S.q(1)+S.q(2))*S.L(2).CoM(1), -sin(S.q(1)+S.q(2))*S.L(2).CoM(1);
            cos(S.q(1))*S.L(1).Length(1)+cos(S.q(1)+S.q(2))*S.L(2).CoM(1), cos(S.q(1)+S.q(2))*S.L(2).CoM(1);
                                0 , 0;
                                0 , 0;
                                0 , 0;
                                1 , 1];
                            
    Jc = [Jc1; Jc2]; % merge Jacobians, think it's used to calculate lin & angular Velocities for every joint
    
    un = zeros(12,1); % wobbly effect, self generated forces 
    
    H = Jc' * M * Jc % calculate Inertia matrix, how "hard" it is to rotate the body. 
    
    S = forward_recursion(S); %to get xsi_dot0 which is lin/angular accelerations
    
    xsi_dot0 = zeros(S.n,1); % Initialize the vector
   
    for i=1:S.n
        % Fill in xsi_dot0 
        xsi_dot0(6*i-5:6*i-3,1) = S.L(i).dv;
        xsi_dot0(6*i-2:6*i,1) = S.L(i).dw;
    end
    
    uc = M*xsi_dot0 + un; % F = ma, uc =  total internal force?,  M = mass(m), xsidot0 = accelerations(a), un = internal forces(in this case 0)   
    c = Jc'*uc;   % is c a force?
    te = Jc'*ue;  % is te a force? why? here we introduce gravity to the mix atleast (ue = force of gravity)
    
    S.dqq = pinv(H)*(te - c);    % a = F/m  pinv(H) = 1/massmatrix ? 
    
    % We got the acceleration - integrate it to get the new dq and q(!)
    S.dq = S.dq + S.dqq * delta_t;
    S.q = S.q + S.dq * delta_t;
    axis([-5 5 -5 5]);
    Draw_System2D(S);
    
    pause(delta_t);
    
end



function M = calcM(S)

tempM = zeros(12,12);
j = 1;
for i =1:12
   if i < 7
       tempM(i,i) = S.L(1).massmatrix(i,i);
   else
       tempM(i,i) = S.L(2).massmatrix(j,j);
       j = j + 1;
   end 
end
M = tempM;
end


function S = calcPos2D(S)


for iL = 1:S.n % update all links positions (CoM) and rotation matrices

   
 % For the first link we only have to use joint 1 to compute the rotation matrix
 if (iL == 1)
    T_prev_link = eye(4);
 end

 %creating homogeneous matrices with rotation matrix & translation
 %vectors, transform_2d_func()
 % and multiplying in order to apply FGM and get positions
 
 T_joint_rotation = transform_2d_func([0 0 0]', S.q(iL));
 T_joint_to_CoM = transform_2d_func(S.L(iL).CoM, 0);
 T_CoM = T_prev_link * T_joint_rotation * T_joint_to_CoM;
 
 T_joint_to_next_link = transform_2d_func(S.L(iL).Length, 0);
 
 % Additional rotations are encoded into Roll,Pitch and Yaw angles. 
 % Add this final rotation matrix in the end.
%  R = rx(S.L(iL).RPY(1))*ry(S.L(iL).RPY(2))*rz(S.L(iL).RPY(3));
%  T_rotation = eye(4);
%  T_rotation(1:3, 1:3) = R;
  
 T_next_link = T_prev_link * T_joint_rotation * T_joint_to_next_link;% * T_rotation;
    
 S.L(iL).p = T_CoM(1:3,4);  % Extract the position
 S.L(iL).R = T_CoM(1:3, 1:3); % Extract the rotation matrix
 

 
 S.L(iL).transform_to_next_link = T_next_link;
 
 
 T_prev_link = T_next_link;
end
end


function T = transform_2d_func(t, th)

% Create an empty matrix of size 3x3.
T = zeros(4,4);

% Add translation values
T(1,4) = t(1);
T(2,4) = t(2);
% T(3,4) = t(3);

% Bottom right 1.
T(4,4) = 1;

% The rotation matrix. (around z)
R = [ cos(th)  -sin(th)   0 ;
      sin(th)   cos(th)   0 ;
           0        0   1 ];
 
% Copy the rotation matrix into T.
T(1:3,1:3) = R;
end


function S = forward_recursion(S)

%updates velocities & accelerations for every link according to previous values. 



  w_prev  = [0 0 0]'; % No rotational speed before the first joint (we'll asume that first joint is connected to something very sturdy)
  v_prev  = [0 0 0]'; % No linear speed before the first joint
  dw_prev = [0 0 0]'; % No rotational accelerations either before the first joint
  dv_prev = [0 0 0]'; % No linear accelerations either before the first joint
  r_prev  = [0 0 0]'; % Previous center of mass position (for the first iteration the w_prev will anyway be zero.
  
  
for iL = 1:S.n % iterate over all links

  k = [0 0 1]'; % joint axis of rotation/translation expressed in the "world" frame
  d = S.L(iL).R*S.L(iL).CoM; % vector from joint iL to CoM of link iL expressed in the "world"
                             % frame  (in case of a revolute joint)
  
  % Compute this once
  kXd = cross(k,d);
  
  %magic happens below 
  S.L(iL).w = w_prev + k*S.dq(iL);
  S.L(iL).v = v_prev + cross(w_prev, (S.L(iL).p-r_prev)) + kXd*S.dq(iL);

  S.L(iL).dw = dw_prev + cross(w_prev,k)*S.dq(iL);
  S.L(iL).dv = dv_prev + cross(dw_prev, (S.L(iL).p-r_prev)) + cross(w_prev, S.L(iL).v-v_prev) + cross(S.L(iL).w, kXd)*S.dq(iL);
    
 
  % Update all the 'previous' values.
  w_prev = S.L(iL).w;   %angular velocity
  v_prev = S.L(iL).v;   %linear velocity
  dw_prev = S.L(iL).dw; %angular acceleration
  dv_prev = S.L(iL).dv; %linear acceleration
  r_prev = S.L(iL).p;   %position of CoM in link i
  
end

end