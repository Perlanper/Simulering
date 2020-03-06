function S = System2D

%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%
S.n = 2; % number of joints, very handy for doing interations such as for i=1:SP.n ...

% definition of joints
% ----------------------------------------
S.J(1).axis = 'z';   % Joint axis (rotation axis or prismatic axis)
S.J(1).type = 'R';   % Value to separate the link from beeing a "rotational" one - R to a prismatic one - P.

S.J(2).axis = 'z';
S.J(2).type = 'R';
% % % %
% S.J(3).axis = 'z';
% S.J(3).type = 'R';
% 
% S.J(4).axis = 'z';
% S.J(4).type = 'R';
% 
% S.J(5).axis = 'z';
% S.J(5).type = 'R';

% definition of links
% ----------------------------------------
S.L(1).m = 1;
S.L(1).I = eye(3);
S.L(1).CoM = [ 1.0 0 0]'; % Position of the CoM
S.L(1).Length = [ 2.5 0]'; % Position of next joint (size of link)

S.L(2).m = 1;
S.L(2).I = eye(3);
S.L(2).CoM = [ 1.0 0 0]'; % Position of the CoM
S.L(2).Length = [ 2.7 0]'; % Position of next joint (size of link)




%%%%%%%%%%%%%%%% Variables %%%%%%%%%%%%%%%%%%%

% Joints positions
S.q   = zeros(S.n,1);
%S.q(1) = 45;
% Joints velocities 
S.dq  = zeros(S.n,1);
% Joints acceleration
S.ddq = zeros(S.n,1);
% Joints torques
S.tau = zeros(S.n,1);



S.tailStart = 1;
S.tailEnd = 1;
S.tail = [];

% Iterate over all links and initiate the variables, add things here to 
% avoid the need to recompute everything all the time.
for iL = 1:S.n
    S.L(iL).R  = eye(2);   % rotation matrix of Link
    S.L(iL).p  = [ 0 0 ]'; % p position of CoM
    S.L(iL).v  = [ 0 0 ]'; %linear velocity
    S.L(iL).dv = [ 0 0 ]'; %linear acceleration
    S.L(iL).w  = [ 0 0 ]'; %angular velocity
    S.L(iL).dw = [ 0 0 ]'; %angular acceleration
    S.L(iL).T  = [ 0 0 ]'; 
    S.L(iL).F  = [ 0 0 ]';
    S.L(iL).u  = [ 0 0 ]';
    S.L(iL).transform_to_next_link = eye(3); % Note - in world coordinates(!)
    S.L(iL).massmatrix = [S.L(iL).m 0 0 0 0 0
                          0 S.L(iL).m 0 0 0 0
                          0 0 S.L(iL).m 0 0 0
                          0 0 0 1 0 0
                          0 0 0 0 1 0
                          0 0 0 0 0 1];
end

%%% EOF