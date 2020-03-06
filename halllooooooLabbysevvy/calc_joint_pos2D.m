function p = calc_joint_pos2D(S)

p = zeros(2,S.n);

% First joint always at 0,0
p(:,1) = [0 0]';

for i = 1:S.n-1

    
    % Joint 1 is connected with link 0 (the world) and link 1 (the first
    % moving link.
    p(:,i+1) = S.L(i).transform_to_next_link(1:2, 4);
  
end

%%%EOF