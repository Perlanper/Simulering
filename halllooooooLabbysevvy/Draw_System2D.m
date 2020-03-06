function S = Draw_System2D(S)

pJ = calc_joint_pos2D(S);

cla
hold on
line_w = 2;

for iL = 1:S.n % iterate over all joints / links
  
  % Plots the position of joint iJ
%   if S.J(iL).type == 'R'
    J(iL) = plot(pJ(1,iL),pJ(2,iL),'ro','MarkerSize',10);
%   else
%     J(iL) = plot(pJ(1,iL),pJ(2,iL),'rd','MarkerSize',7.5);
%   end
  
  % Plots a vector from joint iL to the CoM of link iL
  plot([pJ(1,iL) S.L(iL).p(1)], [pJ(2,iL) S.L(iL).p(2)],'b','LineWidth',line_w);

  % Plots a vector from link iL CoM to the end of link iL.
  end_of_link = S.L(iL).transform_to_next_link(1:2,4);
  plot([S.L(iL).p(1) end_of_link(1)], [S.L(iL).p(2) end_of_link(2)],'g','LineWidth',line_w);

  
end

% Plot a star at every CoM for all links
for iL = 1:S.n % iterate over all links
    plot(S.L(iL).p(1), S.L(iL).p(2), 'rx');
end

% if S.tailEnd > S.tailStart+50
%     S.tail(S.tailStart,:) = [];
%     S.tailStart = S.tailStart+1;
% end
% S.tailEnd = S.tailEnd+1;
% S.tail(S.tailEnd,:) = [S.L(S.n).transform_to_next_link(1,4), S.L(S.n).transform_to_next_link(2,4)];
% 
% for ta = S.tailStart:S.tailEnd
%     plot(S.tail(ta,1),S.tail(ta,2), '.');
% end


%%%EOF