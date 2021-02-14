function [ y,p,r ] = quat2euler( q )
%QUAT2EULER Convert a quaternion to euler angles
%   Convert a quaternion to euler angles

r = atan2(2*(q(:,1).*q(:,2)+q(:,3).*q(:,4)),1-2*(q(:,2).^2+q(:,3).^2));
p = asin(2*(q(:,1).*q(:,3)-q(:,2).*q(:,4)));
y = atan2(2*(q(:,1).*q(:,4)+q(:,2).*q(:,3)),1-2*(q(:,3).^2+q(:,4).^2));

end

