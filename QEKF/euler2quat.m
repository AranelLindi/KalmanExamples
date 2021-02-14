function [ quaternion ] = euler2quat( yaw, pitch, roll )
%EULER2QUAT convert euler angles to a quaternion
%   convert euler angles to a quaternion

cdz = cos(yaw/2);
sdz = sin(yaw/2);
cdy = cos(pitch/2);
sdy = sin(pitch/2);
cdx = cos(roll/2);
sdx = sin(roll/2);

q0l = cdz*cdy*cdx; q0r = sdz*sdy*sdx;
q1l = cdz*cdy*sdx; q1r = sdz*sdy*cdx;
q2l = cdz*sdy*cdx; q2r = sdz*cdy*sdx;
q3l = sdz*cdy*cdx; q3r = cdz*sdy*sdx;

q0 = q0l + q0r;
q1 = q1l - q1r;
q2 = q2l + q2r;
q3 = q3l - q3r;

quaternion = [q0 q1 q2 q3];

end

