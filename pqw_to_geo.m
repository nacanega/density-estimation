function [X] = pqw_to_geo(X, params)
i = params.i;
Om = params.Om;
w = params.w;
rvec = X(1:3);
vvec = X(4:6);
IP = [cos(deg2rad(Om))*cos(deg2rad(w)) - sin(deg2rad(Om))*sin(deg2rad(w))*cos(deg2rad(i)), -cos(deg2rad(Om))*sin(deg2rad(w))-sin(deg2rad(Om))*cos(deg2rad(i))*cos(deg2rad(w)), sin(deg2rad(Om))*sin(deg2rad(i));
      sin(deg2rad(Om))*cos(deg2rad(w)) + cos(deg2rad(Om))*cos(deg2rad(i))*sin(deg2rad(w)), -sin(deg2rad(Om))*sin(deg2rad(w))+cos(deg2rad(Om))*cos(deg2rad(i))*cos(deg2rad(w)), -cos(deg2rad(Om))*sin(deg2rad(i));
      sin(deg2rad(i))*sin(deg2rad(w)), sin(deg2rad(i))*cos(deg2rad(w)), cos(deg2rad(i))];
R = IP*rvec;
V = IP*vvec;
X = [R;V];
end
