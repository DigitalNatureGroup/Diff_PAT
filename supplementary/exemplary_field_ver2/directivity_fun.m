function [D] = directivity_fun(k,a, theta_n)

D =2.* besselj(1, k.*a.*sin(theta_n))./(k.*a.*sin(theta_n));

end

