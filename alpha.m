function a = alpha(t,par)
%alpha(t)is the probability per unit time that a gcp divides into 2 gcps
% t is time, measured in (postnatal) days
% Td = doubling time for gcp's, passed through the "par" structure
    Td=par.Td;
    a=(log(2)/Td)-beta(t,par);
end
