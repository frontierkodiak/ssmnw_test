
% --- Compare different func resps with same param values ---------------------

% Assume a two species system of one prey and one predator
N1 = 0.1:0.1:10; 
a = 1;
b = 1;
g = a/2;
d = b;

% NestProj FR1

fr1 = a .* N1;

% NestProj FR2 (same as FR2 in Okuyama & Holland 2008)

fr2 = a .* N1 ./ (b^-1 + N1);

% NestProj FR Unimodal

fru = a .* N1 ./ (b^-1 + N1) - g .* N1 ./ (d^-1 + N1);

% FR2 in Holland et al. 2002

fr2h = a .* N1 ./ (1 + a .* N1);

% FR Unimodal in Holland et al. 2002

fruh = a .* N1 ./ (1 + a .* N1) - g .* N1 ./ (1 + g .* N1);

% Plot results
figure()
plot(N1, fr1)
hold on
plot(N1,fr2)
plot(N1, fru)
plot(N1,fr2h, '--')
plot(N1, fruh, '--')
legend('fr1', 'fr2', 'fru', 'fr2h', 'fruh')



% ---- Compare different param value sets for unimodal resps -------------------

% Assume a two species system of one prey and one predator
N1 = 0.001:0.1:10; 
a = 1;
b = 1;
g = 1;
d = 1;

aVec = [0.1 1 10];
bVec = [0.1 1 10];
gVec = [0.1 1 10];
dVec = [0.1 1 10];

% NestProj FR2 (same as FR2 in Okuyama & Holland 2008)
fr2 = a .* N1 ./ (b^-1 + N1);

% NestProj FR Unimodal
fru = a .* N1 ./ (b^-1 + N1) - g .* N1 ./ (d^-1 + N1);

% NestProj FR Unimodal paramVec
fru_a = aVec' * N1 ./ repmat((b^-1 + N1),3,1,1) -  repmat(g .* N1 ./ (d^-1 + N1),3,1,1);
fru_b = repmat(a .* N1,3,1,1) ./ ((bVec.^-1)' + repmat(N1,3,1,1)) -  repmat(g .* N1 ./ (d^-1 + N1),3,1,1);
fru_g = repmat(a .* N1 ./ (b^-1 + N1),3,1,1) - gVec' .* N1 ./ repmat((d^-1 + N1),3,1,1);
fru_d = repmat(a .* N1 ./ (b^-1 + N1),3,1,1) - repmat(g .* N1,3,1,1) ./ ((dVec.^-1)' + repmat(N1,3,1,1));

% FR2 in Holland et al. 2002
fr2h = a .* N1 ./ (1 + a .* N1);

% FR Unimodal in Holland et al. 2002
fruh = a .* N1 ./ (1 + a .* N1) - g .* N1 ./ (1 + g .* N1);

% FR Unimodal in Holland et al. 2002 paramVec
fruh_a = aVec' * N1 ./ (1 + aVec' * N1) - repmat(g .* N1 ./ (1 + g .* N1),3,1,1);
fruh_g = repmat(a * N1 ./ (1 + a * N1),3,1,1) - gVec' .* N1 ./ (1 + gVec' .* N1);

% O&H style aVec
figure()
hold on
plot(N1, fru_a)
plot(N1, fr2)
plot(N1, fru, '--')
legend( 'fru_a01', 'fru_a1', 'fru_a10', 'fr2', 'fru',  'location', 'bestoutside')
title('O&H extended. Avec')

% O&H style bVec
figure()
hold on
plot(N1, fru_b)
plot(N1, fr2)
plot(N1, fru, '--')
legend( 'fru_b01', 'fru_b1', 'fru_b10', 'fr2','fru',  'location', 'bestoutside')
title('O&H extended. Bvec')

% O&H style gVec
figure()
hold on
plot(N1, fru_g)
plot(N1, fr2)
plot(N1, fru, '--')
legend( 'fru_g01', 'fru_g1', 'fru_g10', 'fr2','fru',  'location', 'bestoutside')
title('O&H extended. Gvec')

% O&H style dVec
figure()
hold on
plot(N1, fru_d)
plot(N1, fr2)
plot(N1, fru, '--')
legend( 'fru_d01', 'fru_d1', 'fru_d10', 'fr2','fru', 'location', 'bestoutside')
title('O&H extended. Dvec')

% H style aVec
figure()
hold on
plot(N1,fruh_a)
plot(N1, fr2h)
plot(N1, fruh,'--')
legend(  'fruh_a01', 'fruh_a1', 'fruh_a10','fr2h', 'fruh', 'location', 'bestoutside')
title('Holland. Avec')

% H style gVec
figure()
hold on
plot(N1,fruh_g)
plot(N1, fr2h)
plot(N1, fruh, '--')
legend( 'fruh_g01', 'fruh_g1', 'fruh_g10', 'fr2h', 'fruh', 'location', 'bestoutside')
title('Holland. Gvec')


% Comparisons
figure()
hold on
plot(N1, fr2)
plot(N1, fr2h, '--')
plot(N1, fru_b)
plot(N1, fruh_a,'--')
title('O&H Bvec (full) vs Holland Avec (dashed)')
legend( 'fr2 (O&H)','fr2h', 'fru_b01', 'fru_b1', 'fru_b10', 'fruh_a01', 'fruh_a1', 'fruh_a10', 'location', 'bestoutside')

figure()
hold on
plot(N1, fr2)
plot(N1, fr2h, '--')
plot(N1, fru_d)
plot(N1, fruh_g,'--')
title('O&H Dvec (full) vs Holland Gvec (dashed)')
legend( 'fr2 (O&H)','fr2h', 'fru_d01', 'fru_d1', 'fru_d10', 'fruh_g01', 'fruh_g1', 'fruh_g10', 'location', 'bestoutside')



% Play around with Holland style

% Assume a two species system of one prey and one predator
N1 = 0.1:0.1:10; 
a = 1;
g = 0;%0.05;%
T2 = 0;%1;
TU = 0;%1;% 0.7;

% NestProj FR1

fr1 = a .* N1;

% FR2 in Holland et al. 2002

fr2h = a .* N1 ./ (1 + a .* T2 .* N1);

% FR Unimodal in Holland et al. 2002

fruh = a .* N1 ./ (1 + a .* TU .* N1) - g .* N1 ./ (1 + g .* N1);


% Plot results
figure()
plot(N1, fr1)
ylim([0,1])
hold on
plot(N1,fr2h, '--')
plot(N1, fruh, ':')
legend('fr1', 'fr2h', 'fruh')
% legend('fr2h', 'fruh')
