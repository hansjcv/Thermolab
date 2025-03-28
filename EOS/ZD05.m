function ZD05  = ZD05(density,TK)
density = density/1e3; % convert to g/cm3;

% Zhang and Duan 2005 EOS for water
m = 18.01528;
ZD05_R = 83.14467;            %'Gas Constant in units of cm^3 bar/mol/K
ZD05_Vc = 55.9480373;        %'Critical volume in units of cm^3/mol
ZD05_Tc = 647.25;            %'Critical temperature in units of Kelvin

Vr = m ./ density / ZD05_Vc;
Tr = TK / ZD05_Tc;

B = 0.349824207 - 2.91046273 ./ (Tr .* Tr) + 2.00914688 ./ (Tr .* Tr .* Tr);
C = 0.112819964 + 0.748997714 ./ (Tr .* Tr) - 0.87320704 ./ (Tr .* Tr .* Tr);
D = 0.0170609505 - 0.0146355822 ./ (Tr .* Tr) + 0.0579768283 ./ (Tr .* Tr .* Tr);
E = -0.000841246372 + 0.00495186474 / (Tr .* Tr) - 0.00916248538 / (Tr .* Tr .* Tr);
f = -0.100358152 ./ Tr;
g = -0.00182674744 .* Tr;

delta = 1 + B ./ Vr + C ./ (Vr .* Vr) + D ./ Vr.^4 + E ./ Vr.^5 + (f ./ (Vr .* Vr) + g ./ Vr.^4) * exp(-0.0105999998 ./ (Vr .* Vr));


ZD05 = ZD05_R .* TK .* density .* delta / m / 1e3;% converted to kbar
end