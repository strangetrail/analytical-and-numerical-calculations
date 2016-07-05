function modes
  syms e_ro e_phi e_z ro phi z Vro Vphi Vz
  Hro(ro,phi,z) = sym('Hro(ro,phi,z)');
  Hphi(ro,phi,z) = sym('Hphi(ro,phi,z)');
  Hz(ro,phi,z) = sym('Hz(ro,phi,z)');
  rot(Vro,Vphi,Vz) = sym('rot(Vro,Vphi,Vz)');
% Vro(ro,phi,z) = sym('Vro(ro,phi,z)');
% Vphi(ro,phi,z) = sym('Vphi(ro,phi,z)');
% Vz(ro,phi,z) = sym('Vz(ro,phi,z)');
  rot(Vro,Vphi,Vz) = [(diff(Vz,phi)/ro-diff(Vphi,z)) (diff(Vro,z)-diff(Vz,ro)) (diff(ro*Vphi,ro)/ro-diff(Vro,phi))];
  expand(rot(Hro,Hphi,Hz),'ArithmeticOnly',false,'IgnoreAnalyticConstraints',true)
end