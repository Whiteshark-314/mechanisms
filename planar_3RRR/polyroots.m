function Z = polyroots(p)
% *** Solve multiple-root polymonials ***
%     F C Chang  Updated  12/22/2010
    mz = length(p)-max(find(p));
    p0 = p(min(find(p)):max(find(p)));
 if length(p0) < 2, Z = [0,mz]; return, end;
    sr = abs(p0(end)/p0(1));
 if sr < 1,   p0 = p0(end:-1:1);     end;
    q0 = polyder(p0);
    g1 = p0/p0(1); 
    g2 = q0(1:max(find(q0)))/q0(1);
for k = 1:2*length(p0),
    z12 = zeros(1,length(g1)-length(g2));
    z21 = zeros(1,length(g2)-length(g1));
    g3 = [g2,z12]-[g1,z21];
    g3 = g3(min(find(abs(g3)>1.e-8)):max(find(abs(g3)>1.e-8))); 
 if norm(g3,inf)/norm(g2,inf) < 1.e-3,
                             break;  end;
 if length(g1) >= length(g2),  
    g1 = g2;                         end;
    g2 = g3/g3(1);
end
    g0 = g1;
    u0 = deconv(p0,g0);
    v0 = deconv(q0,g0);
    w0 = polyder(u0);
    z0 = roots(u0);
    m0 = polyval(v0,z0)./polyval(w0,z0);
 if sr < 1,  z0 = z0.^-1;            end;
    Z = [z0,round(abs(m0))];
 if mz > 0,  Z = [Z; 0,mz];          end;