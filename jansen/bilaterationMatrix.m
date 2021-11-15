function Z=bilaterationMatrix(ij,jk,ki,sign)
    arguments
        ij double
        jk double
        ki double
        sign double = 1
    end
a=ij.^2-jk.^2+ki.^2;
b=sign*sqrt((ij.^2+jk.^2+ki.^2).^2-2*(ij.^4+jk.^4+ki.^4));
b(abs(imag(b))>0)=NaN;
Z=(1./(2*ij.^2)).*[a,-b;b,a];
end
