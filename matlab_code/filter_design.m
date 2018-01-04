%Design filters (FIR and IIR)
%---------------
Fs=256;
Fn=Fs/2;
f1 = [12.95 13.05];
f2 = [16.9 17.1];
f3 = [20.9 21.1];
Wp1 = f1/Fn; Ws1 = [f1(1)-1 f1(2)+1]/Fn;
Wp2 = f2/Fn; Ws2 = [f2(1)-1 f2(2)+1]/Fn;
Wp3 = f3/Fn; Ws3 = [f3(1)-1 f3(2)+1]/Fn;
Ast = 10;
Ap = 3;
%Design FIR
d1 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Ws1(1),Wp1(1),Wp1(2),Ws1(2),Ast,Ap,Ast); %FIR
Hd1 = design(d1,'equiripple');
fvtool(Hd1)
%--------------------
d2 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Ws2(1),Wp2(1),Wp2(2),Ws2(2),Ast,Ap,Ast); %FIR
Hd2 = design(d2,'equiripple');
fvtool(Hd2)
%--------------------
d3 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Ws3(1),Wp3(1),Wp3(2),Ws3(2),Ast,Ap,Ast); %FIR
Hd3 = design(d3,'equiripple');
fvtool(Hd3)
%--------------------
Hd_fir = {Hd1, Hd2, Hd3};
%Design IIR
Rp = Ap; Rs = Ast;
[n1,Wn1] = buttord(Wp1,Ws1,Rp,Rs);
[n2,Wn2] = buttord(Wp2,Ws2,Rp,Rs);
[n3,Wn3] = buttord(Wp3,Ws3,Rp,Rs);
if mod(n1,2)>0 %the order must be even
    n1=n1+1; 
end
if mod(n2,2)>0 %the order must be even
    n2=n2+1; 
end
if mod(n3,2)>0 %the order must be even
    n3=n3+1; 
end
%----------------------
d1 = fdesign.bandpass('N,F3dB1,F3dB2',n1,Wn1(1)*Fn,Wn1(2)*Fn,Fs);
Hd1 = design(d1,'butter');
fvtool(Hd1)
%----------------------
d2 = fdesign.bandpass('N,F3dB1,F3dB2',n2,Wn2(1)*Fn,Wn2(2)*Fn,Fs);
Hd2 = design(d2,'butter');
fvtool(Hd2)
%----------------------
d3 = fdesign.bandpass('N,F3dB1,F3dB2',n3,Wn3(1)*Fn,Wn3(2)*Fn,Fs);
Hd3 = design(d3,'butter');
fvtool(Hd3)
%----------------------
Hd_iir = {Hd1, Hd2, Hd3};
%Write to file
save('filters','Hd_fir','Hd_iir');