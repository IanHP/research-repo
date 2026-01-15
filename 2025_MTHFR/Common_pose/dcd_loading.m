clear; close all; clc 

addpath('matdcd-1.0')
%datalar yükleniyor

nmc=88
t_cry=readdcd('../../WT_v11/run1/init_coords/tog_step3_lig_align.dcd',1:nmc);
t_cry=t_cry(1,:);

%atom number of cofactor dcd file
t_tot1=readdcd('../../WT_v11/run1/init_coords/tog_step3_lig_align.dcd',1:nmc);
t_tot2=readdcd('../../WT_v11/run2/init_coords/tog_step3_lig_align.dcd',1:nmc);
t_tot3=readdcd('../../WT_v11/run3/init_coords/tog_step3_lig_align.dcd',1:nmc);
t_tot4=readdcd('../../G196F_v1/run1/init_coords/tog_step3_lig_align.dcd',1:nmc);
t_tot5=readdcd('../../G196F_v1/run2/init_coords/tog_step3_lig_align.dcd',1:nmc);
t_tot6=readdcd('../../G196F_v1/run3/init_coords/tog_step3_lig_align.dcd',1:nmc);


t_tot=[t_tot1; t_tot2; t_tot3; t_tot4; t_tot5; t_tot6]; 


% t_tot=[t_tot1; t_tot2; t_tot3;];
% Tot içine datayı ekledi
tot=t_tot;

% Datanın frame sayisini hesapladı.Tot mean aldı(Atom '.'n frame ortalamsi aldi.
% Datayı mean değerinden çıkardı (totm).
st=size(tot,1);
mt=mean(tot);
totm=0*tot;
for i=1:st
    totm(i,:)=tot(i,:)-mt;
end

% Amino asit kordinaat sayısını aldı. Bunun bir 0 matriksini oluşturdu.
% Covariance Matrix oluşturdu. 
sb=size(tot,2);
R=zeros(sb,sb);
for f=1:st
    R=R+totm(f,:)'*totm(f,:);
end
R=R/st;

% Eigen vektörü elde etti
[V,D]=eigs(R,sb);
RD=V;

% V'(Eigen vektörlerün transpozu) ile totm (mean değerinden çıkarlan datayı) çarptı.
drt=V'*totm';

% Kristal mean aldı ve kristal datasından çıkarıldı. Bu aslında kullanılmıyor.
crym=t_cry-mt;
drc=V'*crym';
%control  dot(crym,V(:,2))
save results_tog.mat