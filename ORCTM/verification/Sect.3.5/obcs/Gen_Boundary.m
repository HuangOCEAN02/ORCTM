clear;clc;close all
ie_g=4000;je_g=5;ke=300;
ieee='b';prec='real*8';
load('INI_TS_LS.mat')
SDtime = [1:1:10];
%% Tide Components
T_m2period = 44714.1643;T_k1period = 86164.0906;
U_amp_m2 =  0.05;U_amp_k1 =  0.04;
U_gph_m2 = T_m2period/4;U_gph_k1 = T_k1period/4;
%% East
fid=fopen(['OB_East_u_Amp.bin'],'w',ieee); 
amp(1,1:je_g) = U_amp_m2;
amp(2,1:je_g) = U_amp_k1;
fwrite(fid,amp',prec); fclose(fid);
fid=fopen(['OB_East_u_Gph.bin'],'w',ieee); 
Gph(1,1:je_g) = U_gph_m2;
Gph(2,1:je_g) = U_gph_k1;
fwrite(fid,Gph',prec); fclose(fid);
clear amp Gph

fid=fopen(['OB_East_v_Amp.bin'],'w',ieee); 
amp(1,1:je_g+1) = 0;
amp(2,1:je_g+1) = 0;
fwrite(fid,amp',prec); fclose(fid);
fid=fopen(['OB_East_v_Gph.bin'],'w',ieee); 
Gph(1,1:je_g+1) = 0;
Gph(2,1:je_g+1) = 0;
fwrite(fid,Gph',prec); fclose(fid);
clear amp Gph

ttt = zeros(je_g,ke,size(SDtime,2));
sss = zeros(je_g,ke,size(SDtime,2));
for k=1:ke
    ttt(:,k,:)=t_final_element(k);
    sss(:,k,:)=s_final_element(k);
end
fidt=fopen(['OB_East_T.bin'],'w',ieee); 
fids=fopen(['OB_East_S.bin'],'w',ieee);
for t=1:size(SDtime,2)
    for k=1:ke
     fwrite(fidt,squeeze(ttt(:,k,t)),prec);
     fwrite(fids,squeeze(sss(:,k,t)),prec);
    end
end
fclose(fidt);
fclose(fids);
%% West
fid=fopen(['OB_West_u_Amp.bin'],'w',ieee); 
amp(1,1:je_g) = U_amp_m2;
amp(2,1:je_g) = U_amp_k1;
fwrite(fid,amp',prec); fclose(fid);
fid=fopen(['OB_West_u_Gph.bin'],'w',ieee); 
Gph(1,1:je_g) = T_m2period/4;
Gph(2,1:je_g) = T_k1period/4;
fwrite(fid,Gph',prec); fclose(fid);
clear amp Gph

fid=fopen(['OB_West_v_Amp.bin'],'w',ieee); 
amp(1,1:je_g+1) = 0;
amp(2,1:je_g+1) = 0;
fwrite(fid,amp',prec); fclose(fid);
fid=fopen(['OB_West_v_Gph.bin'],'w',ieee); 
Gph(1,1:je_g+1) = 0;
Gph(2,1:je_g+1) = 0;
fwrite(fid,Gph',prec); fclose(fid);
clear amp Gph

for k=1:ke
    ttt(:,k,:)=t_final_element(k);
    sss(:,k,:)=s_final_element(k);
end
fidt=fopen(['OB_West_T.bin'],'w',ieee); 
fids=fopen(['OB_West_S.bin'],'w',ieee);
for t=1:size(SDtime,2)
    for k=1:ke
     fwrite(fidt,squeeze(ttt(:,k,t)),prec);
     fwrite(fids,squeeze(sss(:,k,t)),prec);
    end
end
fclose(fidt);
fclose(fids);

