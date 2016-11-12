%%%%%%%%% Script to calculate PDF of area occupied by 30dbz contours  %%%%%

TOL=15;
gamma = 0.2857;


%/data/usama/SHEAR_RUNS/zero_mean_fixed_utop_200_zero_v/utop_20_depth_0_zeromean

W = wrf_postprocessing('/data/usama/SHEAR_RUNS/zero_mean_fixed_utop_200_zero_v/utop_20_depth_1500_zeromean/wrfout_d01_0001-04-11_00:00:00',10); 
     
ZZ=W.Z;

Z=squeeze(mean(mean(ZZ,2),3));
Z=Z';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j= [20];
    kkk= [1500 3000 4500 6000 9000];
  for  i=kkk;

    
     qvp = nc_varget(['/data/usama/SHEAR_RUNS/zero_mean_fixed_utop_200_zero_v/utop_' num2str(j) '_depth_' num2str(i) '_zeromean/wrfout_d01_0001-04-11_00:00:00'],['QVAPOR']); 
     qra = nc_varget(['/data/usama/SHEAR_RUNS/zero_mean_fixed_utop_200_zero_v/utop_' num2str(j) '_depth_' num2str(i) '_zeromean/wrfout_d01_0001-04-11_00:00:00'],['QRAIN']); 
     qsn = nc_varget(['/data/usama/SHEAR_RUNS/zero_mean_fixed_utop_200_zero_v/utop_' num2str(j) '_depth_' num2str(i) '_zeromean/wrfout_d01_0001-04-11_00:00:00'],['QSNOW']);   
     qgr = nc_varget(['/data/usama/SHEAR_RUNS/zero_mean_fixed_utop_200_zero_v/utop_' num2str(j) '_depth_' num2str(i) '_zeromean/wrfout_d01_0001-04-11_00:00:00'],['QGRAUP']);  
     
     
     qvp = qvp(10:end,:,:,:);
     qra = qra(10:end,:,:,:);
     qsn = qsn(10:end,:,:,:);
     qgr = qgr(10:end,:,:,:);
     
     
     TT = nc_varget(['/data/usama/SHEAR_RUNS/zero_mean_fixed_utop_200_zero_v/utop_' num2str(j) '_depth_' num2str(i) '_zeromean/wrfout_d01_0001-04-11_00:00:00'],['T']); % perturbation potential temperature
     Theta = TT+300; % potential temperature
    
     
     P = nc_varget(['/data/usama/SHEAR_RUNS/zero_mean_fixed_utop_200_zero_v/utop_' num2str(j) '_depth_' num2str(i) '_zeromean/wrfout_d01_0001-04-11_00:00:00'],['P']);
     PB = nc_varget(['/data/usama/SHEAR_RUNS/zero_mean_fixed_utop_200_zero_v/utop_' num2str(j) '_depth_' num2str(i) '_zeromean/wrfout_d01_0001-04-11_00:00:00'],['PB']);

     Pressure = PB+P;
     prs= Pressure;
     
     prs = prs(10:end,:,:,:);
     
    % PM=squeeze(mean(mean(mean(Pressure,1),3),4))/100; % mean profile of pressure in mb
     
     tmk = Theta.*(Pressure/1e5).^gamma;
     
     tmk = tmk(10:end,:,:,:);
     
     ST=size(tmk);
     
     for t=1:ST(1)
         qvp_t = squeeze(qvp(t,:,:,:));
         qra_t = squeeze(qra(t,:,:,:));
         qsn_t = squeeze(qsn(t,:,:,:));
         qgr_t = squeeze(qgr(t,:,:,:));
         tmk_t = squeeze(tmk(t,:,:,:));
         prs_t = squeeze(prs(t,:,:,:));
         
        
     
         dbz(t,:,:,:) = dbzcalc(qvp_t,qra_t,qsn_t,qgr_t,tmk_t,prs_t,96,96,49,1,1,1,1);
     end
     
      clear qvp qra qsn qgr tmk prs
     
     for t=1:ST(1)
         dbz_z = squeeze(mean(dbz(t,1:13,:,:),2)); %11 % 0-2 km layer
         RL = zeros(size(dbz_z));
         RL(dbz_z>TOL) = 1;
         RLL = bwlabel(RL,4);
         stats = regionprops(RLL,'all');
         A(t)=sum(cat(1,stats.Area)); % total area of dbz >15 in the layer 0-2km
         Num(t)=length(stats); % total number of blubs 
         
         r(t,:,:)= dbz_z;
         rll(t,:,:)= RLL;
         
%          C=contourf(dbz_z,[30 30]);
%          [Area,Centroid,IN]=Contour2Area(C);
%          A(t)=sum(Area);
     end
     
     eval(['A_' num2str(j) '_depth_' num2str(i) '= A;'])
     eval(['Num_' num2str(j) '_depth_' num2str(i) '= Num;'])
     
     eval(['r_' num2str(j) '_depth_' num2str(i) '= r;'])
     eval(['rll_' num2str(j) '_depth_' num2str(i) '= rll;'])
     
     clear A Area dbz_z dbz r rll RL RLL stats Num
          
     
   end
end


%%%%%%%%%%%%%%%%%%%% Blobs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bin_N_1=0:max(Num_20_depth_1500)+5;
N1=hist(Num_20_depth_1500,bin_N_1); 
figure
subplot(5,2,1)
bar(bin_N_1,N1./numel(Num_20_depth_1500));
title('Shear Depth = 1500 m')
%xlabel('Number of blbbs')
%ylabel('Normalized frequency')
axis([-10 50 0 0.3])


bin_N_2=0:2:max(Num_20_depth_3000)+5;
N2=hist(Num_20_depth_3000,bin_N_2); 
subplot(5,2,3)
bar(bin_N_2,N2./numel(Num_20_depth_3000));
title('Shear Depth = 3000 m')
%xlabel('Number of blbbs')
%ylabel('Normalized frequency')
axis([-10 50 0 0.3])

bin_N_3=0:2:max(Num_20_depth_4500)+5;
N3=hist(Num_20_depth_4500,bin_N_3); 
subplot(5,2,5)
bar(bin_N_3,N3./numel(Num_20_depth_4500));
title('Shear Depth = 4500 m')
%xlabel('Number of blbbs')
ylabel('Normalized frequency')
axis([-10 50 0 0.3])

bin_N_4=0:2:max(Num_20_depth_6000)+5;
N4=hist(Num_20_depth_6000,bin_N_4); 
subplot(5,2,7)
bar(bin_N_4,N4./numel(Num_20_depth_6000));
title('Shear Depth = 6000 m')
%xlabel('Number of blbbs')
%ylabel('Normalized frequency')
axis([-10 50 0 0.3])

bin_N_5=0:2:max(Num_20_depth_9000)+5;
N5=hist(Num_20_depth_9000,bin_N_5); 
subplot(5,2,9)
bar(bin_N_5,N5./numel(Num_20_depth_9000));
title('Shear Depth = 9000 m')
xlabel('Number of blobs')
%ylabel('Normalized frequency')
axis([-10 50 0 0.3])


%%%%%%%%%%%%%%%%%%%%%%%  Area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bin_A_1 = 0:20:max(A_20_depth_1500)+10;
aa1 = hist(A_20_depth_1500,bin_A_1);
%figure
subplot(5,2,2)
bar(bin_A_1,aa1./numel(A_20_depth_1500));
title('Shear Depth = 1500 m')
%xlabel('Area of blobs')
%ylabel('Normalized frequency')
axis([-100 700 0 0.2])

bin_A_2 = 0:20:max(A_20_depth_3000)+10;
aa2 = hist(A_20_depth_3000,bin_A_2);
subplot(5,2,4)
bar(bin_A_2,aa2./numel(A_20_depth_3000));
title('Shear Depth = 3000 m')
%xlabel('Area of blobs')
%ylabel('Normalized frequency')
axis([-100 700 0 0.2])

bin_A_3 = 0:20:max(A_20_depth_4500)+10;
aa3 = hist(A_20_depth_4500,bin_A_3);
subplot(5,2,6)
bar(bin_A_3,aa3./numel(A_20_depth_4500));
title('Shear Depth = 4500 m')
%xlabel('Area of blobs')
ylabel('Normalized frequency')
axis([-100 700 0 0.2])

bin_A_4 = 0:20:max(A_20_depth_6000)+10;
aa4 = hist(A_20_depth_6000,bin_A_4);
subplot(5,2,8)
bar(bin_A_4,aa4./numel(A_20_depth_6000));
title('Shear Depth = 6000 m')
%xlabel('Area of blobs')
%ylabel('Normalized frequency')
axis([-100 700 0 0.2])

bin_A_5 = 0:20:max(A_20_depth_9000)+10;
aa5 = hist(A_20_depth_9000,bin_A_5);
subplot(5,2,10)
bar(bin_A_5,aa5./numel(A_20_depth_9000));
title('Shear Depth = 9000 m')
xlabel('Area of blobs')
%ylabel('Normalized frequency')
axis([-100 700 0 0.2])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
subplot(2,6,1)
pcolor(squeeze(r_40_depth_12000(30,:,:))); shading flat;
subplot(2,6,2)
pcolor(squeeze(r_40_depth_12000(31,:,:))); shading flat;
subplot(2,6,3)
pcolor(squeeze(r_40_depth_12000(32,:,:))); shading flat;
subplot(2,6,4)
pcolor(squeeze(r_40_depth_12000(33,:,:))); shading flat;
subplot(2,6,5)
pcolor(squeeze(r_40_depth_12000(34,:,:))); shading flat;
subplot(2,6,6)
pcolor(squeeze(r_40_depth_12000(35,:,:))); shading flat;
subplot(2,6,7)
pcolor(squeeze(rll_40_depth_12000(30,:,:))); shading flat; 
subplot(2,6,8)
pcolor(squeeze(rll_40_depth_12000(31,:,:))); shading flat;
subplot(2,6,9)
pcolor(squeeze(rll_40_depth_12000(32,:,:))); shading flat;
subplot(2,6,10)
pcolor(squeeze(rll_40_depth_12000(33,:,:))); shading flat;
subplot(2,6,11)
pcolor(squeeze(rll_40_depth_12000(34,:,:))); shading flat;
subplot(2,6,12)
pcolor(squeeze(rll_40_depth_12000(35,:,:))); shading flat;

% figure
% subplot(2,6,1)
% pcolor(squeeze(r(30,:,:))); shading flat;
% subplot(2,6,2)
% pcolor(squeeze(r(31,:,:))); shading flat;
% subplot(2,6,3)
% pcolor(squeeze(r(32,:,:))); shading flat;
% subplot(2,6,4)
% pcolor(squeeze(r(33,:,:))); shading flat;
% subplot(2,6,5)
% pcolor(squeeze(r(34,:,:))); shading flat;
% subplot(2,6,6)
% pcolor(squeeze(r(35,:,:))); shading flat;
% subplot(2,6,7)
% pcolor(squeeze(rll(30,:,:))); shading flat; 
% subplot(2,6,8)
% pcolor(squeeze(rll(31,:,:))); shading flat;
% subplot(2,6,9)
% pcolor(squeeze(rll(32,:,:))); shading flat;
% subplot(2,6,10)
% pcolor(squeeze(rll(33,:,:))); shading flat;
% subplot(2,6,11)
% pcolor(squeeze(rll(34,:,:))); shading flat;
% subplot(2,6,12)
% pcolor(squeeze(rll(35,:,:))); shading flat;


















% 
% bin1 = 0:20:max(A_0_depth_12000);
% figure
% subplot(3,2,1)
% hist(A_0_depth_12000,bin1,'linewidth',2)
% title('U top = 0','fontsize',12)
% 
% bin2 = 0:20:max(A_10_depth_12000);
% subplot(3,2,2)
% hist(A_10_depth_12000,bin2,'linewidth',2)
% title('U top = 10','fontsize',12)
% 
% bin3 = 0:20:max(A_20_depth_12000);
% subplot(3,2,3)
% hist(A_20_depth_12000,bin3,'linewidth',2)
% title('U top = 20','fontsize',12)
% 
% 
% bin4 = 0:20:max(A_30_depth_12000);
% subplot(3,2,4)
% hist(A_30_depth_12000,bin4,'linewidth',2)
% title('U top = 30','fontsize',12)
% 
% bin5 = 0:20:max(A_40_depth_12000);
% subplot(3,2,5)
% hist(A_40_depth_12000,bin5,'linewidth',2)
% title('U top =40','fontsize',12)
% 
% 
% 
% 
% 
% 
% figure
% subplot(3,2,1)
% hist(A_0_depth_12000,25,'linewidth',2)
% title('U top = 0','fontsize',12)
% 
% bin2 = 0:20:max(A_10_depth_12000);
% subplot(3,2,2)
% hist(A_10_depth_12000,25,'linewidth',2)
% title('U top = 10','fontsize',12)
% 
% bin3 = 0:20:max(A_20_depth_12000);
% subplot(3,2,3)
% hist(A_20_depth_12000,25,'linewidth',2)
% title('U top = 20','fontsize',12)
% 
% 
% bin4 = 0:20:max(A_30_depth_12000);
% subplot(3,2,4)
% hist(A_30_depth_12000,25,'linewidth',2)
% title('U top = 30','fontsize',12)
% 
% bin5 = 0:20:max(A_40_depth_12000);
% subplot(3,2,5)
% hist(A_40_depth_12000,25,'linewidth',2)
% title('U top =40','fontsize',12)
% 
% 
% figure
% subplot(3,2,1)
% hist(A_0_depth_12000,length(A_0_depth_12000),'linewidth',2)
% title('U top = 0','fontsize',12)
% 
% 
% subplot(3,2,2)
% hist(A_10_depth_12000,length(A_10_depth_12000),'linewidth',2)
% title('U top = 10','fontsize',12)
% 
% 
% subplot(3,2,3)
% hist(A_20_depth_12000,length(A_20_depth_12000),'linewidth',2)
% title('U top = 20','fontsize',12)
% 
% 
% 
% subplot(3,2,4)
% hist(A_30_depth_12000,length(A_30_depth_12000),'linewidth',2)
% title('U top = 30','fontsize',12)
% 
% 
% subplot(3,2,5)
% hist(A_40_depth_12000,length(A_40_depth_12000),'linewidth',2)
% title('U top =40','fontsize',12)
