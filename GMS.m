   # calculate the gross moist stability in double periodic domain
   # i.e. no horizontal advection.
   
    WM=squeeze(mean(W(:,1:19),1));
    TM = squeeze(mean(mean(mean(T,1),3),4));
    QM = squeeze(mean(mean(mean(qv,1),3),4));
    rho_m=squeeze(mean(mean(mean(rho,1),3),4));


    MSE= (cp*TM + Lv*QM + G*Z)/cp;

    DSE= (cp*TM + G*Z)/cp;



    for ii=1:18

            dh_dz(ii) =  (MSE (ii+1) - QM (ii)) / (Z(ii+1) - Z(ii));
        end
            dh_dz (19) =  (MSE (19) - MSE(18)) / (Z(19) - Z(18));




       for ii=1:18

            ds_dz(ii) =  (DSE (ii+1) - DSE (ii)) / (Z(ii+1) - Z(ii));
       end
            ds_dz (19) =  (DSE (19) - DSE (18)) / (Z(19) - Z(18));



    Th= WM.*dh_dz.*rho_m(1:19);  % this is : (w* dh/dz * rho) rho for vertical integral
    Ts= WM.*ds_dz.*rho_m(1:19);


    for jj=1:18
        Th_V(jj) = (0.5)*( Th(jj+1) + Th(jj) )*(Z(jj+1)-Z(jj));
    end

        VI_Th = sum(Th_V);



    for jj=1:18
        Ts_V(jj) =  0.5*( Ts(jj+1) + Ts(jj) )*(Z(jj+1)-Z(jj));
    end


        VI_Ts = sum(Ts_V);




        M = VI_Th./VI_Ts;
