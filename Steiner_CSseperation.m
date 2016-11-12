
# Radar_Steiner for convective-stratiform seperation

dbz_3km = randn(50,50)*15;

 
dx_km = 2;
dy_km = 2;
n_x = 50;
n_y = 50;

% find convective centers

icn = fix(11./dx_km); %background radius in 2.5-km grid cells rather than 11 km (4 grid cells)
dbz_con1 = dbz_3km*0;
dbz_con2 = dbz_3km*0;
dbz_con3 = dbz_3km*0;
% ... step 1: take all cells where dbz > 40
where_high = dbz_3km >= 40.;
dbz_con1(where_high) = 1;

% ... step 2: take cells where dbz > locally-dependent threshold
for ix = 1:n_x
    for iy = 1:n_y
        if( dbz_con1(ix,iy) == 0 && dbz_3km(ix,iy) > 0. )  ; %anything with echo not identified in round 1
            %     ... find locally-dependent threshold value
            dbz_bg = dbz_con1*0.;
            ix_min = max(ix-icn,1);
            ix_max = min(ix+icn,n_x);
            iy_min = max(iy-icn,1);
            iy_max = min(iy+icn,n_y);
            for ix_bg = ix_min:ix_max
                for iy_bg = iy_min:iy_max
                    if( sqrt((dx_km*(ix-ix_bg))^2+(dy_km*(iy-iy_bg))^2) <= 11. )
                        dbz_bg(ix_bg,iy_bg) = 1;
                    end
                end
            end
            dbz_bg(ix,iy) = 0 ; % exclude center point itself
            dbz_bg_mean = 0.;
            
            where_bg_echo =  (dbz_bg > 0) & ( dbz_3km > -90 ); %take only non-zero values into consideration
            
            dbz_bg_mean = 10.*log10(mean(10.^(dbz_3km(where_bg_echo)/10.)));
            
            if( dbz_bg_mean <= 0. )
                del_dbz_thresh = 10.;
            elseif( dbz_bg_mean <= 42.43 )
                del_dbz_thresh = 10. - dbz_bg_mean^2/180.;
            else
                del_dbz_thresh = 0.
            end
            dbz_thresh = dbz_bg_mean + del_dbz_thresh;
            if( dbz_3km(ix,iy) > dbz_thresh );
                dbz_con2(ix,iy) = 1;
            end
        end
        
    end
end

% ... step 3: also take cells within locally-determined radius
for ix = 1:n_x
    for iy = 1:n_y
        if( dbz_con1(ix,iy) == 1 || dbz_con2(ix,iy) == 1 )
            %     ... find background value that will determine local radius
            dbz_bg = dbz_3km*0.;
            ix_min = max(ix-icn,1);
            ix_max = min(ix+icn,n_x);
            iy_min = max(iy-icn,1);
            iy_max = min(iy+icn,n_y);
            for ix_bg = ix_min:ix_max
                for iy_bg = iy_min:iy_max
                    if( sqrt((dx_km*(ix-ix_bg))^2+(dy_km*(iy-iy_bg))^2) <= 11. )
                        dbz_bg(ix_bg,iy_bg) = 1 ;
                    end
                end
            end
            dbz_bg(ix,iy) = 0 ; % exclude center point itself
            dbz_bg_mean = 0.;
            where_bg_echo = dbz_bg > 0 & dbz_3km > -90 ; % take only non-zero values into consideration
            dbz_bg_mean = 10.*log10(mean(10.^(dbz_3km(where_bg_echo)/10.)));
            
            rad_thresh = 0.;
            if( dbz_bg_mean >= 25. && dbz_bg_mean < 30. )
                rad_thresh = 2.;
            elseif( dbz_bg_mean >= 30. && dbz_bg_mean < 35. )
                rad_thresh = 3.;
            elseif( dbz_bg_mean >= 35. && dbz_bg_mean < 40. )
                rad_thresh = 4.;
            elseif( dbz_bg_mean >= 40. )
                rad_thresh = 5.;
            end
            %      ... find values within local radius
            if( rad_thresh > 0. )
                dbz_bg = dbz_3km*0.;
                ix_min = max(ix-icn,1);
                ix_max = min(ix+icn,n_x-1);
                iy_min = max(iy-icn,1);
                iy_max = min(iy+icn,n_y);
                for ix_bg = ix_min:ix_max
                    for iy_bg = iy_min:iy_max
                        if( sqrt((dx_km*(ix-ix_bg))^2+(dy_km*(iy-iy_bg))^2) <= rad_thresh )
                            dbz_bg(ix_bg,iy_bg) = 1;
                        end
                    end
                end
                dbz_bg(ix,iy) = 0 ; % exclude center point itself
                where_high =   dbz_bg > 0 & dbz_con1 == 0 & dbz_con2 == 0  ;
                dbz_con3(where_high) = 1 ;   % close is enough
            end
        end
    end
end
dbz_con = dbz_con1 + dbz_con2 + dbz_con3;

% identify stratiform areas
hi_strat = 0 ; % add additional requirement for stratiform area aloft
hi_km = 6    ; % km elevation required to have dBZ > hi_dbz
hi_dbz = 5   ; % dBZ threshold to apply


dbz_str = dbz_3km*0.;
if( hi_strat )
    where_ech_hi =   (dbz_3km > -90) & dbz_hi > hi_dbz  ;
    dbz_str(where_ech_hi) = 1 ;  %everywhere there is echo low AND sufficient high
else
    where_ech = dbz_3km > 0  ;
    dbz_str(where_ech) = 1 ;  % everywhere there is echo
end
dbz_str( dbz_con > 0 ) = 0 ;   %subtract those identified as convective



subplot(2,2,1)
pcolor(dbz_3km); colorbar;
caxis([0 40])
shading flat

subplot(2,2,3)
pcolor(dbz_con); colorbar
shading flat

subplot(2,2,4)
pcolor(dbz_str); colorbar
shading flat

