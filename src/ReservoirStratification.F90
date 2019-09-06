
program ResStrat

!-----------------------------------------------------------------------
! !MODULE: ResStratMod
!
! Reservoir Stratification Module
! Developed by Wondmagegn Yigzaw at Wonders of Water Lab of HongYi Li (Montana State University)
!
!
  implicit none

	integer,parameter :: r8 = selected_real_kind( 16) ! 8 byte real
	character :: rsrv*22		!reservoir latlon'ed name
	character :: energy*6
	integer :: fr_ln		!total timestep length (forcings @ hourly scale)
	integer :: fl_ln		!total length of reservoir inflow and river temperature @ daily scale
	integer :: sw_ln		!total length of solar radiation data
	integer :: r_leng		!total number of reservoirs to be modeled
	integer :: rsri							! reservoir under simulation
	integer :: ti							! timestep
	integer :: j,i,w,k,ii,l,ww,io,m,n,nn							! indices
	integer :: imixx,jmax,jmin,jmix							! indices

!!!! Constants !!!!
	integer :: d_n = 50              		! Number of vertical depth descritization with d_n + 1 layer areas
	integer, parameter :: d_n_max = 30				! Maximum number of layers
	integer, parameter :: yr_max = 10				! Maximum number of years simulated
	integer :: d_nn = 1000              		! Number of vertical depth descritization to establish depth-area-volume relationship
	real(r8) :: grav   = 9.8062	      				! gravity constant (m/s2)
	real(r8) :: t_frz   = 273.15  				!freezing temperature (K)
	real(r8) :: rho_w = 1.e3     		   			! Water density  (kg/m3)
	real(r8) :: rho_a = 1.177     		   			! Air density  (kg/m3)
	real(r8) :: c_w = 4.188e3     		   			! Water specific heat capacity  (J/kg/k)
	real(r8) :: t_max = 277.      					! temperature of maximum water density (K)
	real(r8) :: st_bl = 5.67e-8             		! Stefan-Boltzmann constant ~ W/m^2/K^4
	integer :: dtime != 3600        				! time step (sec)
	integer :: s_dtime != 3600        				! number of sub hourly time step
	real(r8) :: F = 1.0      						! dimensionless factor for wind sheltering by riparian vegetation, Wu et al, 2012

!!!! Input variables !!!!!!!!!!!!!

	real(r8) :: d_ht           			! Dam hwight (m)
	real(r8) :: d_res           			! Reservoir depth, taken as 0.95*(dam height) in GRanD database, (m)
	real(r8) :: depth						! Reservoir mean depth, taken from in GRanD database, (m)
	real(r8) :: s_w        				! solar radiation (W/m**2)
	real(r8) :: es                      			! Satuated vapor pressure(mb)
	real(r8) :: ea                      			! Atmospheric vapor pressure(mb)
	real(r8) :: pbot                    			! Atmospheric pressure(pa)
	real(r8) :: rh                       			! Relative humidity (%)
	real(r8) :: tsa                       			! 2m air temperature (k)
	real(r8) :: evap      							! evaporation rate (mm/d)
	real(r8) :: M_L              				! Mean lake length (km)
	real(r8) :: M_W              				! Mean lake width (km)
	real(r8) :: C_a           				! Surface area ciefficient (-)
	real(r8) :: C_v           				! Volume coefficient (-)
	real(r8) :: inflow_re       			! reservoir inflow (m3/s)
	real(r8) :: outflow_re       				! reservoir outflow (m3/s)
	real(r8) :: temp_riv       				! inflow temperature (k)
	real(r8) :: in_re(yr_max*365,3)       			! reservoir inflow (m3/s),inflow temperature (k),and reservoir outflow (m3/s)
	real(r8) :: dsr(yr_max*365,1)       			! direct solar radiation (w/m2)
	real(r8) :: forc(yr_max*365*24,6)       			! reservoir inflow (m3/s),inflow temperature (k),and reservoir outflow (m3/s)
	real(r8) :: t_air            			! air temperature
	real(r8) :: coszen             				!cosine of solar zenith angle
	real(r8) :: u_2       						! wind speed at 2m(m/s)
	real(r8) :: gm_j         					! Geometry code: 1='Rhombus_wedge'; 2='Oval_bowl';&
												!3='Rectangular_prism';4='Rectangular_bowl';4='Elliptical_bowl';6='Rectangular_wedge';&
												!7='Oval_wedge';8='Rhombus_bowl';9='Triangular_wedge';10='Parabolic_bowl';11='Parabolic_wedge';12='Triangular_bowl'

!!!! Calculated !!!!
	real(r8) :: kl							! empirical coefficient for the turbulent exchange of water vapor (mm/d/hpa)
	real(r8) :: d_res_sub           			! Reservoir depth, taken as 0.95*(dam height) in GRanD database, (m)
	real(r8) :: df_eff(d_n_max)  			          	! Effective diffusivity (molecular + eddy) [m2/s]
	real(r8) :: rho_z(d_n_max)     		   				! Depth based water density  (kg/m3)
	real(r8) :: rho_r	     		   				! Density of inflow water  (kg/m3)
	real(r8) :: drhodz(d_n_max)                 			! d [rhow] /dz (kg/m**4)
	real(r8) :: d_zs(d_n_max)			           			! Depth at z from surface (m)
	real(r8) :: d_zsb(d_n_max)			           			! Depth at z from surface averaged over sub-timestep(m)
	real(r8) :: d_zf(d_n_max)			           			! Depth at z from surface reversed to top layer(m)
	real(r8) :: d_z(d_n_max)			           			! Depth at z from bottom (m)
	real(r8) :: d_s	= 0.60		           			! Surface layer depth (m)
	real(r8) :: m1(d_n_max),m2(d_n_max),m3(d_n_max)    			! used in tridiagonal matrix
	real(r8) :: AX(d_n_max),BX(d_n_max),CX(d_n_max),DX(d_n_max),FX(d_n_max)	! used in tridiagonal matrix
	real(r8) :: a(d_n_max)    							! "a" left  diagonal of tridiagonal matrix
    real(r8) :: b(d_n_max)    							! "b" diagonal column for tridiagonal matrix
    real(r8) :: c(d_n_max)    							! "c" right  diagonal tridiagonal matrix
	real(r8) :: r(d_n_max)    							! "c" right  diagonal tridiagonal matrix
    real(r8) :: fac_1(d_n_max)                       	! Factor for calculation of triadiagonal matrices elements
	real(r8) :: fac_2(d_n_max)                       	! Factor for calculation of triadiagonal matrices elements
	real(r8) :: ri                      		! Richardson number
	real(r8) :: a_d(d_n_max)			   				! Area at depth z (km2)
	real(r8) :: a_dn(d_n_max)			   				! Adjusted area at depth z (km2)
	real(r8) :: v_z(d_n_max)			   			! Reservoir volume at depth z (m^3)
	real(r8) :: m_zo(d_n_max)			   			! Reservoir beginning mass at depth z (kg)
	real(r8) :: m_intial			   			! Reservoir beginning mass  (kg)
	real(r8) :: m_cal			   				! Reservoir calculated mass (kg)
	real(r8) :: m_mod			   				! Reservoir modeled mass (kg)
	real(r8) :: den			   				!Density (kg/m3)
	real(r8) :: m_cal_sub			   				! Reservoir calculated mass per sub-timestep (kg)
	real(r8) :: m_znn(d_n_max)			   			! Reservoir ending mass at depth z averaged for sub-timestep(kg)
	real(r8) :: m_zno(d_n_max)			   			! Reservoir initial mass at depth z (kg)
	real(r8) :: m_zn(d_n_max)			   			! Reservoir ending mass at depth z (kg)
	real(r8) :: m_zn_old(d_n_max)			   			! Reservoir ending mass at depth z to be used in halving/merging layers (kg)
	real(r8) :: m_zn_sub(d_n_max)			   			! Reservoir ending mass at depth z (kg)
	real(r8) :: m_in(d_n_max)			   			! Reservoir inflow mass at depth z (kg)
	real(r8) :: m_in_sub(d_n_max)			   			! Reservoir inflow mass at sub-timestep depth z (kg)
	real(r8) :: m_ou(d_n_max)			   			! Reservoir outflow mass at depth z (kg)
	real(r8) :: dm_nt(d_n_max)			   			! net mass added at depth z (kg)
	real(r8) :: dm_in(d_n_max)			   			! mass added to depth z (kg)
	real(r8) :: dm_ou(d_n_max)			   			! initial mass removed from depth z (kg)
	real(r8) :: m_ou_sub(d_n_max)			   			! Reservoir outflow mass at depth z (kg)
	real(r8) :: m_ev,e_cal,e_mod			   			! initial evaporation mass (kg)
	real(r8) :: m_ev_sub			   			! Reservoir outflow mass at depth z (kg)
	real(r8) :: dm_z(d_n_max)			   			! Reservoir mass change at depth z (kg)
	real(r8) :: dm_z_sub(d_n_max)			   			! Reservoir mass change at depth z (kg)
	real(r8) :: ds_z(d_n_max)			   			! Reservoir volume change at depth z (kg)
	real(r8) :: ds_z_sub(d_n_max)			   			! Reservoir volume change at depth z (kg)
	real(r8) :: d_v(d_n_max)			   			! Reservoir volume change at layer (m^3)
	real(r8) :: dv_in(d_n_max)			   			! volume increment at layer due to inflow(m^3)
	real(r8) :: dt_in(d_n_max)			   			! temperature increment at layer due to inflow(k)
	real(r8) :: dv_ou(d_n_max)			   			! volume decrease at layer due to inflow(m3)
	real(r8) :: dt_ou(d_n_max)			   			! temperature decrease at layer due to inflow(k)
	real(r8) :: v_t                        		! Total storage (m^3)
	real(r8) :: v_evap                        		! Evaporated volume (m^3)
	real(r8) :: d_evap                        		! Evaporated depth (m)
	real(r8) :: delta_z                        		! depth change to calculate corresponding area/volume(m)
	real(r8) :: top_d                        		! top layer depth
	real(r8) :: v_zt(d_n_max)                        		! Total reservoir volume at depth z from surface(m3)
	real(r8) :: v_mix                        		! Total volume of mixed layer(m3)
	integer :: nmix                        		! Mixed layer(-)
	integer :: Hday,Jday                        		! Hour of the day,Julian day (-)
	real(r8) :: in_tk                        		! Inflow zone thickness (m)
	real(r8)   :: in_f,in_t,ou_f       			! Inflow(m^3/s), inflow temp(k) and outflow(m^3/s)
	real(r8) :: s_tin                        		! Initial total storage (m^3)
	real(r8) :: s_t                        		! Total storage at timestep t (m^3)
	real(r8) :: ds_t                        		! change in storage (m^3)
	real(r8) :: t_s            				! surface temperature (k) -hourly
	real(r8) :: t_sdy            			! surface temperature (k) - for daily conversion
	real(r8) :: t_z(d_n_max)            			! lake layer temperature
	real(r8) :: t_zr(d_n_max)            			! lake layer temperature rounded to three decimal places
	real(r8) :: t_zf(d_n_max)            			! lake layer temperature to be written to a file adjusted so that surface temp will be at layer 50
	real(r8) :: t_zs(d_n_max)            			! lake layer temperature averaged over sub-timestep
	real(r8) :: t_z_old(d_n_max)            			! previous time lake layer temperature
	real(r8) :: t_zsub(d_n_max)            			! sub-timestep lake layer temperature (k)
	real(r8) :: dd_in(1000)			           			! Initial layer depth(m)
	real(r8) :: dd_z(d_n_max)			           			! Layer depth(m)
	real(r8) :: dd_inc(d_n_max)			           			! Depth increment (m)
	real(r8) :: dd_max			           			! Maximum depth change (m)
	real(r8) :: ddz_min,ddz_max						! Minimum and maximum layer thickness limit
	real(r8) :: bv_f                      			! brunt-vaisala frequency (/s**2)
	real(r8) :: eta      							! light extinction coefficient
	real(r8) :: sh_net        				! net short wave radiation
	real(r8) :: bias						! bias correction for solar radiation
	real(r8) :: sh_mix        				! net short wave radiation in mixed layer
	real(r8) :: beta                          		! shortwave absorbtion factor
	real(r8) :: lw_abs  						! atmospheric longwave absorbtion (W/m^2)
	real(r8) :: lw_abr  						! atmospheric longwave absorbtion for sub-timestep (W/m^2)
	real(r8) :: phi_o  								! net surface radiation (W/m^2)
	real(r8) :: phi_z(d_n_max)  							! radiation absorbed by layer (W/m^2)
	real(r8) :: phi_zf(d_n_max)  							! radiation absorbed by layer saved to file so that surface is layer 50
	real(r8) :: phi_x(d_n_max)  							! radiation absorbed by mixed layer (W/m^2)
	real(r8) :: phi_zsub(d_n_max)  							! sub-timestep radiation absorbed by layer (W/m^2)
	real(r8) :: lw_ems       					! longwave emission
	real(r8) :: lt_heat   						! latent heat
	real(r8) :: le= 2.501e6   						! latent heat of vaporaization (kcal/kg)
	real(r8) :: sn_heat   					! sensible heat
	real(r8) :: emt  				            	! average emittance
	real(r8) :: alb_s                 		! surface albedo shortwave radiation
	real(r8) :: k_m 					!molecular diffusivity
	real(r8) :: enr_ol(d_n_max) 					!Starting layer energy
	real(r8) :: enr_nw(d_n_max) 					!Ending energy change
	real(r8) :: enr_in(d_n_max) 					!Layer Energy from inflow
	real(r8) :: enr_ou(d_n_max) 					!Layer energy from outflow
	real(r8) :: enr_0(d_n_max) 					!Initial inner energy
	real(r8) :: enr_1(d_n_max) 					!Inner energy after advection
	real(r8) :: enr_2(d_n_max) 					!Inner energy after stratification
	real(r8) :: enr_2c(d_n_max) 					!Inner energy after stratification and convective mixing
	real(r8) :: enr_0r(d_n_max),enr_inr(d_n_max),enr_our(d_n_max),enr_1r(d_n_max),enr_2r(d_n_max) ! Area rated energy (w/m^2)
	real(r8) :: enr_0sub(d_n_max),enr_insub(d_n_max),enr_ousub(d_n_max),enr_1sub(d_n_max),enr_2sub(d_n_max) ! sub-timestep energy (w)
	real(r8) :: enr_err 					!Energy balance error (W/m^2)
	real(r8) :: enr_phi 					!Enegry from surface net flux before conv mixing
	real(r8) :: enr_phic 					!Enegry from surface net flux after conv mixing
	real(r8) :: enr_err2ca					!Enegry error (w/m2)
	real(r8) :: sn_heatc					!Enegry from sensible heat flux after conv mixing
	real(r8) :: s_err 						!Water balance error (m^3)
	real(r8) :: enr_err1,enr_err2,enr_err2c !Energy error (w) before stratification, after triadiagonal solution, and after convective mixing
	real(r8) :: k_eff,k_ew,k_ec,k_ad(d_n_max)		! Effective,wind,convection,advective kinetic energy (kg.m^2/s^2)
	real(r8) :: c_d			!Drag coefficient
	real(r8) :: cfa,cfw!,cft
	real(r8) :: Fr(d_n_max)  					!Froude number squared and inverted for diffusion coeff. calculation
	real(r8) :: dis_ad(d_n_max)
	real(r8) :: th_en(d_n_max)			!Layer thermal energy (j/s)
	real(r8) :: dis_w
	real(r8) :: l_vel
	real(r8) :: s_vel
	real(r8) :: tau
	real(r8) :: th_ex,enr_bc,enr_ac
	real(r8) :: cntr(d_n_max),cntr1(d_n_max),cntr2(d_n_max)			!Counters to average sub-timestep temperature
	real(r8) :: q_adv(d_n_max) !Layer flow rate (m^3/s)
	real(r8) :: z_str(d_n_max)	!Z* for layer solar radiation calculation (m)
	real(r8) :: dd_z_old(d_n_max),d_v_old(d_n_max),dv_in_old(d_n_max),dv_ou_old(d_n_max)
	real(r8) :: m_diff,m_cor(d_n_max)
	real(r8),dimension(1001) :: d_zi,a_di,v_zti			!Initial depth, area, volume for reservoir geometry
	integer  :: d_n_n								!Adjusted layer number
	real(r8) :: d_zn(d_n_max),d_vn(d_n_max),v_ztn(d_n_max),t_zn(d_n_max)		!Adjusted layer depth, volume,and temperature
	real(r8) :: denmix,denst,tmix,sumvol,tsum,vlmxlw,vlmxtp
	real(r8) :: msmxlw,msmxtp,summas,enmxlw,enmxtp,sumenr
	integer :: mixlow,mixtop
	real(r8) :: sh_neta,lw_absa,lt_heata,lw_emsa,sn_heata,phi_oa,st_sub
	real(r8) :: ta,tb,tab,dv_nt(d_n_max),dd_za,dd_zb,d_va,dv_oua,delta_a,delta_v,e_a,e_b,e_ab
	real(r8) :: dv_inb,dv_ina,dv_oub,dv_ouab,dv_inab,dd_zab,d_vab,d_vb,m_ab,m_a,m_b,num_fac,num_fac1
	real(r8) :: t_zt(d_n_max),rho_zt(d_n_max),rho_ztm
	real(r8) :: zmix,cmz(d_n_max),f_ri,by_e,cmzmix,drhomx,t_zmix
	real(r8) :: s_err_sub,enr_err_sub			!Sub-timestep storage and energy error
	real(r8) :: V_err,Ar_err,V_df,A_df		!Volume error (%), area error(&), volume difference(mcm), and area difference (km2) (+ve difference means underestimation)
	real(r8) :: A_cf,V_cf			!Area and volume correcting factor for error from geometry estimation
	real(r8) :: lat,lon			!Latitude and longitude  of reservoir grid
	real(r8) :: DECL,DELTS,HSR,sina, SUNSET,SUNUP,T1,T2
	integer :: flag, countmax,countmin

!! BOP

open(99, file='/data/ inputs/reservoirs.txt',action='read')

r_leng = 0
do
  read (99,'(A)',end = 2000)rsrv
  r_leng=r_leng+1 !total number or reservoirs to be modeled
end do
2000 rewind(99)

do rsri = 1,r_leng

  	energy='depth_'
  	read(99,'(A)')rsrv
	read(rsrv,99225)lat			!Read latitude of reservoir grid
	if(len(rsrv)==22)read(rsrv,99226)lon			!Read longitude of reservoir grid
	if(len(rsrv)==23)read(rsrv,99227)lon
99225 format(f7.2)
99226 format(9x,f7.2)
99227 format(9x,f8.2)

	open(12, file='/data/inputs/geometry/'//rsrv//'',action='read') !read geometry data
	open(13, file='/data/inputs/flow/'//rsrv//'',action='read') !read inflow, river temperature, and outflow
	open(14, file='/data/inputs/forcing/'//rsrv//'',action='read') !read forcings
	open(15, file='/data/outputs/stratification/'//rsrv//'',action='write') !open stratification result file
	open(199, file='/data/outputs/depth/'//energy//rsrv//'',action='write') !open depth result file

!	Start main program process
	fl_ln = 0
	do
		read (13,*,end = 2400)inflow_re,temp_riv,outflow_re
		fl_ln=fl_ln+1
		in_re(fl_ln,1) = inflow_re
		in_re(fl_ln,2) = temp_riv
		in_re(fl_ln,3) = outflow_re
	end do
2400 rewind(13)

	fr_ln = 0
	do
		read (14,*,end = 22)coszen,lw_abs,s_w,rh,t_air,u_2!
		fr_ln=fr_ln+1
		forc(fr_ln,1)=coszen
		forc(fr_ln,2)=lw_abs
		forc(fr_ln,3)=s_w
		forc(fr_ln,4)=rh
		forc(fr_ln,5)=t_air
		forc(fr_ln,6)=u_2
	end do
22 rewind(14)

	read (12,*)gm_j,depth,d_ht,M_L,M_W,V_err,Ar_err,C_v,C_a,V_df,A_df,d_n ! geometry code,depth,height, length,width, vol. error, area error, vol. coeff., area coeff., vol. difference, area difference, no. of layers

	if (d_ht <= 0.0)d_ht = depth
	d_res=0.95*d_ht


!	Calculate layer depth
	if (M_W <= 0.0)M_W = 1.
	if (M_L <= 0.0)M_L = 1.

	! Uniform subsurface layer depth for initialization	and limit maximum layer thickness
	do j = 1,d_nn
		dd_in(j) = d_res/d_nn !bottom layers evenly descritized
	end do

	! Area and volume correcting factors for relative error as compared to GRanD
	if (A_df>=0) then
		A_cf = 1. + (abs(Ar_err)/100.)
	elseif(A_df<0) then
		A_cf = 1. - (abs(Ar_err)/100.)
	end if

	if (V_df>=0) then
		V_cf = 1. + (abs(V_err)/100.)
	elseif(V_df<0) then
		V_cf = 1. - (abs(V_err)/100.)
	end if
	!	Calculate reservoir geometry to establish depth-area-volume relationship

	call rgeom (M_W,M_L,gm_j,d_res,dd_in,d_nn,C_a,C_v,d_zi,a_di,v_zti)

	!***********************************
	do j = 1, d_nn+1
		a_di(j)  = A_cf*a_di(j) 	!Area corrected for error
		v_zti(j) = V_cf*v_zti(j) 	!Volume corrected for error
	end do
	A_cf = 1.
	V_cf = 1.
	!***********************************

!	Initialize layer thickness if not prescribed
	! if (d_res <=3.) then
		! d_n = 1
	! else if (d_res >3. .and. d_res <=5.) then
		! d_n = 3
	! else if (d_res >5. .and. d_res <=10.) then
		! d_n = 4
	! else if (d_res >10. .and. d_res <=50.) then
		! d_n = 5
	! else if (d_res >50. .and. d_res <=100.) then
		! d_n = 10
	! else if (d_res >100. .and. d_res <=150.) then
		! d_n = 15
	! else if (d_res >150. .and. d_res <=200.) then
		! d_n = 20
	! else if (d_res >200.) then
		! d_n = 25
	! end if
	dd_z(d_n) = d_s		!top layer depth kept constant
	do j = d_n,1,-1
		if (j == d_n .and. d_n == 1) then
			dd_z(j) = d_res
		else if ((d_n>1 .and.j < d_n).and.(d_res - dd_z(d_n))>0.) then
			dd_z(j) = (d_res - dd_z(d_n)) / (d_n - 1) !bottom layers evenly descritized
		end if
	end do

	if (d_n>1 .and. dd_z(d_n-1)<d_s)then !layer thickness too small
		d_n=int(d_res/d_s)+1
		do j=1,d_n
			dd_z(j) = 0.
		end do
	!	Reinitialize layer thickness
		do j = d_n,1,-1
			if (j == d_n) then
				dd_z(j) = d_s      !top layer depth = 0.6m
			else
				dd_z(j) = (d_res - dd_z(d_n))/(d_n - 1) !bottom layers evenly descritized
			end if
		end do
	end if

	!	Calculate maximum and minimum layer thickness
	if (d_n>=15) then
		ddz_min = 1.5
		if (2.5*dd_z(d_n-1)>ddz_min) then
			ddz_max = 5.*dd_z(d_n-1)
		else
			ddz_max = 2.*ddz_min
		end if
	else if(d_n>1 .and. d_n<15) then
		ddz_min = 1.5
		if (2.5*dd_z(d_n-1)>ddz_min) then
			ddz_max = 2.5*dd_z(d_n-1)
		else
			ddz_max = 2.*ddz_min
		end if
	else if(d_n==1) then
		ddz_min = 1.5
		ddz_max = 2.*d_res
	end if
	m_cal = 0.
	Jday = 0
	flag=0
	countmax=0
	countmin=0

!	Calculate sub-time step for numerical stability
	dtime = 60
	s_dtime = 3600/dtime			!	Sub-hourly timestep. Change 3600 if the forcing data is not hourly
! 	************************Start calculation for each timestep	****************************************************
	! print*,'dtime= ',dtime,'sec'
	do ti = 1,fr_ln
	!	Initialize
		sh_neta=0.
		lw_absa=0.
		lt_heata=0.
		lw_emsa=0.
		sn_heata=0.
		phi_oa=0.
		s_err_sub=0.
		enr_err_sub=0.
		m_ev_sub=0.
		dm_z_sub=0.
	    ds_z_sub=0.
		m_cal_sub= 0.
		d_res_sub=0.
		st_sub=0.
		do j = 1,d_n_max
			t_zsub(j) = 0.
			phi_zsub(j) = 0.
			m_zn_sub(j)=0.			   			! Reservoir ending mass at depth z (kg)
			m_in(j)=0.			   			! Reservoir outflow mass at depth z (kg)
			m_in_sub(j)=0.			   			! Reservoir outflow mass at depth z (kg)
			m_ou(j)=0.			   			! Reservoir outflow mass at depth z (kg)
			m_ou_sub(j)=0.			   			! Reservoir outflow mass at depth z (kg)
			cntr(j)=0
			cntr1(j)=0
			cntr2(j)=0
			d_zsb=0.
		end do

		coszen=forc(ti,1)
		lw_abs=forc(ti,2)
		s_w=forc(ti,3)
		rh=forc(ti,4)
		t_air=forc(ti,5)
		u_2=forc(ti,6)

	! 	Index to match daily inflow/outflow value to hourly forcing
		if (mod(ti,24) == 0) then
			io = int(ti/24)
		else
			io = int(ti/24) + 1
		end if

	! 	Calculate julian day
		Jday = mod(io,365)
		Hday = mod(ti,24)
		if (mod(io,365) == 0)Jday = 365
		if (mod(ti,24) == 0)Hday = 24

		in_f=in_re(io,1)
		in_t=in_re(io,2)
		ou_f=in_re(io,3)

		!***************************************************************************************************************
		do ww = 1,s_dtime 	!	Start calculation for each sub-timestep	************************************************

			d_n_n=d_n

			if (ti==1 .and. ww==1 ) then
				do j = 1,d_n_max
					m_zo(j)=0.			   			! Reservoir beginning mass at depth z (kg)
					m_zn(j)=0.			   			! Reservoir ending mass at depth z (kg)
				end do
			end if

			if (ti==1 .and. ww==1 ) then
		!	Calculate layer depth (minimum at bottom)
				d_z(1)=0.
				do j = 2, d_n_max
					if (j<=d_n+1) then
					d_z(j) = d_z(j-1) + dd_z(j-1)
					else
					d_z(j) = 0.
					end if
				end do

		! 	Assign layer area and and volume based on depth-area-volume relationship
				do j = 1, d_n_max
					a_d(j) 	= 0.
					v_zt(j) = 0.
					d_v(j) 	= 0.
					dd_z(j) = 0.
				end do

				a_d(1)=0.1
				v_zt(1)=0.1
				do i=2,d_n+1
					do j=2,d_nn+1
						if (d_z(i)>d_zi(j-1).and.d_z(i)<=d_zi(j))then
							delta_z = (d_z(i)-d_zi(j-1))/(d_zi(j)-d_zi(j-1))
							a_d(i) = delta_z*(a_di(j)-a_di(j-1))+a_di(j-1)
							v_zt(i) = delta_z*(v_zti(j)-v_zti(j-1))+v_zti(j-1)
						else if (d_z(i) > d_zi(d_nn+1) .and. i<= d_n+1)then
							delta_z = (d_z(i)-d_zi(d_nn))/(d_zi(d_nn+1)-d_zi(d_nn))
							a_d(i) = delta_z*(a_di(d_nn+1)-a_di(d_nn))+a_di(d_nn)
							v_zt(i) = delta_z*(v_zti(d_nn+1)-v_zti(d_nn))+v_zti(d_nn)
						end if
					end do
				end do

				!	Calculate layer volume(m^3)
				do j = 1, d_n
					d_v(j) = v_zt(j+1) - v_zt(j)
					dd_z(j) = d_z(j+1) - d_z(j)
				end do
			end if

	! 	Intitialize layer temperature and total storage
			if (ti == 1 .and. ww==1) then
				do j=1,d_n
					t_z(j) = t_air
					rho_z(j) = den(t_z(j))
					m_zo(j) = V_cf*d_v(j)*rho_z(j)
					m_zn(j) = m_zo(j)
				end do
				s_tin = v_zt(d_n+1)
				s_t = s_tin
				m_intial =sum(m_zo)
			else
				do j=1,d_n
					t_z(j) = t_z(j)
				end do
				s_t = v_zt(d_n+1)
			end if

	! 	Allocate layer temperature as old for assigning counter for averaging sub-timestep result
			do j=1,d_n
				t_z_old(j) = t_z(j)
			end do

	! 	Calculation of Surface properties like: albedo, light extinction coefficient, etc
			if (coszen > 0.0) then
				alb_s = 0.05 / (0.15 + coszen)
			else
				alb_s = 0.06
			end if

			eta = 1.1925*(d_res)**(-0.424) ! as used in Subin et al, 2011 (citing Hakanson, 1995) but modified for actual reservoir depth
			beta = 0.175 !
			bias = 1.00 !

	! 	Calculation of Surface fluxes and heat source
		!Net shortwave radiation (w/m^2)
			sh_net = max(bias*s_w*(1 - alb_s),0.)
			t_s=t_z(d_n)		!Initialize surface temperature
			lw_abr = (1. - 0.03)*lw_abs				!longwave radiation (w/m^2)
			es=4.596*exp(17.25*(t_s-273.15)/t_s)		! in mmHg
			ea = 0.01*rh*4.596*exp(17.25*(t_air-273.15)/t_air)									! in mmHg
			lw_ems = 0.97*st_bl*t_s**4    !as used in henderson-sellers, 1984 (w/m^2)
			sn_heat = 1.5701*U_2* (t_s - t_air) !sensible heat (w/m^2)

			! Evaporation calculated as in Wu et al, 2012
			kl = 0.211 + 0.103 * U_2* F
			evap  = max(kl*133.322368*(es - ea)/100.,0.)  ! in mm/d; ea and es converted from mmHg to hpa
			if (t_s > t_frz) then
				lt_heat = rho_w * evap* le/(86.4e6)	!latent heat (w/m^2)
			else
				lt_heat = 0.0
			end if
			if (in_f==0.)evap=0. 	!Avoid continuous water abstraction if there is no inflow. Assumption is precipitation accounts for evaporation
			! net surface heat flux
			phi_o = sh_net + lw_abr - lt_heat - lw_ems - sn_heat

	! 	Calculation of reservoir density at depth z and incoming flow
			do j = 1, d_n_max
				if (j<=d_n)  then
					rho_z(j) = den(t_z(j))
				else
					rho_z(j) = 0.
				end if
			end do

			rho_r = den(in_t)

			! 	Calculation of equivalent evaporated depth (m), and volume(m^3)
			d_evap = evap*dtime/(86.4e6)
			v_evap = max(d_evap*A_cf*a_d(d_n+1),0.)
			m_ev = max(v_evap*rho_z(d_n-1),0.)

	! 	Calculation of flow contibution due to inflow/outflow
			if(s_t <= 0.10*(s_tin + V_df*1e6).and. ou_f > in_f)flag=1
			if((d_res<5. .and. d_n<=3) .and. ou_f > in_f)flag=1
			if(s_t <= 0.10*(s_tin + V_df*1e6).and. ou_f > in_f)ou_f=in_f  ! Avoid extraction of water from reservoir beyond 20% volume of the total storage
			if(s_t >= (s_tin + V_df*1e6) .and. in_f>ou_f)ou_f=in_f  ! Avoid extra inflow of water if the reservoir storage exceeded total dam storage (taken to be initial storage)
			if(d_res>=d_ht .and. in_f>ou_f)ou_f=in_f ! Avoid reservoir level from exceeding dam height
			if((d_res<5. .and. d_n<=3) .and. ou_f > in_f)ou_f=in_f  ! Avoid extraction of water from shallow reservoirs beyond 5m depth for numerical stability until reservoir operation is calibrated
			if(flag==1)v_evap=0.
			if(flag==1)m_ev=0.

			call flowdist(d_n,in_f,in_t,ou_f,dtime,d_v,v_zt,dv_in,dv_ou,dm_in,th_en,c_w)

	!	Resize layer thickness and numbers based on inflow/outflow contribution
		!******************************************************************************
	! 	Calculate initial layer and total mass (kg)
			if (ti==1 .and. ww==1 ) then
				do j = 1,d_n
					m_zo(j) = V_cf*d_v(j)*rho_z(j)
				end do
				m_intial=sum(m_zo)
			else
				do j = 1, d_n_max
					if (j<=d_n)  then
						m_zo(j) = m_zn(j)
					else
						m_zo(j) = 0.
					end if
				end do
			end if

			do j = 1, d_n_max
				if (j<=d_n)  then
					dm_ou(j) = dv_ou(j)*rho_z(j)*dtime
				else
					dm_ou(j) = 0.
				end if
			end do

  999		do j = 1, d_n
				if (j<d_n-1 .and. d_n>1) then
					dm_nt(j)=dm_in(j)-dm_ou(j)
					dv_nt(j)=dv_in(j)-dv_ou(j)
				else if(j==d_n-1 .and. d_n>1) then
					dm_nt(j)=dm_in(j)-dm_ou(j)-m_ev
					dv_nt(j)=dv_in(j)-dv_ou(j)-v_evap
				else if(j==d_n .and. d_n>1) then
					dm_nt(j)=dm_in(j)-dm_ou(j)
					dv_nt(j)=dv_in(j)-dv_ou(j)
				else if(j==d_n .and. d_n==1) then
					dm_nt(j)=dm_in(j)-dm_ou(j)-m_ev
					dv_nt(j)=dv_in(j)-dv_ou(j)-v_evap
				end if
			end do

	! 	Calculate layer mass (kg) and energy (w)
			num_fac=1.e6
			num_fac1=1.e3
			do j = 1, d_n_max
				if (j<=d_n)  then
					m_zn(j) = m_zo(j)+dm_nt(j)
					fac_1(j) = V_cf*d_v(j)*rho_z(j)*c_w/dtime
					enr_0(j) = t_z(j)*fac_1(j)/num_fac
				else
					m_zn(j) = 0.
					enr_0(j) = 0.
				end if
			end do

			if (ti==1 .and. ww==1 ) then
				m_cal=sum(m_zo)+sum(dm_nt)
				m_mod=sum(m_zo)+sum(dm_nt)
			else
				m_cal=m_cal+sum(dm_nt)
				m_mod=sum(m_zo)+sum(dm_nt)
			end if

			d_z(1)=0.
			a_d(1)=0.1
			v_zt(1)=0.1
			do i=1,d_n 			! check layers for available volume to satisfy net outflow
				if(-dm_nt(i) > m_zo(i))then !current layer collapses, hence remaining volume taken from next upper layer

					m_zn(i)=0.
					if (i<d_n-1 .and. d_n>1) then
						m_zn(i+1)=m_zn(i+1) - (-dm_nt(i)-m_zo(i))
					elseif (i==d_n-1 .and. d_n>2) then	!sub-layer collapses, hence remaining mass taken from next lower layer
						m_zn(i-1)=m_zn(i-1) - (-dm_nt(i)-m_zo(i))
						m=i-1
						do k=m,d_n
							v_zt(k+1)=v_zt(k)+m_zn(k)/rho_z(k)
							do j=2,d_nn+1
								if (v_zt(k+1)>v_zti(j-1).and.v_zt(k+1)<=v_zti(j))then
									delta_z  = (d_zi(j)-d_zi(j-1))*(v_zt(k+1)-v_zti(j-1))/(v_zti(j)-v_zti(j-1))
									delta_a  = (a_di(j)-a_di(j-1))*(v_zt(k+1)-v_zti(j-1))/(v_zti(j)-v_zti(j-1))
									d_z(k+1) = d_zi(j-1) + delta_z
									a_d(k+1) = a_di(j-1) + delta_a
								else if (v_zt(k+1)>v_zti(d_nn+1))then
									delta_z  = (d_zi(d_nn+1)-d_zi(d_nn))*(v_zt(k+1)-v_zti(d_nn))&
									/(v_zti(d_nn+1)-v_zti(d_nn))
									delta_a  = (a_di(d_n+1) -a_di(d_n))*(v_zt(k+1)-v_zti(d_nn))&
									/(v_zti(d_nn+1)-v_zti(d_nn))
									d_z(k+1) = d_zi(d_nn) + delta_z
									a_d(k+1) = a_di(d_nn) + delta_a
								end if
							end do
							d_z(k)=d_z(k+1)
							dd_z(k)=d_z(k+1)-d_z(k)
						end do
					elseif ((i==d_n-1 .or. i==d_n) .and. d_n<=2) then	! Top layer collapses, skip outflow
						dm_ou(i)=0.

						go to 999
					end if
					v_zt(i+1)=v_zt(i)+m_zn(i)/rho_z(i)
					do j=2,d_nn+1
						if (v_zt(i+1)>v_zti(j-1).and.v_zt(i+1)<=v_zti(j))then
							delta_z  = (d_zi(j)-d_zi(j-1))*(v_zt(i+1)-v_zti(j-1))/(v_zti(j)-v_zti(j-1))
							delta_a  = (a_di(j)-a_di(j-1))*(v_zt(i+1)-v_zti(j-1))/(v_zti(j)-v_zti(j-1))
							d_z(i+1) = d_zi(j-1) + delta_z
							a_d(i+1) = a_di(j-1) + delta_a
						else if (v_zt(i+1)>v_zti(d_nn+1))then
							delta_z  = (d_zi(d_nn+1)-d_zi(d_nn))*(v_zt(i+1)-v_zti(d_nn))/(v_zti(d_nn+1)-v_zti(d_nn))
							delta_a  = (a_di(d_n+1) -a_di(d_n))*(v_zt(i+1)-v_zti(d_nn))/(v_zti(d_nn+1)-v_zti(d_nn))
							d_z(i+1) = d_zi(d_nn) + delta_z
							a_d(i+1) = a_di(d_nn) + delta_a
						end if
					end do
					d_z(i)=d_z(i+1)
					dd_z(i)=d_z(i+1)-d_z(i)
				else !enough volume, layers don't collapses

					v_zt(i+1)=v_zt(i)+m_zn(i)/rho_z(i)
					do j=2,d_nn+1
						if (v_zt(i+1)>v_zti(j-1).and.v_zt(i+1)<=v_zti(j))then
							delta_z  = (d_zi(j)-d_zi(j-1))*(v_zt(i+1)-v_zti(j-1))/(v_zti(j)-v_zti(j-1))
							delta_a  = (a_di(j)-a_di(j-1))*(v_zt(i+1)-v_zti(j-1))/(v_zti(j)-v_zti(j-1))
							d_z(i+1) = d_zi(j-1) + delta_z
							a_d(i+1) = a_di(j-1) + delta_a
						else if (v_zt(i+1)>v_zti(d_nn+1))then
							delta_z  = (d_zi(d_nn+1)-d_zi(d_nn))*(v_zt(i+1)-v_zti(d_nn))/(v_zti(d_nn+1)-v_zti(d_nn))
							delta_a  = (a_di(d_n+1) -a_di(d_n))*(v_zt(i+1)-v_zti(d_nn))/(v_zti(d_nn+1)-v_zti(d_nn))
							d_z(i+1) = d_zi(d_nn) + delta_z
							a_d(i+1) = a_di(d_nn) + delta_a
						end if
					end do
					dd_z(i)=d_z(i+1)-d_z(i)
				end if
			end do

		! 	Recalculate layer thickness and volume
			do i = 1,d_n
				d_v(i)=m_zn(i)/rho_z(i)
			end do

	!	Recalculare reservoir depth
			d_res=0.
			do i=1,d_n
			   d_res=d_res+dd_z(i)
			end do

! Check if layers are too small
   71   	do i=1,d_n-1
				if(dd_z(i)<ddz_min)go to 75
			end do

			go to 91

! 	Identify lower (small and to be merged) and upper (larger) layer
   75   	if(i<d_n-1) then
				m=i+1
			else
				m=i-1 !Avoid merging to top layer
			end if

			ta =t_z(i)
			tb =t_z(m)
			dd_za=dd_z(i)
			dd_zb=dd_z(m)
			m_a =m_zn(i)
			m_b =m_zn(m)
			e_a =enr_0(i)
			e_b =enr_0(m)
			d_va=d_v(i)
			d_vb=d_v(m)
			dv_oua=dv_ou(i)
			dv_oub=dv_ou(m)
			dv_ina=dv_in(i)
			dv_inb=dv_in(m)

	!	Merge layers
			if(d_va+d_vb > 0.)then
				t_z(i)=(d_va*ta+d_vb*tb)/(d_va+d_vb)
			else
			end if

! 	Adjust new layer volume, thickness, mass and inflow/outflow
			dd_z(i)=dd_za+dd_zb
			d_v(i)=d_va+d_vb
			m_zn(i) =m_a+m_b
			dv_ou(i)=dv_oua+dv_oub
			dv_in(i)=dv_ina+dv_inb
			enr_0(i)=e_a+e_b

!	Re-number layers before collapse
			do j=i,d_n_max-1
				if (j==i .and. i<(d_n-1)) then
					t_z(j)=t_z(i)
					dv_ou(j)=dv_ou(i)
					dv_in(j)=dv_in(i)
					dd_z(j)=dd_z(i)
					d_v(j)=d_v(i)
					m_zn(j) =m_zn(i)
					enr_0(j)=enr_0(i)
					d_z(j+1)=d_z(i+2)
					v_zt(j+1)=v_zt(i+2)
					a_d(j+1)=a_d(i+2)
				else if (j<d_n .and. i<(d_n-1)) then
					t_z(j)=t_z(j+1)
					dv_ou(j)=dv_ou(j+1)
					dv_in(j)=dv_in(j+1)
					dd_z(j)=dd_z(j+1)
					d_v(j)=d_v(j+1)
					m_zn(j) =m_zn(j+1)
					enr_0(j)=enr_0(j+1)
					d_z(j+1)=d_z(j+2)
					v_zt(j+1)=v_zt(j+2)
					a_d(j+1)=a_d(j+2)
				else if (j==i .and. i==(d_n-1)) then
					t_z(j-1)=t_z(i)
					dv_ou(j-1)=dv_ou(i)
					dv_in(j-1)=dv_in(i)
					dd_z(j-1)=dd_z(i)
					d_v(j-1)=d_v(i)
					m_zn(j-1) =m_zn(i)
					enr_0(j-1)=enr_0(i)
					d_z(j)=d_z(j+1)
					v_zt(j)=v_zt(j+1)
					a_d(j)=a_d(j+1)
					t_z(d_n-1)=t_z(d_n)
					dv_ou(d_n-1)=dv_ou(d_n)
					dv_in(d_n-1)=dv_in(d_n)
					dd_z(d_n-1)=dd_z(d_n)
					d_v(d_n-1)=d_v(d_n)
					m_zn(d_n-1) =m_zn(d_n)
					enr_0(d_n-1)=enr_0(d_n)
					d_z(d_n)=d_z(d_n+1)
					v_zt(d_n)=v_zt(d_n+1)
					a_d(d_n)=a_d(d_n+1)
				else if (j>=d_n) then
					t_z(j)=0.
					dv_ou(j)=0.
					dv_in(j)=0.
					dd_z(j)=0.
					d_v(j)=0.
					m_zn(j) =0.
					enr_0(j)=0.
					d_z(j+1)=0.
					v_zt(j+1)=0.
					a_d(j+1)=0.
				end if
			end do
			d_n=d_n-1

			if(d_n <1)go to 4099
			go to 71

! 	Check if layers are too big
   91   	do i=1,d_n
				if(dd_z(i) >= ddz_max)go to 95
			end do
			go to 99

!	Calculate layer geometric properties to be halved
   95   	dd_zab=0.5*dd_z(i)
			m_ab =0.5*m_zn(i)
			e_ab =0.5*enr_0(i)
			d_vab=d_v(i)
			tab=t_z(i)
			dv_ouab=dv_ou(i)
			dv_inab=dv_in(i)

	!	Re-number layers before dividing layer
			d_z(d_n+2)=d_z(d_n+1)
			v_zt(d_n+2)=v_zt(d_n+1)
			a_d(d_n+2)=a_d(d_n+1)
			m=i+1
			do j=m,d_n
				k=d_n-j+i+1
				t_z(k+1)=t_z(k)
				dv_ou(k+1)=dv_ou(k)
				dv_in(k+1)=dv_in(k)
				d_z(k+1)=d_z(k)
				a_d(k+1)=a_d(k)
				v_zt(k+1)=v_zt(k)
				d_v(k+1)=d_v(k)
				dd_z(k+1)=dd_z(k)
				m_zn(k+1) = m_zn(k)
				enr_0(k+1)= enr_0(k)
			end do

	!	Divide layer in half and calculate corresponding properties
			dd_z(i)=dd_zab
			dd_z(i+1)=dd_zab
			d_z(i+1)=d_z(i)+dd_z(i)
			do j=2,d_nn+1
				if (d_z(i+1)>d_zi(j-1).and.d_z(i+1)<=d_zi(j))then
					delta_z   = (d_z(i+1)-d_zi(j-1))/(d_zi(j)-d_zi(j-1))
					a_d(i+1)  = delta_z*(a_di(j)-a_di(j-1))+a_di(j-1)
					v_zt(i+1) = delta_z*(v_zti(j)-v_zti(j-1))+v_zti(j-1)
				else if (d_z(i+1)> d_zi(d_nn+1))then
					delta_z   = (d_z(i+1)-d_zi(d_nn))/(d_zi(d_nn+1)-d_zi(d_nn))
					a_d(i+1)  = delta_z*(a_di(d_nn+1)-a_di(d_nn))+a_di(d_nn)
					v_zt(i+1) = delta_z*(v_zti(d_nn+1)-v_zti(d_nn))+v_zti(d_nn)
				end if
			end do
			d_v(i+1)=d_vab-(v_zt(i+1)-v_zt(i))
			d_v(i)=v_zt(i+1)-v_zt(i)
			dv_ou(i+1)=d_v(i+1)*dv_ouab/(d_v(i)+d_v(i+1))
			dv_ou(i)=d_v(i)*dv_ouab/(d_v(i)+d_v(i+1))
			t_z(i)=tab
			t_z(i+1)=tab
			m_zn(i) = m_ab
			m_zn(i+1) = m_ab
			enr_0(i)= e_ab
			enr_0(i+1)= e_ab
			dv_in(i+1)=d_v(i+1)*dv_inab/(d_v(i)+d_v(i+1))
			dv_in(i)=d_v(i)*dv_inab/(d_v(i)+d_v(i+1))
			d_n=d_n+1
			go to 91

   99		continue
   ! 	Recalculate density after layer change
			do j = 1, d_n_max
				if (j<=d_n)  then
					rho_z(j) = den(t_z(j))
				else
					rho_z(j) = 0.
				end if
			end do

			s_t = v_zt(d_n+1)
			if(s_t <= 0.0050*(s_tin + V_df*1e6) .or. d_res<=1.)go to 4099
			if(d_n ==1 .and.abs(t_z(d_n))>=450. )go to 4099
			if(sum(m_zn) <= 0.)go to 4099
			if(d_n >= d_n_max)go to 4099

	! 	Calculate layer internal energy (w) due to inflow/outflow
			do j = 1, d_n_max
				if (j<=d_n)  then
					enr_in(j) = dv_in(j)*in_t*rho_r*c_w/num_fac		!Energy from inflow
					enr_ou(j) = dv_ou(j)*t_z(j)*rho_z(j)*c_w/num_fac	!Energy loss due to outflow
				else
					enr_in(j) = 0.
					enr_ou(j) = 0.
				end if
			end do

	! 	Calculate layer energy (w)
			do j = 1, d_n_max
				if (j<=d_n)  then
					enr_1(j) = enr_0(j)+ enr_in(j) - enr_ou(j)
				else
					enr_1(j) = 0.
				end if
			end do

	! 	Adjust layer temperature for mixing due to inflow/outflow
			do j = 1, d_n_max
				if (j<=d_n)  then
					t_z(j) = enr_1(j)*num_fac*dtime/(m_zn(j)*c_w)
				else
					t_z(j) = 0.
				end if
			end do

		!	Check energy balance (w) after advective mixing
			enr_err1 = (sum(enr_1) - (sum(enr_0)+ sum(enr_in) - sum(enr_ou)))*num_fac


	!******************************************************************************
	! 	Calculate solar energy absorbed at each layer
			do j=1,d_n_max
				phi_x(j) = 0.
				phi_z(j) = 0.
			end do

			if (sh_net > 0.) then
				k=0
  600 			top_d=d_res-d_z(d_n-k)
				if(top_d < 0.61) then
					k=k+1
					go to 600
				else
				v_mix=V_cf*(v_zt(d_n+1)-v_zt(d_n-k))

			!	Solar radiation energy absorbed at the mixed zone
					sh_mix=sh_net*A_cf*(a_d(d_n+1)-(1.-beta)*a_d(d_n-k))
					j=d_n-k
					do i=j,d_n
					   phi_z(i)=sh_mix*V_cf*d_v(i)/v_mix
					end do
					phi_x(d_n+1)=sh_net
					phi_x(d_n-k)=(1.-beta)*sh_net
					if(k>0) then
						do i=1,k
						   ii=d_n-i+1
						   phi_x(ii)=(A_cf*a_d(ii+1)*phi_x(ii+1)-phi_z(ii))/(A_cf*a_d(ii))
						end do
					end if

		!	Solar radiation energy absorbed at sub layers
					l=d_n-k-1
					do j=1,l
					   i=(d_n-k)-j+1
					   phi_x(i-1)=phi_x(i)*exp(-eta*(d_z(i)-d_z(i-1)))
					end do

		!	Solar radiation energy absorbed in each layer
					j=d_n-k-1
					do i=1,j
					   phi_z(i)=A_cf*(a_d(i+1)*phi_x(i+1)-a_d(i)*phi_x(i))
					end do
				end if
			else
				do j=1,d_n
					phi_z(j) = 0.
				end do
			end if

	! *********************************************************************************************************************************
	! 	Calculation of effective diffusion coefficient
			if (u_2 >= 15.) then
				c_d = 2.6e-3
			else
				c_d = 5.e-4*sqrt(u_2)
			end if
			tau = rho_a*c_d*u_2**2. ! Shear stress at surface
			s_vel = sqrt(tau/rho_w) ! Shear velocity at surface
			k_ew=tau*s_vel*A_cf*a_d(d_n+1)*dtime	! Wind driven kinetic energy at surface
			dis_w = k_ew/(rho_w*V_cf*v_zt(d_n+1)*dtime)  ! rate of dissipation-wind
			cfw = 1.e-02!
			cfa = 1.e-05!
			k_m = 0.57/(c_w*rho_w) 						  !molecular diffusivity
			df_eff(1)= 0.			!bottom interface
			df_eff(d_n+1)= 0.		!air interface
			do j = 2,d_n
				q_adv(j) = max((dv_in(j)+dv_ou(j)),0.)
				k_ad(j)=0.5*rho_w*q_adv(j)*dtime*(q_adv(j)/(M_W*dd_z(j)))**2. ! Advection driven kinetic energy
				dis_ad(j)= k_ad(j)/(rho_w*V_cf*v_zt(j)*dtime)	! rate of dissipation-inflow/outflow

			! Calculate Richardson number
				drhodz(j) = (rho_z(j-1)-rho_z(j))/0.5*(dd_z(j)+dd_z(j-1))
				bv_f = max((grav/rho_w)*drhodz(j),0.)
				if (s_vel <= 0.) ri = 0.
				ri = bv_f/((s_vel/(0.4*d_z(j)))**2.)

			! Calculate Froude number
				l_vel = q_adv(j)*M_L/(A_cf*a_d(j)*dd_z(j))
				if (q_adv(j) <= 0. .or. drhodz(j) <= 0.) Fr(j) = 0.
				Fr(j)= (grav*dd_z(j)*drhodz(j)/rho_w)/l_vel**2.

			! Calculate diffusion coefficients

				df_eff(j)=min(max(dtime**2.*((cfw*dis_w/(1+ri)) &
				+(0.5*cfa*(dis_ad(j)+dis_ad(j-1))/(1+Fr(j)))),k_m),5.56e-03)
			end do

		!*****************************************************************
		! Calculate matrix elements
			do j=1,d_n_max
				a(j)=0.
				b(j)=0.
				c(j)=0.
				r(j)=0.
				AX(j)=0.
				BX(j)=0.
				CX(j)=0.
				DX(j)=0.
				FX(j)=0.
			end do

			do j = 1,d_n
				if (j == 1 .and. d_n>1) then
					m1(j) = 2*dtime/(0.5*A_cf*(a_d(j)+a_d(j+1))*dd_z(j))
					m2(j) = m1(j)*A_cf*a_d(j+1)*df_eff(j+1)/(dd_z(j)+dd_z(j+1))
					m3(j) = 0.
					Fx(j) = dtime*phi_z(j)/(V_cf*d_v(j)*c_w*rho_z(j))
					a(j) = - (m2(j))
					b(j) = 1. + (m2(j) + m3(j))
					c(j) = 0.
					r(j) = t_z(j) + Fx(j) !+ dt_in(j) - dt_ou(j)! bottom boundary condition
				else if (j <= d_n-1 .and. d_n>2) then
					m1(j) = 2*dtime/(0.5*A_cf*(a_d(j)+a_d(j+1))*dd_z(j))
					m2(j) = m1(j)*A_cf*a_d(j+1)*df_eff(j+1)/(dd_z(j)+dd_z(j+1))
					m3(j) = m1(j)*A_cf*a_d(j)*df_eff(j)/(dd_z(j)+dd_z(j-1))
					Fx(j) = dtime*phi_z(j)/(V_cf*d_v(j)*c_w*rho_z(j))
					a(j) = - m2(j)
					b(j) = 1. + m2(j) + m3(j)
					c(j) = - m3(j)
					r(j) = t_z(j) + Fx(j)
				else if (j == d_n) then!top layer
					m1(j) = 2*dtime/(0.5*A_cf*(a_d(j)+a_d(j+1))*dd_z(j))
					m2(j) = 0.
					m3(j) = m1(j)*A_cf*a_d(j)*df_eff(j)/(dd_z(j)+dd_z(j-1))
					Fx(j) = dtime*((phi_o-sh_net)*A_cf*a_d(d_n+1)+phi_z(j))/(V_cf*d_v(d_n)*c_w*rho_z(j))
					a(j) = 0.
					b(j) = 1. + (m2(j) + m3(j))
					c(j) = - (m3(j))
					r(j) = t_z(j) + Fx(j) ! top boundary condition
				end if
			end do

			!	Solve for temperature
				call solve(a,b,c,r,t_z,d_n)
			! end if

			! 	Avoid Numerical instability for multiple reservoir runs
					do j = 1,d_n
            if (isnan(t_z(j)))
						write(*,*),'check reservoir data'
						go to 4099
					end do

	!***********************************************************************************************************************
	! !	Solve convective mixing
		! Recalculate layer density
			do j = 1,d_n
				rho_z(j) = den(t_z(j))
			end do
			enr_bc=0.
			enr_ac=0.

		! Check if instability exists
			k=1
	  501   continue
			if(k >= d_n) go to 560
			if(rho_z(k) < rho_z(k+1)) go to 510
			k=k+1
			go to 501
	  510   continue

		! Start mixing layers
			mixlow=k
			mixtop=mixlow
			sumvol=V_cf*d_v(mixtop)
			summas=m_zn(mixtop)
			sumenr=enr_2(mixtop)
			tsum=t_z(mixtop)*m_zn(mixtop)
	  520   continue
			mixtop=mixtop+1
			vlmxtp=V_cf*d_v(mixtop)
			msmxtp=m_zn(mixtop)
			enmxtp=enr_2(mixtop)
			sumvol=sumvol+vlmxtp
			summas=summas + msmxtp
			sumenr=sumenr + enmxtp
			tsum=tsum+t_z(mixtop)*msmxtp
			tmix=tsum/summas

		! Calculate density of mixed layer
			denmix = den(tmix)
			if(mixtop == d_n) go to 540
			if(denmix < rho_z(mixtop+1)) go to 520
	  540   continue
			if(mixlow <= 1) go to 550

		! Check if instability exists below mixed layer	and mix layers
			if(rho_z(mixlow-1) >= denmix) go to 550
			mixlow=mixlow-1

		! Calculate temperature of mixed layer
			vlmxlw=V_cf*d_v(mixlow)
			msmxlw=m_zn(mixlow)
			enmxlw=enr_2(mixlow)
			sumvol=sumvol+vlmxlw
			summas=summas + msmxlw
			sumenr=sumenr + enmxlw
			tsum=tsum+t_z(mixlow)*msmxlw
			tmix= tsum/summas

		! Calculate density of mixed layer
			denmix = den(tmix)
			go to 540

	  550   continue

		! Set new layer temperature and density
			do j=mixlow,mixtop
				rho_z(j)=denmix
				t_z(j)=tmix
			end do

		! Calculate sum of layer energy in the mixing layer before mixing
			do j=mixlow,mixtop
				enr_ac=enr_ac + t_z(j)*m_zn(j)*c_w/(dtime*num_fac)
			end do
			k=mixtop
			go to 501

	  560    continue

		!   Recalculate layer final internal energy (w) after convective mixing
			do j = 1, d_n_max
				if (j<=d_n)  then
					enr_2c(j) = t_z(j)*m_zn(j)*c_w/(dtime*num_fac)
				else
					enr_2c(j) = 0.
				end if
			end do

		!	Check energy balance (w) after convective mixing
			enr_err2c = sum(enr_2c) - (sum(enr_1) + enr_phi)
			enr_err2c = enr_err2c*num_fac!/a_d(d_n+1)

	!*******************************************************************************************************************
	!	Calculate count for sub-timestep averaging
			do j = 1,d_n_max
				if (t_z_old(j).ne.0. .and. t_z(j).ne.0.) then
					cntr1(j) = cntr1(j) + 1
					cntr2(j) = 0
				else
					cntr1(j) = 0
					cntr2(j) = cntr1(j) + 1
				end if
				cntr(j)= cntr1(j)+cntr2(j)
			end do

	! 	Cummulative sub-timestep fluxes
			sh_neta=sh_neta + sh_net!*a_d(d_n+1)
			lw_absa=lw_absa + lw_abr
			lt_heata=lt_heata + lt_heat
			lw_emsa=lw_emsa + lw_ems
			sn_heata=sn_heata + sn_heat
			phi_oa=phi_oa + phi_o!*a_d(d_n+1)

			m_in_sub(j)=0.
			m_ou_sub(j)=0.
			d_res_sub = d_res_sub + d_res
			do j=1,d_n_max
				if(j<=d_n) then
					phi_z(j) = phi_z(j)
				else
					phi_z(j) = 0.
				end if
			end do
		!	Calculate layer depth for profile plot(minimum at the top)
			if (d_n>=2) then
				d_zs(1)=d_res - 0.5*dd_z(2)
				do j = 2, d_n_max
					if (j<=d_n) then
						d_zs(j) = d_zs(j-1) - 0.5*dd_z(j-1)- 0.5*dd_z(j)
					else
						d_zs(j) = 0.
					end if
				end do
				d_zs(d_n)=0.
			else
				d_zs(1)=d_res
			end if

			! Sum sub-timestep variables
			do j = 1,d_n_max
				t_zsub(j)   = t_zsub(j) + t_z(j)
				! phi_zsub(j) = phi_zsub(j) + phi_z(j)
				d_zsb(j)    = d_zsb(j) + d_zs(j)
			end do
			st_sub=st_sub+s_t

		end do

		d_res = d_res_sub/s_dtime

		if (d_n_n == d_n) then
			do j = 1,d_n_max
				t_zs(j) = t_zsub(j)/cntr(j)
				! phi_z(j) = phi_zsub(j)/cntr(j)
				d_zsb(j)    = d_zsb(j)/cntr(j)
			end do
		else
			do j = 1,d_n_max
				if (j<=d_n)	then
					t_zs(j)  = t_zsub(j)/cntr(j)
					! phi_z(j) = phi_zsub(j)/cntr(j)
					d_zsb(j) = d_zsb(j)/cntr(j)
				else
					t_zs(j)	 = 0.
					! phi_z(j) = 0.
					d_zsb(j) = 0.
				end if
			end do
		end if



		do j = d_n_max,1,-1   ! Reverse layer order so that first layer is at the top
			k =d_n-d_n_max+j
				if (k >= 1) then
					t_zf(j) = t_zs(k)
					! phi_zf(j) = phi_z(k)
					d_zf(j) = d_zsb(k)
				else
					t_zf(j)   = 0.
					! phi_zf(j) = 0.
					d_zf(j)   = 0.
				end if
		end do


		! do j = 1,d_n_max
			! if (isnan(t_zf(j)))t_zf(j)=-999. ! Avoid nan in temp file
		! end do

		write(15,40225),t_zf!
		write(199,50225),ti,d_n,in_f,ou_f,d_res,d_zf !write averaged reservoir depth, inflow/outflow

 40225	format(50f9.1)
 50225	format(i6,1X,i2,1X,f9.2,1X,f9.2,1X,f9.3,1x,50f9.1)

	end do

 4099	continue

end do

end program ResStrat

subroutine solve(a,b,c,r,u,n) !MNLAKE

		implicit none
!	 a - left coefficient , b - middle coefficient, c - right coefficient
!	 r - right hand side (known), u - the answer (unknown), n - number of equations
		integer,parameter :: r8 = selected_real_kind( 16) ! 8 byte real
        integer,intent(in) :: n
        real(r8),dimension(n),intent(in) :: a,b,c,r
        real(r8),dimension(n),intent(out) :: u
        real(r8),dimension(n) :: bp,rp
        real(r8) :: m,tt
        integer i

! 	initialize c-prime and d-prime
		bp(1)=b(1)
		rp(1)=r(1)
		do i=2,n
			tt=a(i)/b(i-1)
			bp(i) = b(i)-c(i-1)*tt
			rp(i) = r(i)-rp(i-1)*tt
		end do
		bp(n)=b(n)
		! rp(n)=r(n)
! 	initialize u
         u(n) = rp(n)/bp(n)

! 	Back substitution
        do i = n-1, 1, -1
			u(i) = (rp(i)-c(i)*u(i+1))/bp(i)
        end do

 end subroutine solve

 subroutine rgeom (M_W,M_L,gm_j,d_res,dd_in,d_nn,C_a,C_v,d_zi,a_di,v_zti)

! Calculate reservoir layer average area (km2)
	implicit none
	integer,parameter :: r8 = selected_real_kind( 16) ! 8 byte real
    real(r8), intent(in)  :: M_W,M_L,gm_j,d_res,dd_in(1000),C_a,C_v!,V_cf,A_cf
	integer, intent(in)  :: d_nn        !
	real(r8),dimension(1001), intent(out) :: d_zi,a_di,v_zti
	real(r8) ::dd_zz(1000),a_dd(1001),a_zi(1001), ar_f = 1.0e6 ,C_aa(1001)            	! Factor to convert area to m^2
	real(r8) :: pi      = 3.14159265358979323846  	! pi
	real(r8) :: dv,da,dz
	integer :: j,k							! indices

	do j = 1,d_nn
		dd_zz(j) = dd_in(j)
	end do

	do j=d_nn+1,1,-1
		C_aa(j) = 1.*C_a!
	end do

!	Calculate depth area
!******************************************* Curved Lake Bottom ****************************************************************
	do j = 1, d_nn
		if (gm_j == 1.0) then
			a_dd(j) = C_aa(j)*M_L*M_W*(1-((dd_zz(j)*(j-1))/d_res)**2.)
		else if (gm_j == 3.0) then
			a_dd(j) = C_aa(j)*M_L*M_W*(1-((dd_zz(j)*(j-1))/d_res)**2.)*((d_res-(dd_zz(j)*(j-1)))/d_res)**0.5
		else if (gm_j == 5.0) then
			a_dd(j) = C_aa(j)*pi*0.25*M_L*M_W*(1-((dd_zz(j)*(j-1))/d_res)**2.)*((d_res-(dd_zz(j)*(j-1)))/d_res)**0.5
		else if (gm_j == 2.0) then
			a_dd(j) = C_aa(j)*M_L*M_W*(1-((dd_zz(j)*(j-1))/d_res)**2.)*(1-((dd_zz(j)*(j-1))/d_res))
		else if (gm_j == 4.0) then
			a_dd(j) = C_aa(j)*(2./3.)*M_L*M_W*(1-((dd_zz(j)*(j-1))/d_res)**2.)*(1-((dd_zz(j)*(j-1)))/d_res)
		end if
		a_dd(d_nn+1) = 0.1			!Bottom area given non-zero value
	end do
!***********************************************************************************************************

	! Reverse indexing so that bottom is 1 and top is d_n+1 and convert to m2
	do j = 1,d_nn+1
		k =d_nn+2-j
		a_di(k) = max(a_dd(j),1.)*ar_f
	end do
	a_di(1) = 0.1			!Bottom Area
	if(a_di(d_nn+1)<a_di(d_nn))a_di(d_nn+1)=a_di(d_nn)

! 	Calculate layer depth from bottom,area,and volume
	d_zi(1) = 0.
	do j = 2, d_nn+1
		d_zi(j) = d_zi(j-1) + dd_in(j-1)
	end do

! 	Calculate layer average area,and total volume from bottom
	v_zti(1) = 0.1
	do j = 2, d_nn+1
		a_zi(j-1) = max(0.5*(a_di(j)+a_di(j-1)),1.) !Area converted to m^2
		v_zti(j) = (v_zti(j-1) + C_v*a_zi(j-1)*dd_in(j-1))
		! v_zti(j) = v_zti(j-1) + C_v*a_zi(j-1)*dd_in(j-1)+V_df/d_nn
	end do
	a_zi(d_nn+1)= a_di(d_nn+1)

	! do j = 1,d_nn+1
		! print*,d_zi(j),a_di(j),v_zti(j)
	! end do
	! stop
  end subroutine rgeom

subroutine flowdist(d_n,in_f,in_t,ou_f,dtime,d_v,v_zt,dv_in,dv_ou,dm_in,th_en,c_w)
!*******************************************************************************************************
! 	Calculation inflow/outflow contribution adopted from CE-QUAL-R1 model
!*******************************************************************************************************
	implicit none
	integer,parameter :: r8 = selected_real_kind( 16) ! 8 byte real
	integer,parameter :: d_n_max = 30
	integer, intent(in)  :: d_n,dtime        !
	real(r8), intent(in)  :: in_f,in_t,ou_f,c_w,d_v(d_n_max),v_zt(d_n_max)
	real(r8),dimension(d_n_max), intent(out) :: dv_in,dv_ou,dm_in,th_en    ! layer inflow/outflow (m3/s)
	real(r8) :: rho_r,in_v,m_in    !
	integer :: j,jmax,jmin							! indices
	! integer, intent(out)  :: jminn,jmaxx

	rho_r = 1000.*( 1.0 - 1.9549e-05*(abs(in_t-277.))**1.68 )
	in_v = in_f*dtime
	m_in = in_v*rho_r
!	Initialize
    do j=1,d_n_max
		dv_in(j)=0.
		dv_ou(j)=0.
		dm_in(j)=0.
		th_en(j)=0.
	end do

!   Layer inflow and energy contribution
    jmin=1
	if(d_n>3)jmax=d_n-3
	if(d_n<=3)jmax=d_n-1

	if(d_n==1) then ! Single layer
		do j=1,d_n_max
			dv_in(j)=0.
			dm_in(j)=0.
            th_en(j)=0.
			dv_ou(j)=0.
		end do
		dv_in(d_n)=in_f
		dm_in(d_n)=m_in
        th_en(d_n)=dv_in(d_n)*rho_r*in_t*c_w
		dv_ou(d_n)=ou_f
	else
		do j=jmin,jmax
            dv_in(j)=in_f*(d_v(j)/(v_zt(jmax+1)-v_zt(jmin)))
			dm_in(j)=m_in*(d_v(j)/(v_zt(jmax+1)-v_zt(jmin)))
            th_en(j)=dv_in(j)*rho_r*in_t*c_w
			dv_ou(j)=ou_f*(d_v(j)/(v_zt(jmax+1)-v_zt(jmin)))
        end do
	end if

end subroutine flowdist

function den(t_z) result(rho)
	! calculate density from temperature
	implicit none
    integer,parameter :: r8 = selected_real_kind( 16) ! 8 byte real
	real(r8), intent(in) :: t_z! Temperature (k)
    real(r8) :: rho! ! density (kg/m3)

	 rho = 1000.*( 1.0 - 1.9549e-05*(abs(t_z-277.))**1.68 ) ! modified from Subin et al, 2011 with lake ice fraction = 0
end function den
