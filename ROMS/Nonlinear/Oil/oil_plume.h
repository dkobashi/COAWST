      SUBROUTINE oil_param_initialize(ng, Lstr, Lend)
!***********************************************************************                                              
!
! Initialize oil parameters that are not changed during the simulation
! 
!***********************************************************************                             
      USE mod_param
      USE mod_floats
      USE mod_scalars
      USE mod_parallel
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
!

#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10, __LINE__, __FILE__)
#endif
      CALL oil_paramini_tile (ng, Lstr, Lend,                          &
     &                      DRIFTER(ng) % frroil0,                     &
     &                      DRIFTER(ng) % wfroil0,                     &
     &                      DRIFTER(ng) % szoil0)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10, __LINE__, __FILE__)
#endif
 

      END SUBROUTINE oil_param_initialize

!***********************************************************************
      SUBROUTINE oil_paramini_tile(ng, Lstr, Lend,                   &
     &                            frroil0, wfroil0, szoil0)
!***********************************************************************
      USE mod_param
      USE mod_parallel
      USE mod_floats
!      USE mod_grid
!      USE mod_iounits
      USE mod_scalars
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY: mp_bcastf
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: frroil0(:,:)  
      real(r8), intent(inout) :: wfroil0(:,:)  
      real(r8), intent(inout) :: szoil0(:)     
#else
      real(r8), intent(inout) :: frroil0(Nfloats(ng),Nocmp)  
      real(r8), intent(inout) :: wfroil0(Nfloats(ng),Nocmp)  
      real(r8), intent(inout) :: szoil0(Nfloats(ng))         
#endif
!
! Local variables
!  rhowt - water density
!  temp, salt - water T/S
!  rhoo - oil particle density
!
!      real(r8), parameter :: pi=3.14159265358979, &
!      real(r8), Parameter :: gg=9.81

      integer :: l

      real(r8) :: rhoo
      real(r8) :: Doil, Docrit, Doil0
      real(r8) :: frsats, frarom, frasph
      real(r8) :: zoil, tsrfo, rsats, rarom, rasph

      IF (MyRank.eq.MyMaster) print *,                                 &
     &    'oil_paramini_tile Lstr=',Lstr,'Lend=',Lend,'MyRank=',MyRank

      DO l=Lstr,Lend
        IF (MyRank.eq.MyMaster) THEN
! Initial Oil density is calculated using mixing rules
! based on mass fraction of 3 types of oil components
! (SARA with resins & asphaltenes components combined)
! saturates, aromatics, and resins+asphalts
! Approach similar to (simplified) to 
! that discussed in Yarranton et al., 2015 "Density and refraction index of petroleum ..."
! I use Reddy et al. mass fractions for Macondo oil (from samples during Deepwater Horiz spill)
! Density of oil components are listed in Yarranton, and here 
! densities are scaled so that 
! the oil density matches the mean Macondo oil rho=858 kg m3 (from North et al., 2015)
! or ~820 reported in Reddy et al., PNAS, 2012
! For each individual float the mass fraction is randomly assigned
! by changing fraction of saturates (resins+asphaltenes are unchanged)
!
! Size - according to Vilcaez et al., 2013 "A new model for the biodegradation kinetics ..."
! in Macondo oil (after the dispersants have been added), oil droplet sizes followed
! Gamma distribution with mean droplet size = 20 micrometers <-- give very slow woil
! Adjusted to mean Doil=300 mcrmeters 
! 
            frsats=-1.0 ! to initialize oil fractions <0
            frarom=-1.0
            frasph=-1.0
            CALL oil_density_ini(rhoo,frsats,frarom,frasph,           &
     &                       rsats,rarom,rasph)
            CALL oil_size_ini(Doil)
            frroil0(l,1)=rsats
            frroil0(l,2)=rarom
            frroil0(l,3)=rasph
            wfroil0(l,1)=frsats
            wfroil0(l,2)=frarom
            wfroil0(l,3)=frasph
            szoil0(l)=Doil
        ENDIF
      ENDDO ! Floats loop
!
      IF (MyRank.eq.MyMaster) print *,'oil initialization is complete '
# ifdef DISTRIBUTE
!
! Broadcast this to all processors
      CALL mp_bcastf(ng, iNLM, frroil0)
      CALL mp_bcastf(ng, iNLM, wfroil0)
      CALL mp_bcastf(ng, iNLM, szoil0)
!      IF (MyRank.eq.MyMaster) print *,'oil_paramini_tile: broadcasting 3'
# endif
!
      END SUBROUTINE oil_paramini_tile


      SUBROUTINE oil_floats (ng, Lstr, Lend, Predictor, my_thread,ifltX)
!
!**************************************************  ***
!***********************************************************************
!                                                                      !
!  This routine calculates vertical velocity of oil partilces          !
!  I use algorithm from Zheng and Yapa, 2000, "Buoynat velocity of     !
!  spherical and nonspherical bubbles/droplets", J. Hydraulic Enginee- !
!  ring, 126(11), 852-854                                              !
! 2 methods are suggested: 2-equation and integrated                   !
! integrated method predicts slower vertical velocities for the oil    !
! droplets, however some extreme cases are not covered (high Reyn.#)   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_floats
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend, ifltX  ! DDTMIRY ifltX

      logical, intent(in) :: Predictor
#ifdef ASSUMED_SHAPE
      logical, intent(in) :: my_thread(Lstr:)
#else
      logical, intent(in) :: my_thread(Lstr:Lend)
#endif
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10, __LINE__, __FILE__)
#endif
      CALL oil_floats_tile (ng, Lstr, Lend,                            &
     &                      nfm3(ng), nfm2(ng), nfm1(ng), nf(ng),      &
     &                      nfp1(ng),                                  &
     &                      Predictor, my_thread,                      &
     &                      DRIFTER(ng) % bounded,                     &
     &                      DRIFTER(ng) % Tinfo,                       &
     &                      DRIFTER(ng) % track,                       &
     &                      DRIFTER(ng) % frroil0,                     &
     &                      DRIFTER(ng) % wfroil0,                     &
     &                      DRIFTER(ng) % szoil0,                      &
     &                      ifltX)                                     
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10, __LINE__, __FILE__)
#endif

      RETURN
      END SUBROUTINE oil_floats
!
!***********************************************************************
      SUBROUTINE oil_floats_tile (ng, Lstr, Lend,                   &
     &                            nfm3, nfm2, nfm1, nf, nfp1,       &
     &                            Predictor, my_thread, bounded,    &
     &                            Tinfo, track,                     &
     &                            frroil0, wfroil0, szoil0,         &  ! DDMITRY
     &                            ifltX)     ! DDMITRY
!***********************************************************************
      USE mod_param
      USE mod_parallel
!      USE mod_oil_floats_param
      USE mod_floats
      USE mod_grid
      USE mod_iounits
      USE mod_scalars
! DDMITRY
      USE nrutil, ONLY: ran1
! END
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
      integer, intent(in) :: nfm3, nfm2, nfm1, nf, nfp1
      integer, intent(in) :: ifltX
      logical, intent(in) :: Predictor
!
#ifdef ASSUMED_SHAPE
      logical, intent(in) :: bounded(:)
      logical, intent(in) :: my_thread(Lstr:)

      real(r8), intent(in) :: Tinfo(0:,:)
      real(r8), intent(inout) :: frroil0(:,:)  
      real(r8), intent(inout) :: wfroil0(:,:)  
      real(r8), intent(inout) :: szoil0(:)     
      real(r8), intent(inout) :: track(:,0:,:)
#else
      logical, intent(in) :: bounded(Nfloats(ng))
      logical, intent(in) :: my_thread(Lstr:Lend)

      real(r8), intent(in) :: Tinfo(0:izrhs,Nfloats(ng))
      real(r8), intent(inout) :: frroil0(Nfloats(ng),Nocmp)  
      real(r8), intent(inout) :: wfroil0(Nfloats(ng),Nocmp)  
      real(r8), intent(inout) :: szoil0(Nfloats(ng))         
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
#endif
!
!  Local variable declarations.
!  rhowt - water density
!  temp, salt - water T/S
!  rhoo - oil particle density
!
!      real(r8), parameter :: pi=3.14159265358979, &
      real(r8), Parameter :: gg=9.81

      integer :: i, i1, i2, j1, j2, l, ii, jj
      integer :: iFltS  ! temporary var for tracking surf float

      real(r8) :: cff1, cff2, cff3, cff4
      real(r8) :: rhowt, temp, salt, rhoo
      real(r8) :: Doil, Docrit, Doil0
      real(r8) :: dwoil, woil, woil_b, dTime
      real(r8) :: drho, gprime, muw, musw, signw
      real(r8) :: skg, A1, B1, sigmOW
      real(r8) :: HalfDT, Fbuoy, Fdrag
      real(r8) :: v1, frsats, frarom, frasph
      real(r8) :: zoil, tsrfo, rsats, rarom, rasph
      real(r8) :: flat, flon
      real(r8) :: Rhoc0(Nocmp), Froil0(Nocmp), Froil(Nocmp)
      real(r8) :: tsrfx   ! max surf time

      logical :: DiagSrf
!
!-----------------------------------------------------------------------
!
!  Calculate vertical velocity from positive buoyancy of the particle
!
!-----------------------------------------------------------------------
!
!  Assign predictor/corrector weights.
!
      IF (Predictor) THEN
        cff1=8.0_r8/3.0_r8
        cff2=4.0_r8/3.0_r8
      ELSE
        cff1=9.0_r8/8.0_r8
        cff2=1.0_r8/8.0_r8
        cff3=3.0_r8/8.0_r8
        cff4=6.0_r8/8.0_r8
      END IF
!
!
      HalfDT=0.5_r8*dt(ng)
! ------
! Checking evap: track a float with max surface exposure
! Will delete later - debugging only
      tsrfx=0.0_r8
      iFltS=0
      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
          IF (Predictor) THEN          
            tsrfo=track(isrfo,nf,l)
          ELSE
            tsrfo=track(isrfo,nfp1,l)
          ENDIF
        ENDIF
        IF (tsrfo.gt.tsrfx+HalfDT) THEN
          iFltS=l
          tsrfx=tsrfo
        ENDIF
      ENDDO

!      iFltS=4116
!
!  ---- END --- 
!
      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
!
!  If newly relased float, initialize oil fields. 
! Oil initial parameters that do not change during the simulation
! should be already initialized and distributed over
! all processors
!  Note that since we need temperature and salinity, we need to initialize
!  their values to all time levels. Otherwise, we will have a parallel
!  bug.
          IF (time(ng)-HalfDT.le.Tinfo(itstr,l).and.                    &
     &        time(ng)+HalfDT.gt.Tinfo(itstr,l)) THEN
            temp=track(ifTvar(itemp),nfp1,l)
            salt=track(ifTvar(isalt),nfp1,l)
! 
            frsats=-1.0 ! to initialize oil fractions <0
            frarom=-1.0
            frasph=-1.0
            CALL oil_density_ini(rhoo,frsats,frarom,frasph,           &
     &                       rsats,rarom,rasph)
            CALL oil_size_ini(Doil)
!
            rsats=frroil0(l,1)
            rarom=frroil0(l,2)
            rasph=frroil0(l,3)
            frsats=wfroil0(l,1)
            frarom=wfroil0(l,2)
            frasph=wfroil0(l,3)
            Doil=szoil0(l)
            flat=track(iflat,nfp1,l)
            flon=track(iflon,nfp1,l)
            zoil=track(idpth,nfp1,l)

            DO i=0,NFT  ! time levels
              track(ifTvar(itemp),i,l)=temp
              track(ifTvar(isalt),i,l)=salt
              track(isats,i,l)=frsats
              track(iarom,i,l)=frarom
              track(isizo,i,l)=Doil
              track(iroil,i,l)=rhoo
              track(isrfo,i,l)=0.0
              track(iflon,i,l)=flon
              track(iflat,i,l)=flat
              track(idpth,i,l)=zoil
            END DO
          END IF

!
! Get T, S, Rho water
!
          temp=track(ifTvar(itemp),nfp1,l)
          salt=track(ifTvar(isalt),nfp1,l)
          CALL water_dens0(temp,salt,rhowt)

! Calculate dynamic viscosity for sea water
! I use: Sharqawy, Mostafa H., Lienhard V, John H.,
! and Zubair, Syed M., 2010, “Thermophysical
! properties of seawater: a review of existing 
! correlations and data,” Desalination
! and Water Treatment, Vol. 16, pp. 354-380.
! To calculate dynamic (abs) sea water visc. 
          skg = salt*1e-3;                          ! salinty -> kg/kg
          A1 = 1.541+1.998e-2*temp-9.52e-5*temp**2;
          B1 = 7.974-7.561e-2*temp+4.724e-4*temp**2;
          muw = 4.2844e-5+(0.157*(temp+64.993)**2-91.296)**(-1); ! dynamic visc. pure water
          musw = muw*(1.0_r8+A1*skg+B1*skg**2);                  ! dyn. visc. sea water
!
! water-oil interfacial tension 
! Eq. from Peters and Arabali, "Interfacial tension between oil
! and water...", Colloids and Surfaces, A, 426, 2013
          sigmOW=(0.1222*temp+32.82)*1e-3
!
! Oil W variables
          IF (Predictor) THEN
            rhoo=track(iroil,nf,l)
            Doil=track(isizo,nf,l)
            frsats=track(isats,nf,l)
            frarom=track(iarom,nf,l)
            zoil=track(idpth,nf,l)
            tsrfo=track(isrfo,nf,l)
          ELSE
            rhoo=track(iroil,nfp1,l)
            Doil=track(isizo,nfp1,l)
            frsats=track(isats,nfp1,l)
            frarom=track(iarom,nfp1,l)
            zoil=track(idpth,nfp1,l)
            tsrfo=track(isrfo,nfp1,l)
          ENDIF
!
! Checking: ++++++++ 
          DO ii=1,Nocmp
            Rhoc0(ii)=frroil0(l,ii)
            Froil0(ii)=wfroil0(l,ii)
            IF (Rhoc0(ii)<1.0 .or. Rhoc0(ii)>1200.0) THEN
              print *,'*** ERR: Rhoc0 invalid',Rhoc0(ii),'ii=',ii
              print *,'*** ERR: Froil=',Froil(ii)
              print *,'*** ERR: float # l=',l
            ENDIF
          ENDDO
!          IF (l.eq.4116) THEN
!            print *,'Time=',time(ng),'l=',l,' frroil0=',frroil0(l,1:3)
!          ENDIF
! +++++++


! Initial oil characteristics
          Doil0=szoil0(l)
!
! Current oil characteristics          
          Froil(1)=frsats
          Froil(2)=frarom
          Froil(3)=1.0_r8-(frsats+frarom)
!
!!!!!!!!!!!!!
! Update rho oil due to weathering processes
! at the surface
          DiagSrf=.FALSE.
          IF (abs(zoil).le.5.0) THEN
            tsrfo=tsrfo+dt(ng)   ! integrated time at the surface, sec
!
            IF (iFltS==0) iFltS=l ! Check evaporation, track surf float
# ifdef OIL_DEBUG
            IF (iFltS==l) THEN 
              DiagSrf=.TRUE.
              print *,'Checking evaporation, Surf float: ',iFltS
            ELSE
              DiagSrf=.FALSE.
            ENDIF
# endif
!
!  Evaporation: update rhoo, Doil, Froil
            CALL oil_evaporation(Rhoc0,Froil0,Doil0,Froil,Doil,rhoo,    &
     &                           tsrfo,temp,DiagSrf,l)
          ENDIF

          drho=rhowt-rhoo
          gprime=gg*drho/rhowt  ! reduced gravity
          signw=drho/abs(drho)
!
! Oil vertical velocity
!
#ifdef WOIL_INTEGRATED
! Integrated approach 
          CALL woil_integrated(rhowt, rhoo, muw, musw, Doil, gg,  &
     &                           woil, signw, sigmOW)
#else
! Two-equation approach following 
! Zheng and Yapa, 2000, J. Hydraulic Engineering, 852-854
! for small Re (small particles)
! Use Stokes law
! for Large Re (large particles)
! use Cd = const = 0.5, different eq. 
! Critical diameter is Docrit
          Docrit=(9.52_r8*musw**(2.0/3.0))/  &
     &            ((gg*rhowt*abs(drho))**(1.0/3.0))

          IF (Doil.lt.Docrit) THEN
            woil=gg*(Doil**2)*drho/(18.0_r8*musw)
          ELSE
            woil=signw*(8.0_r8/3.0_r8*abs(gprime)*Doil)**(1.0/2.0)
          ENDIF
#endif
!
! Update track array:
! surface time, oil fractions, etc.
          track(isrfo,nfp1,l)=tsrfo
          track(iwoil,nfp1,l)=woil
          track(isats,nfp1,l)=Froil(1)
          track(iarom,nfp1,l)=Froil(2)
          track(isizo,nfp1,l)=Doil
          track(iroil,nfp1,l)=rhoo

# ifdef OIL_DEBUG
          IF (l.eq.ifltX) THEN
            print*,'    '
            IF (Predictor) print*,'####  oil_plume: Predictor step'
            IF (.not.Predictor) print*,'#### oil_plume: Corrector step'
            print *,'### oil_plume: l=',l,                             &
     &         'T=',temp, 'S=',salt,                                   &
     &         'RhoWt=',rhowt,'% sats=',frsats,  '% arom=',frarom,     &
     &         'Doil=',Doil,'RhoOil=',rhoo,'Time Surf=',tsrfo/60.0,    &
     &         '  woil=',woil 
            print*,'    '
          ENDIF
# endif

!
!  If newly relased oil float, set vertical w fields for all time
!  levels.
!
          IF (time(ng)-HalfDT.le.Tinfo(itstr,l).and.                    &
     &        time(ng)+HalfDT.gt.Tinfo(itstr,l)) THEN
            woil=track(iwoil,nfp1,l)
            DO i=0,NFT
              track(iwoil,i,l)=woil
            END DO
          END IF


        ENDIF ! if my_thread(l) & bounded(l)
!
! Checking DRIFTER pointer, should be distributed over
! all nodes, same initial values
!        IF (l==4116 .and. bounded(l)) THEN
!          print*,'Check DRIFT: MyRank=',MyRank,'my_thread=',my_thread(l)
!          print*,'::: x=',Tinfo(ixgrd,l),' frroil(1)=',frroil0(l,1), &
!     &           'Doil0=',szoil0(l),  &
!     &           ' Wfroil0=',wfroil0(l,1)
!        ENDIF

      ENDDO  ! l=1,N floats

      RETURN
      END SUBROUTINE oil_floats_tile
!
!
!***********************************************************************
      SUBROUTINE oil_evaporation(Rhoc0,Froil0,Doil0,Froil,Doil,rhoo,    &
     &                           tsrfo,temp,flgdg,l)
!***********************************************************************
! SARA-component evaporation
! Evaporation is calculated for each component
! oil size and density is adjusted
! accordingly
! Adopted approach from
! from Marv Fingas book 2011
! Ch. 9 Evaporation Modeling
! SARA mass fractions are not percents (as in Fingas formulae)
! Note that in my approach model resins and asphaltenes are combined
! 2 steps: calculate parameter Evap for each SARA (in Fingas, eq. 39)
!          using densities and viscosities (use eq. from Sanchez-Minero et al., 2014)
!          Then, use eq. 38 (Fingas) to calculate evaporated % of oil components
!
! Steps: 1) Calculate parts of evaporated oil component (pevp)
! 2) Update change in weight by components (mfr) -> change in overall weight (moil)
! 3) Calculate new weight fractions for SARA components after evaporation
! 4) Update oil density and droplet size
!
! INPUT:
! Rhoc0 - initial oil dens by components
! Froil0 - initial wegith fraction (0</= parts </= 1) by components
! Froil  -  current Weight fractions of saturates, aromatics, asphaltenes
! Doil   - current oil size (m)
! rhoo   - current total oil droplet density
! tsrfo - time (sec) since float is at the surface (Note need sec -> min)
!  temp - ocean T (however, surface oil T may be higher)
! 
! OUTPUT - updated state after evaporation: 
! Froil, Doil, rhoo
!
      USE mod_param
      USE mod_scalars
      USE mod_floats

      integer, intent(in) :: l

      real(r8), intent(in) :: temp, tsrfo
      real(r8), intent(in) :: Rhoc0(Nocmp)  
      real(r8), intent(in) :: Froil0(Nocmp)  
      real(r8), intent(inout) :: Doil, rhoo
      real(r8), intent(inout) :: Froil(Nocmp)

      logical, intent(inout) :: flgdg
! Local variables:      
      real(r8) :: Evap, rhowt0
      real(r8) :: TK, aa, bb, api, SG, rhoxx
      real(r8) :: muOil, dmm
      real(r8) :: tsrfo1, pevp
      real(r8) :: rhoo0, Doil0, voil0, moil0
      real(r8) :: SPIg(Nocmp), PEVC(Nocmp)
      real(r8) :: Mfr0(Nocmp), Mfr1(Nocmp), Froil1(Nocmp)
      real(r8) :: moil1, rhoo1, voil1, Doil1

      integer :: ii

! Given information by components (weight fraction and density)
! Get oil droplet weight (moil0), density, etc
      CALL oil_comp2total(Froil0,Rhoc0,Mfr0,rhoo0,Doil0,voil0,moil0)

! calculate SPI gravity
      rhowt0=1000.0_r8
      SPIg=Rhoc0/rhowt0

      SG=-1.e23
      DO ii=1,Nocmp ! evaporation of individual oil components
        SG=SPIg(ii)
        rhoxx=Rhoc0(ii)        
        api=141.5_r8/SG-131.5_r8  !API gravity
!
! Calc viscosity based on Sanchez-Minero et al, 2014
        TK=temp+273.15_r8
        aa=3.9e-5_r8*api**3-4.0e-3_r8*api**2+0.1226_r8*api-0.7626_r8
        IF (aa<0.001) aa=0.001 ! super-viscous and dense oil (>1010) may give aa<0
        bb=9.1638e9_r8*api**(-1.3257_r8)
        muOil=aa*exp(bb/TK**3) ! oil component viscosity, cPoise
!
! Fingas, 2011, ch. 9
! Evaporation based on dens and viscosity
! Parameter Evap for each component
        Evap=15.4_r8-14.5_r8*rhoxx/1000.0_r8+2.58_r8/muOil

! Finga's formula calculates fraction evaporated (fr_evap) over the total time from t=0 -> t1
        tsrfo1=tsrfo/60.0_r8 ! sec->min
        IF (tsrfo1<1.e-10) tsrfo1=1.e-10

        pevp=((Evap+0.675_r8)+0.045_r8*(temp-15.0_r8))*   &
     &         log(tsrfo1)/100.0_r8                            !  frac. evaporated over total time

        PEVC(ii)=pevp

# ifdef OIL_DEBUG
        IF (isnan(pevp)) THEN
          flgdg=.true.
          print*,'*** NaN pevp *** Float l=',l,'Evap=',Evap
          print*,'*** rhoxx=',rhoxx,'muOil=',muOil
          print*,'aa=',aa,'bb=',bb,'api=',api,' SG=',SG
          print*,'*** ERR  pevp*** temp=',temp,' tsrfo1=',tsrfo1
          print*,'   ********** '
       ENDIF
# endif        

      ENDDO  ! evaporation by components

! Update oil mass fraction reduction due to evaporation
      Mfr1=(1.0_r8-PEVC)*Mfr0  
      moil1=sum(Mfr1)            ! total weight oil droplet
      Froil1=1.0_r8/moil1*Mfr1   ! new weight fractions
      dmm=sum(Froil1/Rhoc0)
      rhoo1=1.0_r8/dmm           ! new oil density
      voil1=moil1/rhoo1
      Doil1=2.0_r8*(3.0_r8/4.0_r8*voil1/pi)**(1./3.)

      rhoo=rhoo1
      Doil=Doil1
      Froil=Froil1

# ifdef OIL_DEBUG
      IF (flgdg) THEN
        print*,' ---    Evaporation Subr, surf time, min:',tsrfo1
        print*,' Evaporated:',PEVC
        print*,'Weight Oil comp, initial:',Mfr0
        print*,'Weight Oil comp, new:',Mfr1
        print*,'Weight Fraction, initial:',Froil0
        print*,'Weight Fraction, new:',Froil1
        print*,'Oil Rho, initial:',rhoo0
        print*,'Oil Rho, new:',rhoo1
        print*,'Oil Vol,m3, initial:',voil0
        print*,'Oil Vol,m3, new:',voil1
        print*,'Oil Size,m, initial:',Doil0
        print*,'Oil Size,m, new:',Doil1
        print*,' ------------'
        print*,' '
      ENDIF
# endif

      RETURN
      END SUBROUTINE oil_evaporation
!
!
!  
!***********************************************************************
      SUBROUTINE oil_comp2total(Froil,Rhoc,Mfr,rhoo,Doil,voil,moil) 
!***********************************************************************
! Input: given oil size Doil
! and oil component fractions (Froil) and oil component densities (Rhoc)
! Output:
! into total weight (moil), volume (voil), density (rhoo)
! also weight by components (Mfr)
!
      USE mod_param
      USE mod_scalars
      USE mod_floats

      real(r8), intent(in) :: Froil(Nocmp), Rhoc(Nocmp)
      real(r8), intent(inout) :: rhoo, Doil, voil, moil
      real(r8), intent(inout) :: Mfr(Nocmp) 

      real(r8) :: dmm

      dmm=sum(Froil/Rhoc) 
      rhoo=1.0_r8/dmm

      voil=4.0_r8/3.0_r8*pi*(Doil/2.0_r8)**3
      moil=rhoo*voil
      Mfr=moil*Froil

      END SUBROUTINE oil_comp2total

#ifdef WOIL_INTEGRATED
!
!***********************************************************************
      SUBROUTINE woil_integrated(rhowt, rhoo, muw, musw, Doil, gg, &
     &                           woil, signw, sigmOW)
!***********************************************************************
      USE mod_param
!
! Imported variables declaration
!
      real(r8), intent(in) :: rhowt, rhoo, muw, musw
      real(r8), intent(in) :: Doil, signw, sigmOW, gg
      real(r8), intent(inout) :: woil

      real(r8) :: ND, logW, Re, E0, MW, HH, JJ
      real(r8) :: logRe, vbound, drho

      vbound = 1.11111e-6 ! out of bound value
      drho=rhowt-rhoo

      IF (Doil<1.1e-3) THEN ! Small size spherical range
        ND=4.0_r8*rhowt*drho*gg*(Doil**3)/(3.0_r8*musw**2)
        logW=log10(ND)

        IF (ND.le.73.0) THEN
          Re=ND/24.0_r8-1.7569e-4*(ND**2)+6.9252e-7*(ND**3)-            & 
     &        2.3027e-10*ND**4
        ELSEIF (ND.gt.73.0 .and. ND.le.580.0) THEN
          logRe=-1.7095+1.33438*logW-0.11591*(logW**2)
          Re=10**(logRe)
        ELSEIF (ND.gt.580.0 .and. ND.le.1.55e7) THEN
          logRe=-1.81391+1.34671*logW-0.12427*(logW**2)+                & 
     &            0.006344*(logW**3)
          Re=10**(logRe)
        ELSE  ! out of range
          Re=-1.0
        ENDIF

        IF (Re>0.0) THEN
          woil=signw*Re*musw/(rhowt*Doil)
        ELSE
          woil=vbound   ! out of range value
        ENDIF
      
      ELSE
        E0=gg*abs(drho)*(Doil**2)/sigmOW
        MW=gg*(musw**4)*abs(drho)/((rhowt**2)*(sigmOW**3))

        IF (MW.lt.1.0e-3 .and. E0.le.40.0) THEN
          HH=4.0/3.0*E0*MW**(-0.149)*(muw/musw)**(-0.14)

          IF (HH.gt.2.0 .and. HH.le.59.3) THEN
            JJ=0.94*(HH**0.757)
          ELSE
            JJ=3.42*(HH**0.441)
          ENDIF

          woil=signw*musw/(rhowt*Doil)*(MW**(-0.149))*(JJ-0.857)

        ELSEIF (E0.gt.40.0) THEN
          woil=signw*0.711*sqrt(gg*Doil*drho/rhowt)
        ELSE  ! out of bound case
          woil=vbound
        ENDIF

      ENDIF

      RETURN
      END SUBROUTINE woil_integrated
!
#endif

!                                                                                                                                 
!***********************************************************************
      SUBROUTINE oil_density_ini(rhoo,frsats,frarom,frasph, &
     &                           rsats,rarom,rasph)
!***********************************************************************
! Initialize oil density for Macondo oil
! Using a 3- component (~SARA with resins and asphaltenes combined)
! oil structure
! In the future, this can be modified by an N-component oil
! specified in an input file 
!
      USE mod_param
      USE nrutil, ONLY : ran1
      USE mod_floats

!
! Imported variables declaration
!
      real(r8), intent(inout) :: rhoo, frsats, frarom,frasph
      real(r8), intent(inout) :: rsats, rarom, rasph
!
! Local variables 
!
      real(r8) :: v1
!
! Macondo oil density estimates varies in the literature
! North et al,. 2015: mean roil = 858 kg m3
! Reddy et al., 2013: sampled oil roil = 820 kg m3
! I'm trying to match Reddy's mean Rho oil = 820 kg m3
! North's estimates result in low woil and late
! surfacing of oil 
!
! This will need to be changed
! For now, 3-component SAR+A oil composition
! is assumed, Numb of components is
! specified in oil*.in however current code
! works for 3 components
      rsats=RhoOilComp(1)
      rarom=RhoOilComp(2)
      rasph=RhoOilComp(3)

! Uniform Randomly assign fraction of 
! aromatics and resins+asphalts fractions
! sats - the rest
! for initializing 
      IF (frsats.lt.0.0 .or. frarom.lt.0.0 .or. frasph.lt.0.0) then
        CALL ran1 (v1)
        frarom=0.16_r8+(v1-0.5_r8)*0.15_r8
        CALL ran1 (v1)
        frasph=0.1_r8+(v1-0.5_r8)*0.08_r8
        frsats=1.0_r8-(frarom+frasph)
!
!    Make sure that the sum of fractions = 1 - will need to add check and adjustment
! 
      ENDIF

      v1=frsats/rsats+frasph/rasph+frarom/rarom
      rhoo=1.0_r8/v1      

      END SUBROUTINE oil_density_ini

!*********************************************************************** 
      SUBROUTINE oil_size_ini(osize)
!*********************************************************************** 
! Gamma distribution with mean size=20 micrometer
! Based on observations/lab data, Vilcaez et al., 2013
! small droplet sizes result in very slow woil
! Other papers suggest much larger droplet sizes
! e.g., Zhao et al., 2016 "Underwater oil jet ..."
! Doil ~200-700 microm with different PDF
! 
! To match timing of oil surfacing 
! I increase mean Doil to 300 and 500
! 
! I use algorithm for generating Gamma Random variable
! described in 
! Marsaglia and Tsang, 2000, A simple method for
! generating Gamma variables
! Transactions on Mathematical Software, 26(3), 363-372
! adopted for alfa>1 and alfa<1
! Gamma distr = 1/(betta^alf*Gamma(alf))*x^(alf-1)*exp(-x/betta)
! I use Gamma formulae with lambda=1/betta

      USE mod_param
      USE mod_scalars
      USE mod_iounits
      USE mod_floats
      USE nrutil, ONLY : ran1, gasdev
!
! Imported variables declaration 
!
      real(r8), intent(inout) :: osize

! 
! Local variables
!
      real(r8), parameter :: alf0=4.94   ,  &     ! Gamma parameter - controls St Dev spread, decreases with alf0++
     &                       Dmin=2.0e-7          ! minimum droplet size, oil dissolves <0.19 mcrm
      real(r8) :: v1, sgmm, lmbd1, alf, lmm, Doil0
      real(r8) :: dd1, cc1, xnorm, xunif, vv1, rgm1
      integer :: cc

      logical :: lalf1

! Mean oil droplet size is user-specified in
! oil*.in
      Doil0=DoilMn     

      lalf1=.false.
      alf=alf0
      IF (alf0.le.1.0_r8) THEN
        lalf1=.true.
        alf=alf0+1.0_r8
      ENDIF

      sgmm=Doil0/(alf**0.5)
      lmbd1=alf/Doil0        ! Gamma parameter

      dd1=alf-1.0_r8/3.0_r8
      cc1=1.0_r8/((9.0_r8*dd1)**0.5)

      rgm1=-1.0
      cc=0
      DO WHILE (rgm1.lt.0.)
        xnorm=-100./cc1
        cc=cc+1
        IF (cc>10) THEN
          WRITE(stdout,*) 'Gamma distr oil size, endless loop cc=',cc
          exit_flag=8
          RETURN
        ENDIF          

        DO WHILE (xnorm.le.-1.0/cc1)
          CALL gasdev(xnorm)
        ENDDO
        vv1=(1.0_r8+cc1*xnorm)**3


        CALL ran1 (xunif)
        lmm=0.5*xnorm**2+dd1*(1.0_r8-vv1+log(vv1))
        IF (log(xunif).lt.lmm) THEN
          rgm1=dd1*vv1/lmbd1

          IF (lalf1) THEN
            CALL ran1 (v1)
            rgm1=rgm1*v1**(1.0_r8/alf0)
          ENDIF

          osize=rgm1
          IF (osize.lt.Dmin) osize=Dmin
          EXIT 

        ENDIF

      ENDDO



      END SUBROUTINE oil_size_ini
!
!***********************************************************************
      SUBROUTINE water_dens0(temp,salt,rhowt)
!***********************************************************************
      USE mod_param
      USE mod_scalars
!
! Imported variables declaration
!      
      real(r8), intent(in) :: temp, salt
      real(r8), intent(inout) :: rhowt
!
! Local variables declaration
!
      real(r8), parameter :: b0=8.24493e-1, b1=-4.0899e-3
      real(r8), parameter :: b2=7.6438e-5, b3=-8.2467e-7
      real(r8), parameter :: b4=5.3875e-9, c0=-5.72466e-3
      real(r8), parameter :: c1=1.0227e-4, c2=-1.6546e-6
      real(r8), parameter :: d0=4.8314e-4
      real(r8), parameter :: a0=999.842594, a1= 6.793952e-2
      real(r8), parameter :: a2=-9.095290e-3, a3=1.001685e-4
      real(r8), parameter :: a4=-1.120083e-6, a5=6.536332e-9

      real(r8) :: t68, rhow0s ! standard water dens, S=0

      t68=temp*1.00024
      rhow0s=a0+(a1+(a2+(a3+(a4+a5*t68)*t68)*t68)*t68)*t68
      rhowt=rhow0s+(b0+(b1+(b2+(b3+b4*t68)*t68)*t68)*t68)*salt+&
                  (c0+(c1+c2*t68)*t68)*salt*sqrt(salt)+d0*salt**2 

      RETURN
      END SUBROUTINE water_dens0



