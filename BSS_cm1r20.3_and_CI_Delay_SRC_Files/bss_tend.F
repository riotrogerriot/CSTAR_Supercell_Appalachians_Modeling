 MODULE bss_tend_module

! This subroutine is designed to read in the original base-state environment (via
! input_sounding) and any subsequent environments (thermo and wind separately, to
! assist with any controlled tests). The total difference in u, v, th, etc. is 
! calculated between two subsequent profiles (initial vs. 2nd, 2nd vs. 3rd, etc.)
! and then divided into increments to determine the rate of change (i.e., the 
! tendency of each base-state variable. The actual separation of perturbation
! from the base-state and the actual base-state substitution takes place in the
! time loop (see cm1.F).

! This is a sizeable update from the original BSS formulation (see Letkewicz et al.
! 2013 in Mon. Wea. Rev.)--April 2017 (see Davenport et al. 2019 in MWR)

   implicit none

   CONTAINS

      subroutine bss_tend(bss_start,bss_end,zh,sndg2time,sndg3time, &
                      windchange,thermochange,u_tend1,u_tend2,      &
                      v_tend1,v_tend2,pi_tend1,pi_tend2, prs_tend1, &
                      prs_tend2,rho_tend1,rho_tend2,th_tend1,       &
                      th_tend2,qv_tend1,qv_tend2,udiff,vdiff,       &
                      pidiff,prsdiff,rhodiff,thdiff,qvdiff,         &
                      u_start,v_start,pi_start,prs_start,rho_start, &
                      thv_start,th_start,qv_start,                  &
                      u_end,v_end,pi_end,prs_end,rho_end,           &
                      thv_end,th_end,qv_end)

      use input
      use constants
      use bc_module
      use comm_module
      use goddard_module, only : T0K,T00K,RT0
      use cm1libs , only : rslf,rsif
      use turb_module , only : ramp_up_turb
      use getcape_module


#ifdef MPI
      use mpi
#endif
      implicit none

      real bss_start, bss_end  !BSS start/end times (in sec) 
      real sndg2time, sndg3time  !Model run time (since beginning) to implement 2nd, 3rd sndg 
                                 !NOTE: current formulation below says 3rd sndg is essentially 
                                 !when BSS ends (i.e., changes from orig to 2nd, 2nd to 3rd only)
      integer windchange, thermochange !check to change wind/thermo profile for BSS

      real :: total_time

      double precision, dimension(ib:ie+1,jb:je,kb:ke) :: u_tend1, u_tend2
      double precision, dimension(ib:ie,jb:je+1,kb:ke) :: v_tend1, v_tend2
!      real, dimension(ib:ie+1,jb:je,kb:ke) :: u_tend1, u_tend2
!      real, dimension(ib:ie,jb:je+1,kb:ke) :: v_tend1, v_tend2


      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, dimension(ib:ie+1,jb:je,kb:ke) :: udiff
      real, dimension(ib:ie,jb:je+1,kb:ke) :: vdiff

      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,prs0,rho0,thv0,th0,rth0,qv0
      real, dimension(ib:ie,jb:je,kb:ke) :: prs 
!      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa, qavorig
      real, dimension(ib:ie,jb:je,kb:ke) :: rho
!      real, dimension(ib:ie,jb:je,kb:ke) :: rhorig, rhnew
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u_start
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v_start
      real, dimension(ib:ie,jb:je,kb:ke) :: pi_start,prs_start,rho_start
      real, dimension(ib:ie,jb:je,kb:ke) :: thv_start,th_start
      real, dimension(ib:ie,jb:je,kb:ke) :: qv_start
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u_end
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v_end
      real, dimension(ib:ie,jb:je,kb:ke) :: pi_end,prs_end,rho_end
      real, dimension(ib:ie,jb:je,kb:ke) :: thv_end,th_end
      real, dimension(ib:ie,jb:je,kb:ke) :: qv_end

      double precision, dimension(ib:ie,jb:je,kb:ke) :: pi_tend1,pi_tend2
      double precision, dimension(ib:ie,jb:je,kb:ke) :: prs_tend1,prs_tend2
      double precision, dimension(ib:ie,jb:je,kb:ke) :: rho_tend1,rho_tend2
      double precision, dimension(ib:ie,jb:je,kb:ke) :: th_tend1,th_tend2
      double precision, dimension(ib:ie,jb:je,kb:ke) :: qv_tend1,qv_tend2

!      real, dimension(ib:ie,jb:je,kb:ke) :: pi_tend1,pi_tend2
!      real, dimension(ib:ie,jb:je,kb:ke) :: prs_tend1,prs_tend2
!      real, dimension(ib:ie,jb:je,kb:ke) :: rho_tend1,rho_tend2
!      real, dimension(ib:ie,jb:je,kb:ke) :: th_tend1,th_tend2
!      real, dimension(ib:ie,jb:je,kb:ke) :: qv_tend1,qv_tend2

      real, dimension(ib:ie,jb:je,kb:ke) :: pidiff,prsdiff,rhodiff
      real, dimension(ib:ie,jb:je,kb:ke) :: thdiff, qvdiff

      real, dimension(ib:ie,jb:je,kb:ke) :: zh
 

!-----------------------------------------------------------------------
! Vars related to reading in external soundings and base.F 
!-----------------------------------------------------------------------
      integer i,j,k,m,n,nsnd,kbot,ktop,tflag,nmax
      real zu,zv
      real z_trop,th_trop,th_sfc,t_trop,prs_sfc,qv_pbl,pi_sfc,t_sfc,rh_sfc,rh_pbl
      real qv_sfc,thv_sfc,psurf,tsurf,qsurf,thsurf,thvsurf,p_sfc
      real alpha,umax,nm,dudz,dvdz,rinum
      real, dimension(:), allocatable :: zsnd,thsnd,qvsnd,usnd,vsnd
      real, dimension(:), allocatable :: thvsnd,pisnd,psnd,tsnd,rhsnd
      real, dimension(:), allocatable :: thinterp,qvinterp,uinterp
      real, dimension(:), allocatable :: vinterp,pinterp,tinterp
      real, dimension(:), allocatable :: rhinterp
      integer :: kk,kup,kdn
      real :: interp_frac
!      real rslf,rsif

      real :: pisfc,the_sfc,thv1,thv2,pi1,pi2,z1,z2,p2,theq,qt_sfc
      real :: du,dv

      integer :: k1,k2,status


!------------------------------------------------------------------


      if(dowr) write(outfile,*) ' INSIDE BSS TEND CALCULATION'
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
! First determine "tend1" between original starting sounding and first
! sounding replacement.
! Read in the start and end profiles; wind and thermo separately,
! unless it's input_sounding file
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code here determines the starting
! environment (input_sounding)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Code below here is cribbed from base.F
      
        nmax = 100000 !should be enough!

        allocate(   zsnd(nmax) )
        allocate(  thsnd(nmax) )
        allocate(  qvsnd(nmax) )
        allocate(   usnd(nmax) )
        allocate(   vsnd(nmax) )
        allocate( thvsnd(nmax) )
        allocate(  pisnd(nmax) )
        allocate(   psnd(nmax) )
        allocate(   tsnd(nmax) )
        allocate(  rhsnd(nmax) )

        allocate( thinterp(nk) )
        allocate( qvinterp(nk) )
        allocate(  uinterp(nk) )
        allocate(  vinterp(nk) )
        allocate(  pinterp(nk) )
        allocate(  tinterp(nk) )
        allocate( rhinterp(nk) )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Trying updates to this based on some mpi calls...

      prc0:  &
      IF( myid.eq.0 )THEN
        ! only processor 0 does this:

        open(unit=40,file='input_sounding',status='unknown')

        ! read surface parameters:
        read(40,*) p_sfc, th_sfc, qv_sfc
!        if(dowr) write(outfile,*) 'Orig for BSS: p_sfc, th_sfc, qv_sfc = ', p_sfc, th_sfc, qv_sfc

        ! now, read entire sounding until end of file is discovered
        nsnd=1
        do k=1,nmax
            read(40,*,IOSTAT=status) zsnd(k+1),thsnd(k+1),qvsnd(k+1),usnd(k+1),vsnd(k+1)
            if(status /= 0) exit
!            if(dowr) write(outfile,*) 'Orig for BSS data:', zsnd(k+1),thsnd(k+1),qvsnd(k+1),usnd(k+1),vsnd(k+1)
          ! put qv in g/g
          qvsnd(k+1) = qvsnd(k+1)/1000.0
          nsnd=nsnd+1
        enddo
445     continue
        if(dowr) write(outfile,*) '  Found ',nsnd,'  levels (including surface)'
        if(dowr) write(outfile,*)
        close(unit=40)

        if( nsnd.le.2 )then
          print *
          print *,'  There is only one u,v level.  At least two u,v levels are needed. '
          print *
          print *,'  ..... stopping cm1 ..... '
          print *
          call stopcm1
        endif

      ENDIF  prc0

#ifdef MPI
        call MPI_BCAST(nsnd    ,  1 ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(p_sfc   ,  1 ,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(th_sfc  ,  1 ,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qv_sfc  ,  1 ,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(zsnd(1) ,nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(thsnd(1),nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(qvsnd(1),nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(usnd(1) ,nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(vsnd(1) ,nsnd,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
#endif

        if(imoist.eq.0) qv_sfc = 0.0

        p_sfc = p_sfc * 100.0
        ! put qv in g/g
        qv_sfc = qv_sfc / 1000.0
        pi_sfc  = (p_sfc/p00)**(rd/cp)
        thv_sfc = th_sfc*(1.0+qv_sfc*reps)/(1.0+qv_sfc)

        psurf  = p_sfc
        thsurf = th_sfc
        tsurf  = th_sfc * pi_sfc
        qsurf  = qv_sfc

        zsnd(1) = 0.0
        thsnd(1) = th_sfc
        qvsnd(1) = qv_sfc


        
!        open(unit=50,file='input_sounding',status='old')
!
!        ! read surface parameters:
!        read(50,*) p_sfc, th_sfc, qv_sfc
!        
!        p_sfc = p_sfc * 100.0
!        ! put qv in g/g
!        qv_sfc = qv_sfc / 1000.0
!        pi_sfc  = (p_sfc/p00)**(rd/cp)
!        thv_sfc = th_sfc*(1.0+qv_sfc*reps)/(1.0+qv_sfc)
!
!        psurf  = p_sfc
!        thsurf = th_sfc
!        tsurf  = th_sfc * pi_sfc
!        qsurf  = qv_sfc
!
!        zsnd(1) = 0.0
!        thsnd(1) = th_sfc
!        qvsnd(1) = qv_sfc
!
!        ! now, read entire sounding until end of file is discovered
!        nsnd=1
!        do k=1,nmax
!          read(50,*,end=345) zsnd(k+1),thsnd(k+1),qvsnd(k+1),usnd(k+1),vsnd(k+1)
!          ! put qv in g/g
!          qvsnd(k+1) = qvsnd(k+1)/1000.0
!          nsnd=nsnd+1
!        enddo
!345     continue
!
!        close(unit=50)

        if(imoist.eq.0) qvsnd = 0.0

! GHB had this commented out ("don't remember where this formulation comes from")
!!!        usnd(1) = 1.75*usnd(2)-usnd(3)+0.25*usnd(4)
!!!        vsnd(1) = 1.75*vsnd(2)-vsnd(3)+0.25*vsnd(4)

        ! 190324:  use simple extrapolation
        usnd(1) = usnd(2) - (zsnd(2)-zsnd(1))*(usnd(3)-usnd(2))/(zsnd(3)-zsnd(2))
        vsnd(1) = vsnd(2) - (zsnd(2)-zsnd(1))*(vsnd(3)-vsnd(2))/(zsnd(3)-zsnd(2))
        
!--------------------------------------------------------------------
!  Added by GHB, 061021:
!  Get thv and prs ... check if qv is too small.  If so, set rh to 5%
!  (This code has no effect on the sounding if qv > 1e-12 everywere)
!  (It was added to deal with the 0 g/kg qv values in the Trier sounding.)
        do k=1,nsnd
          thvsnd(k)=thsnd(k)*(1.0+reps*qvsnd(k))/(1.0+qvsnd(k))
        enddo
        pisnd(1)=pi_sfc
        do k=2,nsnd
          pisnd(k)=pisnd(k-1)-g*(zsnd(k)-zsnd(k-1))   &
                               /(cp*0.5*(thvsnd(k)+thvsnd(k-1)))
        enddo
        do k=1,nsnd
          psnd(k)=p00*(pisnd(k)**(cp/rd))
          tsnd(k)=thsnd(k)*pisnd(k)
        enddo
      if(imoist.eq.1)then
        do k=1,nsnd
          if(qvsnd(k).lt.1.0e-12)then
            qvsnd(k)=0.01*rslf(psnd(k),thsnd(k)*pisnd(k))
          endif
          rhsnd(k)=qvsnd(k)/rslf(psnd(k),tsnd(k))
        enddo
      endif
!--------------------------------------------------------------------
! interpolating sounding to grid

      DO j=jb,je
      DO i=ib,ie

        DO k=1,nk

            kk = 1
            do while( zsnd(kk) .lt. zh(i,j,k) )
              kk = kk+1
            enddo
            kdn = kk-1
            kup = kk

            zu=0.5*(zh(max(ib,i-1),j,k)+zh(i,j,k))
            interp_frac = (   zu        - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            uinterp(k) =  usnd(kdn) + ( usnd(kup)- usnd(kdn))*interp_frac

            zv=0.5*(zh(i,max(jb,j-1),k)+zh(i,j,k))
            interp_frac = (   zv        - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            vinterp(k) =  vsnd(kdn) + ( vsnd(kup)- vsnd(kdn))*interp_frac

!!!            ! if this is first grid point, utilize surface values of th,qv:
!!!            if( k.eq.1 ) kdn = 1

            interp_frac = (   zh(i,j,k) - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )

            thinterp(k) = thsnd(kdn) + (thsnd(kup)-thsnd(kdn))*interp_frac
            qvinterp(k) = qvsnd(kdn) + (qvsnd(kup)-qvsnd(kdn))*interp_frac
             pinterp(k) =  psnd(kdn) + ( psnd(kup)- psnd(kdn))*interp_frac
             tinterp(k) =  tsnd(kdn) + ( tsnd(kup)- tsnd(kdn))*interp_frac
            rhinterp(k) = rhsnd(kdn) + (rhsnd(kup)-rhsnd(kdn))*interp_frac

        ENDDO

!------------------------------------------------------
! convert read-in values to start (pre-BSS) values
!------------------------------------------------------

        do k=1,nk
           u_start(i,j,k) =  uinterp(k)
           v_start(i,j,k) =  vinterp(k)
          qv_start(i,j,k) = qvinterp(k)
          th_start(i,j,k) = thinterp(k)
         prs_start(i,j,k) =  pinterp(k)
        enddo

        ! get pi0 and prs0 from thv0, using hydrostatic equation

!        do k=1,nk
!          ! get qv from linear interpolation of rh:
!          qv_start(i,j,k) = rh_start(i,j,k)     &
!                            *rslf(prs_start(i,j,k),t_start(i,j,k))
!        enddo

        do k=1,nk
          thv_start(i,j,k)=th_start(i,j,k)*(1.0+reps   &
                           *qv_start(i,j,k))/(1.0+qv_start(i,j,k))
        enddo

          pi_start(i,j,1)=pi_sfc-g*zh(i,j,1)/(cp*0.5*    &
                         (thv_sfc+thv_start(i,j,1)))
        do k=2,nk
          pi_start(i,j,k)=pi_start(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))  &
                         /(cp*0.5*(thv_start(i,j,k)+thv_start(i,j,k-1)))
        enddo

        do k=1,nk
          prs_start(i,j,k)=p00*(pi_start(i,j,k)**(cp/rd))
        enddo

        do k=1,nk
           rho_start(i,j,k)=prs_start(i,j,k)/(rd*th_start(i,j,k)  &
                           *pi_start(i,j,k)*(1.0+qv_start(i,j,k)*reps))
        enddo

        ! rh, just in case we want/need it later

!      if(imoist.eq.1)then
!        do k=1,nk
!          rh_start(i,j,k)=qv_start(i,j,k)/   &
!                (rslf(prs_start(i,j,k),th_start(i,j,k)*pi_start(i,j,k)))
!        enddo
!      endif

      ENDDO    ! enddo for i loop
      ENDDO    ! enddo for j loop


!-----------------------------------------------------
! get end points for u- and v-winds
!-----------------------------------------------------


      DO j=jb,je

         do k=1,nk
           u_start(ie+1,j,k) =  u_start(ie,j,k)
         enddo

       ENDDO    ! enddo for j loop


      DO i=ib,ie

         do k=1,nk
           v_start(i,je+1,k) =  v_start(i,je,k)
         enddo

       ENDDO    ! enddo for i loop



        ! deallocate temporary 1D arrays
        deallocate(   zsnd )
        deallocate(  thsnd )
        deallocate(  qvsnd )
        deallocate(   usnd )
        deallocate(   vsnd )
        deallocate( thvsnd )
        deallocate(  pisnd )
        deallocate(   psnd )
        deallocate(   tsnd )
        deallocate(  rhsnd )
        deallocate( thinterp )
        deallocate( qvinterp )
        deallocate(  uinterp )
        deallocate(  vinterp )
        deallocate(  pinterp )
        deallocate(  tinterp )
        deallocate( rhinterp )



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code here determines the ending environment
! (wind_replace1 and/or thermo_replace1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (windchange.eq.1) then

        call read_wind('wind_replace1',u_end,v_end,zh)

      endif


      if (thermochange.eq.1) then

        call read_thermo('thermo_replace1',th_end,thv_end,qv_end, &
                          pi_end,prs_end,rho_end,zh)

      endif

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now, determine 'tend1' values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        total_time = sndg2time - bss_start


!---------------------------------------------------------------
! Start with wind profile tendencies
!---------------------------------------------------------------

      if (windchange.eq.1) then


       DO k=1,ke
       DO j=jb,je
       DO i=ib,ie+1

! Calculate wind profile differences

           udiff(i,j,k) = u_end(i,j,k) - u_start(i,j,k)
  
! Calculate wind profile tendency for each time step


           if(total_time.gt.0.0)then !make sure you don't divide by zero...
              
              u_tend1(i,j,k) = udiff(i,j,k)/total_time

           else

              u_tend1(i,j,k) = 0.0

           endif
 
        ENDDO
        ENDDO
        ENDDO


       DO k=1,ke
       DO j=jb,je+1
       DO i=ib,ie

! Calculate wind profile differences

           vdiff(i,j,k) = v_end(i,j,k) - v_start(i,j,k)

!     Calculate wind profile tendency for each time step

         if(total_time.gt.0.0)then !make sure you don't divide by zero...

            v_tend1(i,j,k) = vdiff(i,j,k)/total_time

         else

            v_tend1(i,j,k) = 0.0

         endif

        ENDDO
        ENDDO
        ENDDO

       end if


!---------------------------------------------------------------
! Next, thermo profile tendencies
!---------------------------------------------------------------

      if(thermochange.eq.1) then

       DO k=1,ke
       DO j=jb,je
       DO i=ib,ie
 

! Calculate thermo profile differences

           thdiff(i,j,k) = th_end(i,j,k) - th_start(i,j,k)
!           thvdiff(i,j,k) = thv_end(i,j,k) - thv_start(i,j,k)
           qvdiff(i,j,k) = qv_end(i,j,k) - qv_start(i,j,k)
           pidiff(i,j,k) = pi_end(i,j,k) - pi_start(i,j,k)
           prsdiff(i,j,k) = prs_end(i,j,k) - prs_start(i,j,k)
           rhodiff(i,j,k) = rho_end(i,j,k) - rho_start(i,j,k)

! Calculate thermo profile tendency for each time step

        if(total_time.gt.0.0)then !make sure you don't divide by zero... 

           th_tend1(i,j,k) = thdiff(i,j,k)/total_time
!           thv_tend(i,j,k) = thvdiff(i,j,k)/total_time
           qv_tend1(i,j,k) = qvdiff(i,j,k)/total_time
           pi_tend1(i,j,k) = pidiff(i,j,k)/total_time
           prs_tend1(i,j,k) = prsdiff(i,j,k)/total_time
           rho_tend1(i,j,k) = rhodiff(i,j,k)/total_time

        else

           th_tend1(i,j,k) = 0.0
!           thv_tend(i,j,k) = 0.0
           qv_tend1(i,j,k) = 0.0
           pi_tend1(i,j,k) = 0.0
           prs_tend1(i,j,k) = 0.0
           rho_tend1(i,j,k) = 0.0


        endif

        ENDDO
        ENDDO
        ENDDO

       
       end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
! Now determine "tend2" between first replacement sounding and second
! sounding replacement.
! Read in the start and end profiles; wind and thermo separately.
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      if (windchange.eq.1) then !read in previous base-state (wind only)

         call read_wind('wind_replace1',u_start,v_start,zh)

      end if  ! end if-check for reading in starting wind profile



      if (thermochange.eq.1) then !read in previous base-state (thermo only) 

       call read_thermo('thermo_replace1',th_start,thv_start,qv_start, &
                          pi_start,prs_start,rho_start,zh)

      end if  ! end if-check for reading in starting thermo profile


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code here determines the ending
! environment (read in external wind vs. 
! thermo files separately)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Read in second sounding (1st sounding is original base-state) as "end" (or new) environment to be incorporated

      if (windchange.eq.1) then !read in 1st final base-state (wind only)

         call read_wind('wind_replace2',u_end,v_end,zh)

      end if  ! end if-check for reading in final wind profile



      if (thermochange.eq.1) then !read in 2nd final base-state (thermo only) 

         call read_thermo('thermo_replace2',th_end,thv_end,qv_end,  &
                          pi_end,prs_end,rho_end,zh)

      end if  ! end if-check for reading in final thermo profile



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now, determine 'tend2' values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        total_time = bss_end - sndg2time


!---------------------------------------------------------------
! Start with wind profile tendencies
!---------------------------------------------------------------

      if (windchange.eq.1) then


       DO k=1,ke
       DO j=jb,je
       DO i=ib,ie+1

! Calculate wind profile differences

           udiff(i,j,k) = u_end(i,j,k) - u_start(i,j,k)
  
! Calculate wind profile tendency for each time step

           if(total_time.gt.0.0)then    !make sure you don't divide by zero...
              
              u_tend2(i,j,k) = udiff(i,j,k)/total_time

           else

              u_tend2(i,j,k) = 0.0

           endif
 
        ENDDO
        ENDDO
        ENDDO


       DO k=1,ke
       DO j=jb,je+1
       DO i=ib,ie

! Calculate wind profile differences

           vdiff(i,j,k) = v_end(i,j,k) - v_start(i,j,k)

! Calculate wind profile tendency for each time step

          if(total_time.gt.0.0)then    !make sure you don't divide by zero...
           
             v_tend2(i,j,k) = vdiff(i,j,k)/total_time

          else

             v_tend2(i,j,k) = 0.0

          endif

        ENDDO
        ENDDO
        ENDDO

       end if

 

!---------------------------------------------------------------
! Next, thermo profile tendencies
!---------------------------------------------------------------

      if(thermochange.eq.1) then

       DO k=1,ke
       DO j=jb,je
       DO i=ib,ie
 

! Calculate thermo profile differences

           thdiff(i,j,k) = th_end(i,j,k) - th_start(i,j,k)
!           thvdiff(i,j,k) = thv_end(i,j,k) - thv_start(i,j,k)
           qvdiff(i,j,k) = qv_end(i,j,k) - qv_start(i,j,k)
           pidiff(i,j,k) = pi_end(i,j,k) - pi_start(i,j,k)
           prsdiff(i,j,k) = prs_end(i,j,k) - prs_start(i,j,k)
           rhodiff(i,j,k) = rho_end(i,j,k) - rho_start(i,j,k)

! Calculate thermo profile tendency for each time step

           if(total_time.gt.0.0)then !make sure you don't divide by zero...
                            
             th_tend2(i,j,k) = thdiff(i,j,k)/total_time
!            thv_tend(i,j,k) = thvdiff(i,j,k)/total_time
             qv_tend2(i,j,k) = qvdiff(i,j,k)/total_time
             pi_tend2(i,j,k) = pidiff(i,j,k)/total_time
             prs_tend2(i,j,k) = prsdiff(i,j,k)/total_time
             rho_tend2(i,j,k) = rhodiff(i,j,k)/total_time

            else

             th_tend2(i,j,k) = 0.0
!            thv_tend(i,j,k) = 0.0
             qv_tend2(i,j,k) = 0.0
             pi_tend2(i,j,k) = 0.0
             prs_tend2(i,j,k) = 0.0
             rho_tend2(i,j,k) = 0.0

            endif
               

        ENDDO
        ENDDO
        ENDDO

       end if



!misc addition in GHB base.F       
       !--------------------------------------------------------------
!
!9000  print *
!      print *,'  9000: error reading BSS input_sounding file '
!      print *,'    ... stopping cm1 ... '
!      print *
!      call stopcm1
!
!9100  print *
!      print *,'  9100: error opening BSS input_sounding file '
!      print *,'    ... stopping cm1 ... '
!      print *
!      call stopcm1

    !--------------------------------------------------------------
 
      return
      end subroutine bss_tend




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read-in external file subrountines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_thermo(filename,theta,thetav,qv,piread,pres, &
                             density,zh)

      use input
      use constants
      use bc_module
      use comm_module
      use goddard_module, only : T0K,T00K,RT0
      use cm1libs , only : rslf,rsif
      use getcape_module

#ifdef MPI
      use mpi
#endif

      implicit none

      character (len=15) :: filename
      real, dimension(ib:ie,jb:je,kb:ke) :: piread,pres,density
      real, dimension(ib:ie,jb:je,kb:ke) :: thetav,theta
      real, dimension(ib:ie,jb:je,kb:ke) :: qv


      integer i,j,k,m,n,nsnd,kbot,ktop,tflag,nmax
      real zu,zv
      real th_sfc,prs_sfc,pi_sfc,t_sfc,rh_sfc
      real qv_sfc,thv_sfc,psurf,tsurf,qsurf,thsurf,thvsurf
      real p_sfc
      real alpha,umax,nm,dudz,dvdz,rinum
      real, dimension(:), allocatable :: zsnd,thsnd,qvsnd,usnd,vsnd,   &
                                         thvsnd,pisnd,psnd,tsnd,rhsnd
      real, dimension(:), allocatable :: thinterp,qvinterp, &
                                         pinterp,tinterp,rhinterp
      integer :: kk,kup,kdn
      real :: interp_frac
 !     real rslf,rsif

      real :: du,dv

      integer :: k1,k2
      real, dimension(ib:ie,jb:je,kb:ke) :: zh




! Code below here is cribbed from base.F

        nmax = 100000 !should be enough!

        allocate(   zsnd(nmax) )
        allocate(  thsnd(nmax) )
        allocate(  qvsnd(nmax) )
        allocate( thvsnd(nmax) )
        allocate(  pisnd(nmax) )
        allocate(   psnd(nmax) )
        allocate(   tsnd(nmax) )
        allocate(  rhsnd(nmax) )

        allocate( thinterp(nk) )
        allocate( qvinterp(nk) )
        allocate(  pinterp(nk) )
        allocate(  tinterp(nk) )
        allocate( rhinterp(nk) )

        open(unit=40,file=filename,status='old')

        ! read surface parameters:
        read(40,*) p_sfc, th_sfc, qv_sfc
        
        p_sfc = p_sfc * 100.0
        ! put qv in g/g
        qv_sfc = qv_sfc / 1000.0
        pi_sfc  = (p_sfc/p00)**(rd/cp)
        thv_sfc = th_sfc*(1.0+qv_sfc*reps)/(1.0+qv_sfc)

        psurf  = p_sfc
        thsurf = th_sfc
        tsurf  = th_sfc * pi_sfc
        qsurf  = qv_sfc

        zsnd(1) = 0.0
        thsnd(1) = th_sfc
        qvsnd(1) = qv_sfc

        ! now, read entire sounding until end of file is discovered
        nsnd=1
        do k=1,nmax
          read(40,*,end=445) zsnd(k+1),thsnd(k+1),qvsnd(k+1)
          ! put qv in g/g
          qvsnd(k+1) = qvsnd(k+1)/1000.0
          nsnd=nsnd+1
        enddo
445     continue

        close(unit=40)

        if(imoist.eq.0) qvsnd = 0.0
 

!--------------------------------------------------------------------
!  Added by GHB, 061021:
!  Get thv and prs ... check if qv is too small.  If so, set rh to 5%
!  (This code has no effect on the sounding if qv > 1e-12 everywere)
!  (It was added to deal with the 0 g/kg qv values in the Trier sounding.)
        do k=1,nsnd
          thvsnd(k)=thsnd(k)*(1.0+reps*qvsnd(k))/(1.0+qvsnd(k))
        enddo
        pisnd(1)=pi_sfc
        do k=2,nsnd
          pisnd(k)=pisnd(k-1)-g*(zsnd(k)-zsnd(k-1))   &
                               /(cp*0.5*(thvsnd(k)+thvsnd(k-1)))
        enddo
        do k=1,nsnd
          psnd(k)=p00*(pisnd(k)**(cp/rd))
          tsnd(k)=thsnd(k)*pisnd(k)
        enddo
      if(imoist.eq.1)then
        do k=1,nsnd
          if(qvsnd(k).lt.1.0e-12)then
            qvsnd(k)=0.01*rslf(psnd(k),thsnd(k)*pisnd(k))
          endif
          rhsnd(k)=qvsnd(k)/rslf(psnd(k),tsnd(k))
        enddo
      endif
!--------------------------------------------------------------------
! interpolating sounding to grid

      DO j=jb,je
      DO i=ib,ie

        DO k=1,nk

            kk = 1
            do while( zsnd(kk) .lt. zh(i,j,k) )
              kk = kk+1
            enddo

            kdn = kk-1
            kup = kk

!!!            ! if this is first grid point, utilize surface values of th,qv:
!!!            if( k.eq.1 ) kdn = 1

            interp_frac = (   zh(i,j,k) - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )

            thinterp(k) = thsnd(kdn) + (thsnd(kup)-thsnd(kdn))*interp_frac
            qvinterp(k) = qvsnd(kdn) + (qvsnd(kup)-qvsnd(kdn))*interp_frac
             pinterp(k) =  psnd(kdn) + ( psnd(kup)- psnd(kdn))*interp_frac
             tinterp(k) =  tsnd(kdn) + ( tsnd(kup)- tsnd(kdn))*interp_frac
            rhinterp(k) = rhsnd(kdn) + (rhsnd(kup)-rhsnd(kdn))*interp_frac


        ENDDO

!------------------------------------------------------
! convert read-in values to start (pre-BSS) values
!------------------------------------------------------

        do k=1,nk
            qv(i,j,k) = qvinterp(k)
         theta(i,j,k) = thinterp(k)
          pres(i,j,k) =  pinterp(k)
        enddo

        ! get pi0 and prs0 from thv0, using hydrostatic equation

!        do k=1,nk
!          ! get qv from linear interpolation of rh:
!          qv_start(i,j,k) = rh_start(i,j,k)     &
!                            *rslf(prs_start(i,j,k),t_start(i,j,k))
!        enddo

        do k=1,nk
          thetav(i,j,k)=theta(i,j,k)*(1.0+reps   &
                           *qv(i,j,k))/(1.0+qv(i,j,k))
        enddo

          piread(i,j,1)=pi_sfc-g*zh(i,j,1)/(cp*0.5*    &
                         (thv_sfc+thetav(i,j,1)))
        do k=2,nk
          piread(i,j,k)=piread(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))  &
                          /(cp*0.5*(thetav(i,j,k)+thetav(i,j,k-1)))
        enddo

        do k=1,nk
          pres(i,j,k)=p00*(piread(i,j,k)**(cp/rd))
        enddo

        do k=1,nk
           density(i,j,k)=pres(i,j,k)/(rd*theta(i,j,k)  &
                           *piread(i,j,k)*(1.0+qv(i,j,k)*reps))
        enddo

        ! rh, just in case we want/need it later

!      if(imoist.eq.1)then
!        do k=1,nk
!          rh_start(i,j,k)=qv_start(i,j,k)/   &
!                (rslf(prs_start(i,j,k),th_start(i,j,k)*pi_start(i,j,k)))
!        enddo
!      endif

      ENDDO    ! enddo for i loop
      ENDDO    ! enddo for j loop



        ! deallocate temporary 1D arrays
        deallocate(   zsnd )
        deallocate(  thsnd )
        deallocate(  qvsnd )
        deallocate( thvsnd )
        deallocate(  pisnd )
        deallocate(   psnd )
        deallocate(   tsnd )
        deallocate(  rhsnd )
        deallocate( thinterp )
        deallocate( qvinterp )
        deallocate(  pinterp )
        deallocate(  tinterp )
        deallocate( rhinterp )

     
      end subroutine read_thermo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      subroutine read_wind(filename,ucomp,vcomp,zh)

      use input
      use constants
      use bc_module
      use comm_module
      use goddard_module, only : T0K,T00K,RT0
      use cm1libs , only : rslf,rsif
      use getcape_module

#ifdef MPI
      use mpi
#endif

      implicit none

      character (len=13) :: filename



      integer i,j,k,m,n,nsnd,kbot,ktop,nmax
      real zu,zv
      real, dimension(:), allocatable :: zsnd,usnd,vsnd
      real, dimension(:), allocatable :: uinterp,vinterp
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ucomp
      real, dimension(ib:ie,jb:je+1,kb:ke) :: vcomp
      integer :: kk,kup,kdn
      real :: interp_frac

      integer :: flag,ttype
      real :: du,dv

      integer :: k1,k2
      real, dimension(ib:ie,jb:je,kb:ke) :: zh


! Code below here is cribbed from base.F

        nmax = 100000 !should be enough!

        allocate(   zsnd(nmax) )
        allocate(   usnd(nmax) )
        allocate(   vsnd(nmax) )
        allocate(  uinterp(nk) )
        allocate(  vinterp(nk) )

        ucomp = 0.0
        vcomp = 0.0


        open(unit=70,file=filename,status='old')

         ! now, read entire sounding until end of file is discovered
        nsnd=1
        do k=1,nmax
          read(70,*,end=745) zsnd(k+1),usnd(k+1),vsnd(k+1)
          nsnd=nsnd+1
        enddo
745     continue

        close(unit=70)

 
        usnd(1) = 1.75*usnd(2)-usnd(3)+0.25*usnd(4)
        vsnd(1) = 1.75*vsnd(2)-vsnd(3)+0.25*vsnd(4)

!--------------------------------------------------------------------
! interpolating sounding to grid

      DO j=jb,je
      DO i=ib,ie

        DO k=1,nk

            kk = 1
            do while( zsnd(kk) .lt. zh(i,j,k) )
              kk = kk+1
            enddo
            kdn = kk-1
            kup = kk

            zu=0.5*(zh(max(ib,i-1),j,k)+zh(i,j,k))

            interp_frac = (   zu        - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            uinterp(k) =  usnd(kdn) + ( usnd(kup)- usnd(kdn))*interp_frac

            zv=0.5*(zh(i,max(jb,j-1),k)+zh(i,j,k))
 
            interp_frac = (   zv        - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            vinterp(k) =  vsnd(kdn) + ( vsnd(kup)- vsnd(kdn))*interp_frac


        ENDDO

  
!------------------------------------------------------
! convert read-in values
!------------------------------------------------------

        do k=1,nk
           ucomp(i,j,k) =  uinterp(k)
           vcomp(i,j,k) =  vinterp(k)
        enddo

  
      ENDDO    ! enddo for i loop
      ENDDO    ! enddo for j loop



!-----------------------------------------------------
! get end points for u- and v-winds
!-----------------------------------------------------


      DO j=jb,je

         do k=1,nk
           ucomp(ie+1,j,k) =  ucomp(ie,j,k)
         enddo

       ENDDO    ! enddo for j loop


      DO i=ib,ie

         do k=1,nk
           vcomp(i,je+1,k) =  vcomp(i,je,k)
         enddo

       ENDDO    ! enddo for i loop

 


        ! deallocate temporary 1D arrays
        deallocate(   zsnd )
        deallocate(   usnd )
        deallocate(   vsnd )
        deallocate(  uinterp )
        deallocate(  vinterp )



      end subroutine read_wind

!
! ---------------------------------------------------------------------
! THIS FUNCTION CALCULATES THE LIQUID SATURATION VAPOR MIXING RATIO AS
! A FUNCTION OF TEMPERATURE AND PRESSURE
!
      REAL FUNCTION RSLF(P,T)

      use input, only : ptype
      use constants
      IMPLICIT NONE
      REAL ESL,X,T,P,C0,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (C0= .611583699E03)
      PARAMETER (C1= .444606896E02)
      PARAMETER (C2= .143177157E01)
      PARAMETER (C3= .264224321E-1)
      PARAMETER (C4= .299291081E-3)
      PARAMETER (C5= .203154182E-5)
      PARAMETER (C6= .702620698E-8)
      PARAMETER (C7= .379534310E-11)
      PARAMETER (C8=-.321582393E-13)

!  Note to self ... this should be changed, somehow, in the future.
!  GHB 060806

    if(ptype.eq.1.or.ptype.eq.2.or.ptype.eq.3.or.ptype.eq.5.or.ptype.eq.6)then

      ! from Bolton (1980, MWR)
      esl=611.2 * EXP( 17.67 * ( T  - 273.15 ) / ( T  - 29.65 ) )
      ! 171023 (fix for very cold temps):
      esl = min( esl , p*0.5 )
      rslf= eps * ESL /(P-ESL)

    elseif(ptype.eq.4)then

      rslf=380.00*exp(17.2693882-4097.8531/(t-35.86))/p

!    elseif(ptype.eq.3)then
!
!      X=MAX(-80.,T-273.16)
!      ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
!      RSLF=eps*ESL/(P-ESL)

    else

      esl=611.2 * EXP( 17.67 * ( T  - 273.15 ) / ( T  - 29.65 ) )
      ! 171023 (fix for very cold temps):
      esl = min( esl , p*0.5 )
      rslf= eps * ESL /(P-ESL)

    endif

      END FUNCTION RSLF

!
! ---------------------------------------------------------------------
! THIS FUNCTION CALCULATES THE ICE SATURATION VAPOR MIXING RATIO AS A
! FUNCTION OF TEMPERATURE AND PRESSURE
!
      REAL FUNCTION RSIF(P,T)

      use input, only : ptype
      use constants
      IMPLICIT NONE
      REAL ESI,X,T,P,C0,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (C0= .609868993E03)
      PARAMETER (C1= .499320233E02)
      PARAMETER (C2= .184672631E01)
      PARAMETER (C3= .402737184E-1)
      PARAMETER (C4= .565392987E-3)
      PARAMETER (C5= .521693933E-5)
      PARAMETER (C6= .307839583E-7)
      PARAMETER (C7= .105785160E-9)
      PARAMETER (C8= .161444444E-12)

!  Note to self ... this should be changed, somehow, in the future.
!  GHB 060806

    if(ptype.eq.1.or.ptype.eq.2.or.ptype.eq.3.or.ptype.eq.5.or.ptype.eq.6)then

      ! from Tao et al (1989, MWR)
      esi=611.2 * EXP( 21.8745584 * ( T  - 273.15 ) / ( T  - 7.66 ) )
      ! 171023 (fix for very cold temps):
      esi = min( esi , p*0.5 )
      rsif= eps * ESI /(P-ESI)

    elseif(ptype.eq.4)then

      rsif=380.00*exp(21.87455-5807.4743/(t-7.66))/p

!    elseif(ptype.eq.3)then
!
!      X=MAX(-80.,T-273.16)
!      ESI=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
!      RSIF=eps*ESI/(P-ESI)

    else

      esi=611.2 * EXP( 21.8745584 * ( T  - 273.15 ) / ( T  - 7.66 ) )
      ! 171023 (fix for very cold temps):
      esi = min( esi , p*0.5 )
      rsif= eps * ESI /(P-ESI)

    endif

      END FUNCTION RSIF



 END MODULE bss_tend_module
