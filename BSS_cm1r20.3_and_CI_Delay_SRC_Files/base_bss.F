 MODULE base_bss_module

! This subroutine is applying the base-state tendency (as determined in the
! bss_tend subroutine) after separating out the perturbation from the 
! current base-state. It also updates some other derived base-state variables.

! This is a sizeable update from the original BSS formulation (see Letkewicz et al.
! 2013 in Mon. Wea. Rev.)--April 2017 (see Davenport et al. 2019 in MWR)

      implicit none

      private
      public :: base_bss

   CONTAINS

      subroutine base_bss(bss_start,bss_end,u0,v0,ua,va,th0,qv0,thv0,  &
                      pi0,tha,qa,ppi,prs,prs0,rho,rho0,zh,c1,c2,       &
                      rho0s,pi0s,prs0s,rth0s,rth0,qc0,qi0,rr0,rf0,rrf0,&
                      mtime,sndg2time,sndg3time,windchange,            &
                      thermochange,dt,                                 &
                      upert,vpert,pipert,prspert,rhopert,              &
                      thpert,qvpert,                                   &
                      qavorig,tempc,rhorig,rhnew,rh0orig,              &
                      qvsorig,qvsiorig,satworig,                       &
                      satiorig,ssatworig,ssatiorig,                    &
                      u_tend1,u_tend2,v_tend1,v_tend2,pi_tend1,        &
                      pi_tend2, prs_tend1,prs_tend2,rho_tend1,         &
                      rho_tend2,th_tend1,th_tend2,qv_tend1,qv_tend2,   &
                      reqs_u,reqs_v,reqs_s,nw1,nw2,ne1,ne2,sw1,sw2,se1,&
                      se2,n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2,     &
                      uw31,uw32,ue31,ue32,us31,us32,un31,un32,         &
                      vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,         &
                      sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)


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

      real bss_start, bss_end  !BSS start/end times (in sec) 
      real sndg2time, sndg3time  !Model run time (since beginning) to implement 2nd, 3rd sndg 
                                 !NOTE: current formulation below says 3rd sndg is essentially 
                                 !when BSS ends (i.e., changes from orig to 2nd, 2nd to 3rd only)
      integer windchange, thermochange !check to change wind/thermo profile for BSS
!      integer num_incr  !number of increments (based on time step) that environ will be modified over

      double precision :: mtime

      double precision, dimension(ib:ie+1,jb:je,kb:ke) :: u_tend1,u_tend2
      double precision, dimension(ib:ie,jb:je+1,kb:ke) :: v_tend1,v_tend2

!      real, dimension(ib:ie+1,jb:je,kb:ke) :: u_tend1,u_tend2
!      real, dimension(ib:ie,jb:je+1,kb:ke) :: v_tend1,v_tend2

      real dt

      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
!      real, dimension(ib:ie+1,jb:je,kb:ke) :: udiff
!      real, dimension(ib:ie,jb:je+1,kb:ke) :: vdiff

      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,prs0,rho0,thv0,th0,rth0,qv0
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,tha,prs 
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa, qavorig
!     real, dimension(ib:ie,jb:je,kb:ke) :: rho, rh0, tempc
      real, dimension(ib:ie,jb:je,kb:ke) :: rho
      real, dimension(ib:ie,jb:je,kb:ke) :: tempc, rhorig, rhnew

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


      real, dimension(ib:ie+1,jb:je,kb:ke) :: upert
      real, dimension(ib:ie,jb:je+1,kb:ke) :: vpert
      real, dimension(ib:ie,jb:je,kb:ke) :: pipert,prspert,rhopert
      real, dimension(ib:ie,jb:je,kb:ke) :: thpert
      real, dimension(ib:ie,jb:je,kb:ke) :: qvpert


!saturation check variables
      real, dimension(ib:ie,jb:je,kb:ke) :: qvsorig,qvsiorig,satworig
      real, dimension(ib:ie,jb:je,kb:ke) :: satiorig,ssatworig,ssatiorig
!      real, dimension(ib:ie,jb:je,kb:ke) :: qvsnew,qvsinew,satwnew
!      real, dimension(ib:ie,jb:je,kb:ke) :: satinew,ssatwnew,ssatinew

!rh control variables 
      real, dimension(ib:ie,jb:je,kb:ke) :: rh0orig
!hydrostatic balance variables
!      real qt,thvnew(nk),pinew(nk)
      real, dimension(ib:ie,jb:je,kb:ke) :: pp3d

    
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,c1,c2
      real, dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s,rth0s
      real, dimension(ib:ie,jb:je,kb:ke) :: qc0,qi0,rr0,rf0,rrf0


      integer i,j,k,m,n
      real th_sfc,t_trop,prs_sfc,qv_pbl,pi_sfc,t_sfc,rh_sfc,rh_pbl
      real qv_sfc,thv_sfc,psurf,tsurf,qsurf,thsurf,thvsurf,p_sfc
 !     real, dimension(:), allocatable :: zsnd,thsnd,qvsnd,usnd,vsnd,   &
 !                                        thvsnd,pisnd,psnd,tsnd,rhsnd
 !     real, dimension(:), allocatable :: thinterp,qvinterp,uinterp,vinterp, &
 !                                        pinterp,tinterp,rhinterp
      integer :: kk,kup,kdn
      real :: interp_frac
 !     real rslf,rsif


 
!------------------------------------------------------------------
! Vars for updating boundary conditions
!------------------------------------------------------------------

      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_s
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,cmp,kmt+1) :: n3w1,n3w2,n3e1,n3e2,s3w1,s3w2,s3e1,s3e2
      real, intent(inout), dimension(cmp,jmp,kmp) :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp) :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp) :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp) :: ss31,ss32,sn31,sn32

!------------------------------------------------------------------

 !     if(dowr) write(outfile,*)'1)Inside BSS where BS changes are applied'
      if (mtime.eq.bss_start) write(outfile,*) '2)Inside BSS where BS changes are applied'




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
! First, separate out perturbation from the *current* base-state
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if(windchange.eq.1)then

!Separate perturbation from 3D field
      do k=kb,ke
        do j=jb,je+1
         do i=ib,ie
           vpert(i,j,k)=va(i,j,k)-v0(i,j,k)
         enddo
        enddo
      enddo

      do k=kb,ke
        do j=jb,je
         do i=ib,ie+1
           upert(i,j,k)=ua(i,j,k)-u0(i,j,k)
         enddo
        enddo
      enddo


      end if  



      if(thermochange.eq.1)then

!Separate perturbation from full 3D field


!QV
      do k=kb,ke
        do j=jb,je
         do i=ib,ie
           qvpert(i,j,k)=qa(i,j,k,nqv)-qv0(i,j,k)
         enddo
        enddo
      enddo



!prs
      do k=kb,ke
        do j=jb,je
         do i=ib,ie
           prspert(i,j,k)=prs(i,j,k)-prs0(i,j,k)
         enddo
        enddo
      enddo


!rho
      do k=kb,ke
        do j=jb,je
         do i=ib,ie
           rhopert(i,j,k)=rho(i,j,k)-rho0(i,j,k)
         enddo
        enddo
      enddo


      endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
! Second, apply computed tendency to the base-state and recombine to
! create a new 3D field
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! WIND MODIFICATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (windchange.eq.1) THEN

! Apply computed tendency to winds

      do k=kb,ke
        do j=jb,je+1
         do i=ib,ie

           if (mtime.le.sndg2time) then
             v0(i,j,k)=v0(i,j,k)+(v_tend1(i,j,k)*dt)
           elseif ((mtime.gt.sndg2time).and.(mtime.le.bss_end)) then
             v0(i,j,k)=v0(i,j,k)+(v_tend2(i,j,k)*dt)
           endif

         enddo
        enddo
      enddo

      do k=kb,ke
        do j=jb,je
         do i=ib,ie+1

           if (mtime.le.sndg2time) then
             u0(i,j,k)=u0(i,j,k)+(u_tend1(i,j,k)*dt)
           elseif ((mtime.gt.sndg2time).and.(mtime.le.bss_end)) then
             u0(i,j,k)=u0(i,j,k)+(u_tend2(i,j,k)*dt)
           endif

         enddo
        enddo
      enddo


!  Recombine 3D winds  

      do k=kb,ke
        do j=jb,je+1
         do i=ib,ie
           va(i,j,k)=v0(i,j,k)+vpert(i,j,k)
         enddo
        enddo
      enddo

      do k=kb,ke
        do j=jb,je
         do i=ib,ie+1
           ua(i,j,k)=u0(i,j,k)+upert(i,j,k)
         enddo
        enddo
      enddo

         ENDIF  !mod winds on/off
!!!!!!!!!!!!!!!!!!!!!!!!!! End Wind Modification !!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!THERMO MODIFICATION !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (thermochange.eq.1) THEN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!First, saturation check: compute original (current) state of saturation 
!!!cribbed from thompson.F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do k=kb,ke
        do j=jb,je
         do i=ib,ie

             tempc(i,j,k)=((th0(i,j,k)+tha(i,j,k)*(pi0(i,j,k)    &
                          +ppi(i,j,k))))-273.15
             qavorig(i,j,k,nqv)=qa(i,j,k,nqv)
             qvsorig(i,j,k)=rslf(prs(i,j,k),(th0(i,j,k)+tha(i,j,k))*   &
                            (pi0(i,j,k)+ppi(i,j,k)))

             if (tempc(i,j,k).le.0.0) then
              qvsiorig(i,j,k)=rsif(prs(i,j,k),(th0(i,j,k)+tha(i,j,k))* &
                              (pi0(i,j,k)+ppi(i,j,k)))
             else
              qvsiorig(i,j,k)=qvsorig(i,j,k)
             endif
             
             satworig(i,j,k)=qa(i,j,k,nqv)/qvsorig(i,j,k)
             satiorig(i,j,k)=qa(i,j,k,nqv)/qvsiorig(i,j,k)
             ssatworig(i,j,k)=satworig(i,j,k) - 1.
             ssatiorig(i,j,k)=satiorig(i,j,k) - 1.
     
         enddo
        enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!compute original relative-humidity   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do k=kb,ke
        do j=jb,je
         do i=ib,ie
          tempc(i,j,k)=((th0(i,j,k)+tha(i,j,k)*(pi0(i,j,k)    &
                        +ppi(i,j,k))))-273.15
          if(tempc(i,j,k).ge.0.0)then
! compute base state RH too
          rh0orig(i,j,k)=qv0(i,j,k)/(rslf(prs0(i,j,k),  &
                        th0(i,j,k)*pi0(i,j,k)))

           rhorig(i,j,k)=qa(i,j,k,nqv)/(rslf(prs(i,j,k),      &
                         (th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)  &
                          +ppi(i,j,k))))
          else
! compute base state RH too
          rh0orig(i,j,k)=qv0(i,j,k)/(rsif(prs0(i,j,k),  &
                        th0(i,j,k)*pi0(i,j,k)))
           rhorig(i,j,k)=qa(i,j,k,nqv)/(rsif(prs(i,j,k),      &
                         (th0(i,j,k)+tha(i,j,k))*(pi0(i,j,k)  &
                         +ppi(i,j,k))))
          endif
         enddo
        enddo
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply computed tendency to qv first
! Recombine qv field to check saturation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do k=kb,ke
        do j=jb,je
         do i=ib,ie

           if (mtime.le.sndg2time) then
             qv0(i,j,k)=qv0(i,j,k)+(qv_tend1(i,j,k)*dt)
             qa(i,j,k,nqv)=qv0(i,j,k)+qvpert(i,j,k)
           elseif ((mtime.gt.sndg2time).and.(mtime.le.bss_end)) then
             qv0(i,j,k)=qv0(i,j,k)+(qv_tend2(i,j,k)*dt) 
             qa(i,j,k,nqv)=qv0(i,j,k)+qvpert(i,j,k)
           endif

         enddo
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Okay, now apply tendency to other fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!potential temp
      do k=kb,ke
        do j=jb,je
         do i=ib,ie

           if (mtime.le.sndg2time) then
             th0(i,j,k)=th0(i,j,k)+(th_tend1(i,j,k)*dt)
           elseif ((mtime.gt.sndg2time).and.(mtime.le.bss_end)) then
             th0(i,j,k)=th0(i,j,k)+(th_tend2(i,j,k)*dt) 
           endif

         enddo
        enddo
      enddo

!non-dimensional pressure
      do k=kb,ke
        do j=jb,je
         do i=ib,ie

           if (mtime.le.sndg2time) then
             pi0(i,j,k)=pi0(i,j,k)+(pi_tend1(i,j,k)*dt)
           elseif ((mtime.gt.sndg2time).and.(mtime.le.bss_end)) then
             pi0(i,j,k)=pi0(i,j,k)+(pi_tend2(i,j,k)*dt) 
           endif
 
         enddo
        enddo
      enddo


!prs
      do k=kb,ke
        do j=jb,je
         do i=ib,ie

           if (mtime.le.sndg2time) then
             prs0(i,j,k)=prs0(i,j,k)+(prs_tend1(i,j,k)*dt)
           elseif ((mtime.gt.sndg2time).and.(mtime.le.bss_end)) then
             prs0(i,j,k)=prs0(i,j,k)+(prs_tend2(i,j,k)*dt) 
           endif

         enddo
        enddo
      enddo


!rho
      do k=kb,ke
        do j=jb,je
         do i=ib,ie

           if (mtime.le.sndg2time) then
             rho0(i,j,k)=rho0(i,j,k)+(rho_tend1(i,j,k)*dt)
           elseif ((mtime.gt.sndg2time).and.(mtime.le.bss_end)) then
             rho0(i,j,k)=rho0(i,j,k)+(rho_tend2(i,j,k)*dt) 
           endif

         enddo
        enddo
      enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Before recombining other fields, compute new relative-
!!!!humidity and difference from before.... 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do k=kb,ke
        do j=jb,je
         do i=ib,ie
          tempc(i,j,k)=((th0(i,j,k)+tha(i,j,k)*(pi0(i,j,k)    &
                        +ppi(i,j,k))))-273.15
          if(tempc(i,j,k).ge.0.0)then
           rhnew(i,j,k)=qa(i,j,k,nqv)/(rslf((prs0(i,j,k)+     &
                        prspert(i,j,k)),(th0(i,j,k)+          &
                        tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k))))
          else
           rhnew(i,j,k)=qa(i,j,k,nqv)/(rsif((prs0(i,j,k)+     &
                        prspert(i,j,k)),(th0(i,j,k)+          &
                        tha(i,j,k))*(pi0(i,j,k)+ppi(i,j,k))))
          endif
         enddo
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!adjust saturation so there isn't evap of cloud water/ice 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do k=kb,ke
        do j=jb,je
         do i=ib,ie
          tempc(i,j,k)=(th0(i,j,k)*pi0(i,j,k))-273.15

!adjust if gridpoint was saturated but now isn't
          if ((qa(i,j,k,nqc).gt.0.0).and.(rhnew(i,j,k).le.1.0))then

           if(tempc(i,j,k).ge.0.0)then

            qa(i,j,k,nqv)=qavorig(i,j,k,nqv)
            qv0(i,j,k)=qa(i,j,k,nqv)-qvpert(i,j,k)   !adjust base-state accordingly
            rho0(i,j,k)=prs0(i,j,k)/(rd*th0(i,j,k)*pi0(i,j,k)*(1.0+qv0(i,j,k)*reps))     !adjust base-state accordingly
            thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))   !adjust base-state accordingly
           endif

          endif

          if(tempc(i,j,k).lt.0.0)then

           if ((qa(i,j,k,nqi).gt.0.0).and.(rhnew(i,j,k).le.1.0))then

            qa(i,j,k,nqv)=qavorig(i,j,k,nqv)
            qv0(i,j,k)=qa(i,j,k,nqv)-qvpert(i,j,k)   !adjust base-state accordingly
            rho0(i,j,k)=prs0(i,j,k)/(rd*th0(i,j,k)*pi0(i,j,k)*(1.0+qv0(i,j,k)*reps))     !adjust base-state accordingly
            thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))   !adjust base-state accordingly
           endif
           if ((qa(i,j,k,nqs).gt.0.0).and.(rhnew(i,j,k).le.1.0))then

            qa(i,j,k,nqv)=qavorig(i,j,k,nqv)
            qv0(i,j,k)=qa(i,j,k,nqv)-qvpert(i,j,k)   !adjust base-state accordingly
            rho0(i,j,k)=prs0(i,j,k)/(rd*th0(i,j,k)*pi0(i,j,k)*(1.0+qv0(i,j,k)*reps))     !adjust base-state accordingly
            thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))   !adjust base-state accordingly
           endif

          endif

         enddo
        enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! Recombine 3D fields

!QV
      do k=kb,ke
        do j=jb,je
         do i=ib,ie
            qa(i,j,k,nqv)=qv0(i,j,k)+qvpert(i,j,k)
         enddo
        enddo
      enddo


!prs
      do k=kb,ke
        do j=jb,je
         do i=ib,ie
          prs(i,j,k)=prs0(i,j,k)+prspert(i,j,k)
         enddo
        enddo
      enddo


!rho
      do k=kb,ke
        do j=jb,je
         do i=ib,ie
           rho(i,j,k)=rho0(i,j,k)+rhopert(i,j,k)
         enddo
        enddo
      enddo



         ENDIF  !mod thermo on/off
!!!!!!!!!!!!!!!!!!!!!!!!!! End thermo Modification !!!!!!!!!!!!!!!!



!---------------------------------------------------------------------
! Not sure if needed, but this is stuff at end of base.F that also
! feeds into solve that's based on the base-state...


!------------------------------------------------------------------
!  thv0
!  Assumes th0, qv0 are accurate

      do j=jb,je
      do i=ib,ie

       do k=kb,ke
        if(imoist.eq.1)then
          thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))
        else
          thv0(i,j,k)=th0(i,j,k)
        endif
       enddo

      enddo
      enddo


!---------------------------------------------------------------------
! rr0, rf0, rho0s, rrf0

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rr0(i,j,k)=1.0/rho0(i,j,k)
      enddo
      enddo
      enddo

      do k=2,nk
      do j=jb,je
      do i=ib,ie
        rf0(i,j,k)=c1(i,j,k)*rho0(i,j,k-1)+c2(i,j,k)*rho0(i,j,k)
      enddo
      enddo
      enddo


      do j=jb,je
      do i=ib,ie
        ! cm1r17, 2nd-order extrapolation:
        rf0(i,j,1) = cgs1*rho0(i,j,1)+cgs2*rho0(i,j,2)+cgs3*rho0(i,j,3)
        rf0(i,j,0)=rf0(i,j,1)
        rho0s(i,j) = rf0(i,j,1)
        ! cm1r17, 2nd-order extrapolation:
        rf0(i,j,nk+1) = cgt1*rho0(i,j,nk)+cgt2*rho0(i,j,nk-1)+cgt3*rho0(i,j,nk-2)
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rrf0(i,j,k)=1.0/rf0(i,j,k)
      enddo
      enddo
      enddo


!-------------------------------------------------------------------------
!  values at surface: prs0s, pi0s, rth0s, rth0

      ! Get surface p/T/q (for surface models and for CAPE calculation):

        psurf = prs0(1,1,1)
        pi_sfc  = (psurf/p00)**(rd/cp)
        thsurf = th0(1,1,1)
        tsurf = thsurf*pi_sfc
        do j=jb,je
        do i=ib,ie
          prs0s(i,j) = psurf
          pi0s(i,j) = (psurf*rp00)**rovcp
          rth0s(i,j) = thsurf**(-1)
        enddo
        enddo

!-------------------------------------------------------------------------
! rth0

       do k=kb,ke
        do j=jb,je
        do i=ib,ie
          rth0(i,j,k)=1.0/th0(i,j,k)
        enddo
        enddo
        enddo



!--------------------------------------------------------------------------
! Update boundary conditions, as seen at end of base.F
!--------------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Fill in ghost cells

      call bcu(u0)
      call bcv(v0)

      !--------
#ifdef MPI
      call comm_3u_start(u0,uw31,uw32,ue31,ue32,   &
                            us31,us32,un31,un32,reqs_u)
      call comm_3u_end(u0,uw31,uw32,ue31,ue32,   &
                          us31,us32,un31,un32,reqs_u)
      call getcorneru3(u0,n3w1(1,1,1),n3w2(1,1,1),n3e1(1,1,1),n3e2(1,1,1),  &
                          s3w1(1,1,1),s3w2(1,1,1),s3e1(1,1,1),s3e2(1,1,1))
      call bcu2(u0)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=jb,je
      do i=ib,ie+1
        u0(i,j,0)    = cgs1*u0(i,j,1)+cgs2*u0(i,j,2)+cgs3*u0(i,j,3)
        u0(i,j,nk+1) = cgt1*u0(i,j,nk)+cgt2*u0(i,j,nk-1)+cgt3*u0(i,j,nk-2)
      enddo
      enddo
      !--------
#ifdef MPI
      call comm_3v_start(v0,vw31,vw32,ve31,ve32,   &
                            vs31,vs32,vn31,vn32,reqs_v)
      call comm_3v_end(v0,vw31,vw32,ve31,ve32,   &
                          vs31,vs32,vn31,vn32,reqs_v)
      call getcornerv3(v0,n3w1(1,1,1),n3w2(1,1,1),n3e1(1,1,1),n3e2(1,1,1),  &
                          s3w1(1,1,1),s3w2(1,1,1),s3e1(1,1,1),s3e2(1,1,1))
      call bcv2(v0)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=jb,je+1
      do i=ib,ie
        v0(i,j,0)    = cgs1*v0(i,j,1)+cgs2*v0(i,j,2)+cgs3*v0(i,j,3)
        v0(i,j,nk+1) = cgt1*v0(i,j,nk)+cgt2*v0(i,j,nk-1)+cgt3*v0(i,j,nk-2)
      enddo
      enddo




      if (mtime.eq.bss_end) print *, 'Leaving BSS'

      return
      end subroutine base_bss

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



 END MODULE base_bss_module




