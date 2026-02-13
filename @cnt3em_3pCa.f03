!*--------------------------------------------------------------------*
!  @cnt3ems_3pCa.f03                      Dec.24, 2016, Nov.1, 2024   !   
!                                                                     !
!  ## Molecular Dynamics in Relativistic Electromagnetic Fields ##    !
!                                                                     !
!   MPI 3: subroutine /forces/; real atomic mass mp/me=1836           !
!   READ_conf is used in nZ, nZA, intens, lambda                      !
!                                                                     !
!  Used files :                                                       !
!  1. @cnt3-3pCa.f03:  Molecular dynamics simulation code             !
!  2. param_em3p7_Ca.h: Common parameters of this simulation          !
!  3. Cntemp_config.STARTC: Define simulation time, box sizes,        !
!             atom names and pellet size, injected laser parameters   !
!  4. p_config_ss.xyz_D150, P135: Initially loaded particles of       !
!             H,C,Au atoms and electrons                              !             
!                                                                     !
!  > An explicit simulation code is strictly bound by the Courant     !
!    condition, Dx(length)/Dt(time) > c, the speed of light.          !
!  > The Gauss's law is used to find the product of current and       !
!    the longitudinal electric field.                                 !
!                                                                     !
!   Reference:                                                        !
!      Computer Physics Commun., vol.241, pp.56-63 (2019).            !
!      Fortran 2003 and MPI Packages   Nov. 3 (2020)                  !
!                                                                     !
!   Author/maintainer: Motohiko Tanaka, Ph.D.,Professor,              !
!      Graduate School of Chubu University, Kasugai 487-8501, Japan.  !
!                                                                     !
!*--------------------------------------------------------------------*
!                                                                     !
!   GPL-3.0 License, at https://github.com/Mtanaka77/                 !
!     /Relativistic_Molecular_Dynamics_Simulation/                    !
!                                                                     !
!  Initial version by Fujitsu FX100 Supercomputer, 2015-2019.         !
!    mpifrtpx -Kfast,openmp -o aq4.out @cnt3emq.f03 -lfftw3           !
!    (Ez,Bx) in envelop: (sin,cos)*exp(-(t/tq)**2)                    !
!                                                                     !
!*--------------------------------------------------------------------*
!                                                                     !
!   The CGS system:                                                   !
!      a_unit= 1.0000d-08 cm, Angstrom                                !
!      t_unit= 1.0000d-15 sec, ato second                             !
!      m_unit= 0.9110d-27 g, electron mass                            ! 
!      e_unit= 4.8032d-10 esu, unit charge                            !
!                                                                     !
!   In subroutine /READ_CONF/, nZ, nZA, intens, lambda are used.      !
!    3D filled H(+) by ranff(0.d0) ... in /moldyn/                      !
!    No gap between r>rout and r<rout  ... in /forces/                !
!                                                                     !
!   Subroutine descriptions                                           !
!   A) Program cnt3emp                                                !  
!        - /RUN_MD/ - /READ_CONF/                                     !
!                     formation of pellets (for kstart=0)             !
!                        or read read(12) at restart (for kstart=1)   !
!                     definition of several constants in /init/:      !
!                        constants and others                         !
!                    - /RUN_MD/ for start, and restart write(12)      !
!                    - /moldyn/: the main subroutine                  !
!                      final time /lplots/                            !
!                                                                     !
!   B) Inside of subroutine /moldyn/                                  ! 
!       initialization for it=1 case (kstart=0)                       !
!       Labels                                                        !
!       particle index: ncel(...)                                     !
!       initialization for writes: FT.13, FT.23                       !
!       start label 1000                                              !
!         FFTW initialization (one time at a start/restart time)      !
!         Maxwell equation of electric and magnetic fields            !
!         current density                                             !   
!         Poisson equation of electric field (longitudinal)           !
!         call subroutine /forces/                                    !  
!         advance position and momentum                               !
!         diagnosis (in every iwrt1 or iwrt2 steps)                   !
!       end of label 1000                                             !
!                                                                     !
!   C) Many subroutines including /edges/, /sendrecv1,2/, /filt3/...  !
!      when parallelization is utilized.                              !   
!                                                                     !   
!    Fortran 2003 /Fortran 2008 Write output (finished)               !   
!       write(11,'("This run uses ",i3," ranks",/)') size             !
!                                                                     !
!*--------------------------------------------------------------------*
! Parallel Fortran 2003
! $ mpif90 -mcmodel=medium -fpic -O2 @cnt3em_3pCa.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3 
!*--------------------------------------------------------------------
! 
      program  cnt3ems
      implicit  none
!
      include   'param_em3p7_Ca.h'
      include   'mpif.h'
!
      integer*4  size,rank,ierror,ipar,igrp,ifDebye
      real*8     ctime,wtime,ctime1,ctime2,walltime1,walltime2
      common/ps_time/ ctime1,ctime2,walltime1,walltime2 
!
      logical    if_start
! ---------------------------------------------
!**********************************************
!*  A global communication group.             *
!**********************************************
!  Fujitsu FX100: 52 node, 104 rank
!
      call mpi_init (ierror)
!
      call mpi_comm_rank (MPI_COMM_WORLD,rank,ierror)
      call mpi_comm_size (MPI_COMM_WORLD,size,ierror)
!
      ipar = 1 +rank  ! rank's name
      igrp = kgrp
!  
      if(ipar.eq.1) then   ! performs i/o if io_pe= 1.
        ionode= .true.
      else
        ionode= .false.
      end if
!
      suffix2= numbr2  ! Ca - new files
      suffix1= numbr1  ! Ca - old file
      suffix0= numbr0  ! C  - one character
!
      if(iflinx) then
        praefixs = '/home2/mtanaka/cntem3/'//sname
        praefixi = '/lv01/mtanaka/cntem3/'//cname  ! read(12) - common by NFS
        praefixc = '/lv01/mtanaka/cntem3/'//cname  ! write(13)
        praefixe = '/lv01/mtanaka/cntem3/'//sname  ! WRITE_CONF
      else
        praefixs = '/home/tanakam/cntem3A/'//sname      ! Simulator FX100
        praefixi = '/data/sht/tanakam/'//cname  ! read(12) - common by NFS
        praefixc = '/data/sht/tanakam/'//cname  ! write(13)
        praefixe = '/data/sht/tanakam/'//sname  ! WRITE_CONF
      end if
!
      ifDebye= 0   ! =1 for Debye-screening on
!
      if(kstart.eq.0) then
        if_start= .true.
      else
        if_start= .false.
      end if
!
!
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2, &
              status='unknown',form='formatted')
!
        write(11,100) size
  100   format('This run uses ',i3,' ranks',/)
!c      
        write(11,103) ipar,igrp
  103   format('My process is #',i3,' of group #',i3)
!
        close(11)
      end if
!
!*******************************************
!*  A group for intra-task communication.  *
!*******************************************
!
      call clocks (ctime1,walltime1,size)
! ---------------------------------------------------------------
      call RUN_MD (size,ipar,igrp,if_start,ifDebye)
! ---------------------------------------------------------------
      call clocks (ctime2,walltime2,size)
!
      ctime= ctime2 -ctime1
      wtime= walltime2 -walltime1
!
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,693) ipar,ctime,wtime
  693   format(/,'*ipar, ctime, wtime(sec)=',i3,1p2d15.7)
!
        close(11)
      end if
!*
      call mpi_finalize  (ierror)
!*
      stop
      end program cnt3ems
!
!
!------------------------------------------------------------------------
      subroutine date_and_time_7 (date_now,time_now)
!------------------------------------------------------------------------
      implicit none
!
      integer, dimension(8) :: ipresent_time
      character(len=8) :: time_now
      character(len=10) :: date_now

      call date_and_time(values=ipresent_time)

      write(time_now,'(i2,":",i2,":",i2)') ipresent_time(5:7)
      write(date_now,'(i4,"/",i2,"/",i2)') &
               ipresent_time(1),ipresent_time(2),ipresent_time(3)
!
      return
      end subroutine date_and_time_7
!
!
!------------------------------------------------------
      subroutine clocks (cputime,walltime,size)
!------------------------------------------------------
!*  Measure both cpu and elapsed times 
!
      implicit none
!
      include    'param_em3p7_Ca.h'
      include    'mpif.h'
!  
      real*8     cputime,walltime,buffer1(2),buffer2(2), &
                 cputime0,walltime0,ct(2),tm(2)
      integer*4  size,MPIerror
!  
      logical  first
      data     first/.true./
      save     first,cputime0,walltime0

!* Get the cpu and elapsed times on each node.

      if(iflinx) then
!       ct = etime(tm)   ! Pentium and Unix general.
!       cputime = tm(1)
!       buffer1(1) = cputime
!
        buffer1(2) = mpi_wtime()  ! in sec
        buffer1(1) = mpi_wtime()
      else
        buffer1(2) = mpi_wtime()  ! in sec
        buffer1(1) = buffer1(2)   ! SR16000
      end if

      if(first) then
        first= .false.
        cputime0 = buffer1(1)
        walltime0= buffer1(2)
      endif
!
      buffer1(2) = buffer1(2) -walltime0
!
!* The cpu and wall times, averaged over brother nodes.
!  size
      if(size.ne.1) then
        call mpi_allreduce (buffer1,buffer2,2,mpi_real8,   &
                            mpi_sum,MPI_COMM_WORLD,MPIerror)
        cputime  = buffer2(1)/size
        walltime = buffer2(2)/size
      else
        cputime  = buffer1(1)
        walltime = buffer1(2)
      end if
!
! When the date line is passed:
      if(cputime.lt.0.) cputime = cputime + 86400.d0

      return
      end subroutine clocks
!
!
!----------------------------------------------------------------------
      subroutine RUN_MD (size,ipar,igrp,if_start,ifDebye)
!----------------------------------------------------------------------
      implicit  none
!
      include    'param_em3p7_Ca.h'
      include    'mpif.h'
!
      integer*4  size,ipar,igrp,ifDebye
      logical    if_start
!
      integer*4  ns,np,nq,nCLp                     ! defined in l.330
      common/mainsp/ ns,np,nq,nCLp
!
      real*8     xg,yg,zg,px,py,pz,ch,am,ag        ! xg-zg in l.330
      common/mainda/ xg(npq0),yg(npq0),zg(npq0), & ! px-ag in /init/
                     px(npq0),py(npq0),pz(npq0), &
                     ch(npq0),am(npq0),ag(npq0)
      real*8     x0,y0,z0                          ! defied at l.1030
      common/initpos/ x0(npq0),y0(npq0),z0(npq0)
!
      integer*4  nZ,nZA,itab,nipl0,lipl0,itabs
      logical    cr_table
      real*8     fchar,fcharA,heavy,heavyA
      real*8     r_sp,d_sp,n_sp,lambda_d,massi,  &
                 ch_ion,wt_ion,rd_cp,rd_hp,      &
                 ch_el,wt_el,rd_el,rdebye
!
      common/iroha/ nZ,nZA                          ! in READ_CONF at l.316
      common/charg/ fchar,fcharA                    !
      common/hev_el/ heavy,heavyA                   ! defied at l.397
      common/ionsiz/ r_sp,d_sp,n_sp,massi,lambda_d, &  ! l.310
                     ch_ion,wt_ion,rd_cp,rd_hp,     &  ! l.557
                     ch_el,wt_el,rd_el
!
      common/srflst0/ cr_table,itab,nipl0(n00),lipl0(nbxs,n00)
      common/plupdat/ itabs
!
!  Single precision for the diagnostic routines.
!
      integer*4     i,j,k,it,is,istop,iwa,iwb,iwc,iwd,       &
                    ifrefl,nframe,ierror,ns1,np1,            &
                    nh1,nh2,nes1,nes2,nes3
      real*8        pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax, &
                    kJoule,kcal,mol,eV,                      &
                    tmax0,xx,yy,zz,r1,aau,ranff
      real*8        a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
      real*4        cptot,phi,tht,rgmax
!
      integer*4       iiwrt1,iiwrt2,iiwrt3,iiwrt4
      common/iiwrite/ iiwrt1,iiwrt2,iiwrt3,iiwrt4
!
      common/parm1/  it,is
      common/parm2/  pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm4/  phi,tht,rgmax
      common/physc/  a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
!
      common/parm9/  cptot
      common/abterm/ istop
      common/imemo/  iwa,iwb,iwc,iwd
!
      real*8         rcut_Clf,rcutlj,Temp,Temp_erg,epsCLJ,epsLJ, &
                     W_1p,Nele0
      common/cutoffrd/ rcut_Clf,rcutlj
      COMMON/ELSTA/  Temp,epsCLJ,epsLJ
      common/energ0/ W_1p,Nele0
!
      real*8          R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
      common/cntubes/ R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
      real*4        ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,          &
                    Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,            &
                    sex,sez,sbx,sbz,time
      common/ehist/ ekin(3000),ppot(3000),ekn1(3000),ekn2(3000), &
                    etot(3000),Rgyi(3000),Rgye(3000),Rgyf(3000), &
                    Rgyc(3000),Rion(3000),Nele(3000),dPot(3000), &
                    ecr(3000),elj(3000),sex(3000),sez(3000),     &
                    sbx(3000),sbz(3000),time(3000)
!+++
!     integer*4    mx,my,mz  ! main parameter statement
!
      integer*4    item3ab,item3
      real*8       Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
      real*8       etx,ety,etz,btx,bty,btz
      common/itre/ item3ab,item3
      common/hx2l/ Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
      common/feld/ etx(mx,my,mz),ety(mx,my,mz),etz(mx,my,mz),   &
                   btx(mx,my,mz),bty(mx,my,mz),btz(mx,my,mz) 
      integer*4    wrt8
!+++
!* added
      integer, dimension(8) :: ipresent_time
      character(len=8) :: time_now
      character(len=10) :: date_now
!
      CHARACTER*8    label,dummy*2
      COMMON/HEADR1/ label,date_now
!
      real*4         t,xp_leng
      COMMON/HEADR2/ t,xp_leng
!
      namelist/inp1/ phi,tht,cptot
      namelist/inp2/ dt,tmax,iiwrt1,iiwrt2,iiwrt3,iiwrt4
!
!
      label='CNT expl'
      call date_and_time_7 (date_now,time_now)
!
!**************************************************************
      xp_leng= 7.
      nframe= 4
!
      call FLOPEN (nframe,igrp)
!**************************************************************
!
      pi = 4.d0*datan(1.d0)
!* cgs
      a_unit= 1.0000d+00   ! cm
      m_unit= 0.9110d-27   ! g, electron mass
      e_unit= 4.8032d-10   ! esu, unit charge
      t_unit= 1.0000d+00   ! sec
!
!     ++++++++++++
!     sconf= 1.d-8  ! conversion from Angstrom to cm
!     ++++++++++++
      call READ_CONF (np,ifrefl)
!
! Carbon nanotube - str-cnt-mors.c: y-axis is a major axis
!   Use @dat2xyz.f to make CNT along z-axis
!
! 1040
!                 in Angstrom
! C    -8.30518   -83.83900   -41.48620
! C    -9.00426   -83.12850   -40.47340
!
      if(if_start) then
        if(ionode) then
!
!* CNT  60 cycles
        OPEN (unit=21,file='p_config_ss.xyz_D150',form='formatted')
!
        read(21,'(i6)') ns1
        read(21,'(a2)') dummy
!
        R_cnt1=  75.d-8  ! cm
        Z_cnt1= 153.d-8
!
        do i= 1,ns1
        read(21,'(a2,3f12.5)') dummy,xx,yy,zz
!
        xg(i)= sconv*xx
        yg(i)= sconv*yy
        zg(i)= sconv*zz   ! - long axis
        end do
!
        close(21)
!
!* Attach Au radially outside of each C
!
        aau= 2.0d-8  ! C-Au distance (cm)
!
        do i= 1,ns1
        r1= sqrt(xg(i)**2 +yg(i)**2)
        xg(ns1+i)= (r1 +aau)*xg(i)/r1
        yg(ns1+i)= (r1 +aau)*yg(i)/r1
        zg(ns1+i)= zg(i)
        end do
!
!       ---------
        ns= 2*ns1   ! C + Au
!       ---------
!
!* Two identical pellets
!    P135: 5060(half number); sqrt(2.) size of CNT length
!       OPEN (unit=21,file='p_config_ss.xyz_P063',form='formatted')
        OPEN (unit=21,file='p_config_ss.xyz_P135',form='formatted')
!
        read(21,'(i6)') np1
        read(21,'(a2)') dummy
!
        R_cnt2 = 30.0d-8
        Z_cnt2a= 10.0d-8   ! low of CNT
        Z_cnt2b= Z_cnt2a + 135.0d-8  ! 23 cycles
!
        do i= ns+1,ns+np1
        read(21,'(a2,3f12.5)') dummy,xx,yy,zz
!
        xg(i)= sconv*xx  ! cm
        yg(i)= sconv*yy
        zg(i)= Z_cnt2a +sconv*(zz +68.3d0)  ! - half length of long axis
!                                  +++++ Ang
        xg(i+np1)=  xg(i)
        yg(i+np1)=  yg(i)
        zg(i+np1)= -zg(i)
        end do
!
        close(21)
!
!       ---------
        np= 2*np1   ! two identical pellets
!       ---------
!
!* Associated electrons
!
!  READ_conf is used
!!      nZ = 1             ! heavy electrons for C(+Z)
!!      nZA= 4             ! heavy electrons for Au(+ZA/nZA)
        heavy= fchar/nZ    ! 6/1 in lump, used in /hev_el/
        heavyA= fcharA/nZA  ! 20/4, 40/4, 60/4
!
        nq = np + (nZ +nZA)*(ns/2)   ! ionized from H, C and Au
        nCLp= ns +np +nq   ! defined here
!       ++++++++++++++++
!
!-----------------
!* 2. Positions
!-----------------
!  It is important to generate positions of protons and electrons
!  in a mixed fashion so that their initial positions are close.
!
! Step 1: Generate protons and electrons in a pellet
!    Already done after read(21)
!
! Step 2: Electrons around H, carbon and Au
!
      do i= ns+np+1,ns+2*np
      xg(i)= xg(i-np) +0.1d-8*ranff(0.d0) !sgn(1,k)*(rd_CP +rd_el)
      yg(i)= yg(i-np) +0.1d-8*ranff(0.d0) !sgn(2,k)*(rd_CP +rd_el)
      zg(i)= zg(i-np) +0.1d-8*ranff(0.d0) !+sgn(3,k)*(rd_CP +rd_el)
      end do
        nes1= i
!
!
      i= ns +2*np
!* C
      do j= 1,ns/2 ! C(+Z)
      do k= 1,nZ   ! Heavy +(heavy) electrons per C(+Z)
      i= i +1
!
      xg(i)= xg(j) +0.1d-8*ranff(0.d0) !sgn(1,k)*(rd_CP +rd_el)
      yg(i)= yg(j) +0.1d-8*ranff(0.d0) !sgn(2,k)*(rd_CP +rd_el)
      zg(i)= zg(j) +0.1d-8*ranff(0.d0) !+sgn(3,k)*(rd_CP +rd_el)
      end do
      end do
        nes2= i
!
!* Au
      do j= ns/2+1,ns ! Au(+ZA)
      do k= 1,nZA  ! Heavy +(heavyA) electrons per Au
!             +++++
      i= i +1
!
      xg(i)= xg(j) +0.1d-8*ranff(0.d0) !sgn(1,k)*(rd_CP +rd_el)
      yg(i)= yg(j) +0.1d-8*ranff(0.d0) !sgn(2,k)*(rd_CP +rd_el)
      zg(i)= zg(j) +0.1d-8*ranff(0.d0) !sgn(3,k)*(rd_CP +rd_el)
      end do
      end do
        nes3= i
!
        if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'electron1=',nes1
        write(11,*) 'electron2=',nes2
        write(11,*) 'electron3=',nes3
        close(11)
        end if
!
!
!**
        OPEN (unit=22,file='ft22.data',form='formatted')
!
        do i= 1,ns/2
        write(22,'(a6,3f12.8)') ' C(cm)',xg(i),yg(i),zg(i)
        end do
!
        do i= ns/2+1,ns
        write(22,'(a6,3f12.8)') 'Au(cm)',xg(i),yg(i),zg(i)
        end do
!
        do i= ns+2*np+1,ns+2*np+nZ*(ns/2)
        write(22,'(a6,3f12.8)') 'elec 1',xg(i),yg(i),zg(i)
        end do
!
        do i= ns+2*np+nZ*(ns/2)+1,nCLp
        write(22,'(a6,3f12.8)') 'elec 2',xg(i),yg(i),zg(i)
        end do
!
        close (22)
!**
        end if
!
        call mpi_bcast (ns,1,mpi_integer,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (np,1,mpi_integer,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (nq,1,mpi_integer,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (nCLp,1,mpi_integer,0,MPI_COMM_WORLD,ierror)
!
!!      call mpi_bcast (fchar,1,mpi_real8,0,MPI_COMM_WORLD,ierror)
!!      call mpi_bcast (fcharA,1,mpi_real8,0,MPI_COMM_WORLD,ierror)
!!      call mpi_bcast (nZ,1,mpi_integer,0,MPI_COMM_WORLD,ierror)
!!      call mpi_bcast (nZA,1,mpi_integer,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (heavy,1,mpi_real8,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (heavyA,1,mpi_real8,0,MPI_COMM_WORLD,ierror)
!
        call mpi_bcast (R_cnt1,1,mpi_real8,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (Z_cnt1,1,mpi_real8,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (R_cnt2,1,mpi_real8,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (Z_cnt2a,1,mpi_real8,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (Z_cnt2b,1,mpi_real8,0,MPI_COMM_WORLD,ierror)
!
        call mpi_bcast (xg,npq0,mpi_real8,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (yg,npq0,mpi_real8,0,MPI_COMM_WORLD,ierror)
        call mpi_bcast (zg,npq0,mpi_real8,0,MPI_COMM_WORLD,ierror)
      else
!
!* Restart data.
!    These data overwrite definitions made above.
!
        if(ionode) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,630) 
  630     format(' Restart data are loaded from FT10.....',/)
          close(11)
        end if
!
!------------------------------------------
!* FT12 must be mounted on NSF volume.
!------------------------------------------
        OPEN (unit=12,file=praefixi//'.12'//suffix1,   &
                        status='old',form='unformatted')
!
        read(12) it,is,ns,np,nq,nCLp,itab,item3           ! nCLp
        read(12) xg,yg,zg,px,py,pz,ch,am,ag
        read(12) etx,ety,etz,btx,bty,btz
! 
        read(12) x0,y0,z0
        read(12) fchar,fcharA,heavy,heavyA,massi,nZ,nZA   ! nZ,nZA
        read(12) R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
        read(12) pi,tg,dt,dth,pthe,tmax0        ! dth; prefC_LJ,pref_LJ
        read(12) t,phi,tht,iiwrt1,iiwrt2,iiwrt3,iiwrt4
        read(12) ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,        &
                 Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,          &
                 sex,sez,sbx,sbz,time
!
        read(12) iwa,iwb,iwc,iwd
        read(12) r_sp,d_sp,n_sp,ch_ion,wt_ion,rd_cp,rd_hp,  &
                 ch_el,wt_el,rd_el 
        read(12) W_1p,Nele0
        close(12)
!
!          
!!      OPEN (unit=17,file=praefixi//'.17'//suffix2,        &
!!            status='old',form='unformatted')
!!                            +++
!!      read(17) nipl0,lipl0           ! is the same... 
!       close(17)
!                                        *****
!!      -----------------
!!      heavy= fchar/nZ     ! added to READ_CONF 
!!      heavyA= fcharA/nZA  !    1/01/2017
!!      nCLp= ns +np +nq    !
!!      -----------------
!
        if(ionode) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,640) t,it,is
  640     format(' Starting: t=',1pd11.3,'   it,is=',i7,i5)
!
          nCLp= ns +np +nq  ! just for write out here
          write(11,643) ns,np,nq,nCLp
  643     format(' ns,np,nq; nCLp=',i8,i6,i8,2x,i8,/)
!
          write(11,645) iwa,iwb,iwc,iwd
  645     format(' iwa, iwb, iwc, iwd=',4i7,/)
          close(11)
!          
!         OPEN (unit=50,file='ft50.data',               &
!               status='unknown',form='formatted')
!         write(50,640) t,it,is
!         write(50,643) ns,np,nq,nCLp
!         write(50,645) iwa,iwb,iwc,iwd
!         close(50)
        end if
!
!       call mpi_broadcast ()
      end if
!
!  -----------------------------
!
!     N_sp = (4.d0*pi/3.d0)* (1.d-8*R_sp)**3 * D_sp  ! R_sp in micron
!                ++
!     if(N_sp.lt.nq) then
!       if(ionode) 
!    *     write(11,*) 'Natural N_sp is used... N_sp, old nq=',
!    *                   N_sp,nq
!   ++ (nq,np) defined here ++
!       nq = N_sp                 ! do not subdivide (e,m)
!     end if
!  ------------------------
!     np= nq - fchar*ns
!  ------------------------
!
      c1= 2.9979d+10
      c2= 2.9979d+10**2
!
      ch_ion = e_unit
      wt_ion = massi*m_unit 
      rd_CP = 0.92d-8  ! cm, C(+Z)
      rd_HP = 0.50d-8  ! cm, H(+)
!
      ch_el = -e_unit
      wt_el =  m_unit
      rd_el = 0.50d-8  ! cm, e
!
!   statv= 2.998d+2 volt
!     pref_CL = e_unit**2/a_unit
!     pref_EM = e_unit*a_unit
!
      Temp_erg= Temp*1.6022d-12
      pthe = sqrt(Temp_erg)  ! cm/sec
!
!     kJoule  = 1.d10         ! erg
!     kcal= 4.1868d0 *kJoule  ! 4.18 J/cal
!     mol = 6.0220d23
      eV  = e_unit/299.79d0   ! eV to erg, 1.6022d-12
!
      pref_LJ = 48.0d0*epsLJ*eV
      prefC_LJ = 48.0d0*epsCLJ*eV
!***
      Lx3= xmax3 -xmin3  ! EM in cm
      Ly3= ymax3 -ymin3
      Lz3= zmax3 -zmin3
!
      hx= Lx3/mx
      hy= Ly3/my
      hz= Lz3/mz
!
      hx2= 2.d0*hx
      hy2= 2.d0*hy
      hz2= 2.d0*hz
! 
!   per statV/cm, per 1/(cm)^3
      p4dt= 4*pi*5.0d-19
      dV=  hx*hy*hz
      cdt= 2.9979d+10*5.0d-19
!***
      Lambda_D= sqrt(Temp_erg /(4*pi*D_sp*e_unit**2)) ! cm
!
!
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,600) label,date_now,time_now
  600   format(/,' ## Pellet in CNT (isolated system) ## ',a8,  &
               /,'   date=',a10,'  time=',a8,/)
!
        write(11,601) ns,fchar,fcharA
  601   format(' ns(C+Au)=',i7,'  fchar(C),fchar(Au)=',2f8.2,/)
!
        write(11,*) '&inp1  phi,tht.....'
        write(11,inp1)
!
!*  System size (-L/2, L/2)
!
        if(ifrefl.eq.0) then
          write(11,*) ' '
          write(11,*) ' Open cell boundary: mp=',massi,'.......'
          write(11,*) ' '
        else if(ifrefl.eq.2) then
          write(11,*) ' '
          write(11,*) ' Open/2-D closed boundary.........'
          write(11,*) ' '
        else
          write(11,*) ' '
          write(11,603) xmax3,ymax3,zmax3
  603     format('# Lx/2, Ly/2, Lz/2 (cm) = ',3f12.5,/)
        end if
!
        close(11)
      end if
!
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,605) dt,itabs
  605   format(' time step: dt=',1pd13.5,/,         &
               ' particle list update: itabs=',i4)
!
        write(11,607) ns,np,nq,fcharA,fchar
  607   format(' #carbon atoms =',i7,' (C+Au)',/,   &
               ' #protons      =',i7,' (H)   ',/,   &
               ' #electrons    =',i7,' (el)  ',//,  &
               ' Charge Au =',f8.2,'  (197*mp)',/,  &
               ' Charge C  =',f8.2,'  ( 12*mp)',/,  &
               ' Proton H  =    1.0   (    mp)',/)
!
        write(11,608) R_sp,D_sp,nq,ch_ion,ch_el,   &
                       wt_ion,wt_el,rd_CP,rd_HP,rd_el
  608   format(' R_sp(cm)=',1pd12.3,/,              &
               ' D_sp(/cm^3)=',d12.3,/,             &
               ' N_sp=',i10,/,                      &
               ' ',/,                               &
               ' ch_ion/e, ch_el/e=',2d15.7,/,      &
               ' wt_ion/m, wt_el/m=',2d15.7,/,      &
               ' rd_CP,rd_HP, rd_el/a=',3d15.7,/)
!
        write(11,609) R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
  609   format(' R_cnt1, Z_cnt1(cm)=',1p2d12.3,/,   &
               ' R_cnt2, Z_cnt2a,Z_cnt2b=',3d12.3,/)
!
        write(11,610) Temp,Temp_erg  !/Wrest
  610   format(' Temp(eV), Temp (erg)=',1p2d15.5,/)
!        
        write(11,611) rcut_Clf,prefC_LJ,pref_LJ,pthe,      &
                       p4dt,dV,p4dt/dV,cdt
  611   format('---------------------------------------',/, &
               ' Coulomb(short range): ch(i)*ch(j)/r2  ',/, &
               '   rcut_Clf=',1pd15.5,/,                    &
               '  ',/,                                      &
               '   prefC_LJ =',d15.5,/,                     &
               '   pref_LJ =',d15.5,/,                      &
               '  ',/,                                      &
               ' pthe (cm/s)=',d15.5,/,                     &
               ' ',/,                                       &
               ' p4dt =',d15.5,'   dV =',d15.5,/,           &
               ' p4dt/dV =',d15.5,/,                        &
               ' cdt =',d15.5,/,                            & 
               '---------------------------------------',/)
!
        write(11,613) mx,my,mz,Lx3,Ly3,Lz3,xmax3,xmin3,    &
                       ymax3,ymin3,zmax3,zmin3
  613   format('----------------------------------------------',/,  &
               ' 3D-EM field: mx, my, mz=',3i6,/,                   &
               '   with Lx3,Ly3,Lz3 (cm)=',3f13.8,/,                &
               '                                              ',/,  &
               '   Lx3: xmax3, xmin3=',2f13.8,' (cm)',/,            &
               '   Ly3: ymax3, ymin3=',2f13.8,/,                    &
               '   Lz3: zmax3, zmin3=',2f13.8,/,                    &
               '                                              ',/,  &
               '   special sinusoidal field: (Ez,Bx)          ',/,  &
               '     ak= 2*pi/lambda                          ',/,  &
               '     omega= ck= 2.3d+15 sec (for 0.8 micron)  ',/,  &
               '----------------------------------------------')
        write(11,*) '  '
        write(11,*) ' Eta= E0*sin(omg*t-k*y)*exp( ),'
        write(11,*) ' ffr( ): Eta*ch(i)'
        write(11,*) '  '
        write(11,*) '  Intensity: 6*10^17 W/cm^2 -> 2.1d+12 V/m (base)'
        write(11,*) '  Electric field: E0=7.d+7 statV/cm (base) '
        write(11,*) '  '
!
        close(11)
      end if
!--------------------------------------------------------------
!****************************************
!*  Prepare for graphic output.         *
!****************************************
!*  For 3-D plot of particles /pplt3d/.
!
      phi= -60.e0
      tht=  15.e0
!
      call ggauss 
!
!-----------------------------------------------------
!************************************
!*   Step 1 : Molecular dynamics.   *
!************************************
!
      istop= 0
      call init (xg,yg,zg,px,py,pz,ch,am,ag,        &
                 ipar,istop,if_start,ns,np,nq,nCLp)
!                **** ***** no igrp  ++ ++ ++ ++++
      if(istop.ne.0) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        if(ionode) write(11,*) ' Stop: call to init...'
        close(11)
!
        stop
      end if
!
      if(if_start) then
        if(ionode) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'ns,np,nq,nCLp=',ns,np,nq,nCLp
!
          write(11,*) '   '
          write(11,*) 'write: i ch am ag xg yg zg px py pz...'
          write(11,*) '# carbon atoms #'
          write(11,*) ' fchar(C),fchar(Au)=',fchar,fcharA
!  
          write(11,*) ' ch(esu), am(g), ag(cm)... '
          write(11,621) (i,ch(i),am(i),ag(i),xg(i),yg(i),zg(i),  &
                          px(i),py(i),pz(i),i=1,5)
          write(11,621) (i,ch(i),am(i),ag(i),xg(i),yg(i),zg(i),  &
                          px(i),py(i),pz(i),i=ns/2+1,ns/2+5)
!
          write(11,*) '   '
          write(11,*) '# protons #'
          write(11,621) (i,ch(i),am(i),ag(i),xg(i),yg(i),zg(i),  &
                          px(i),py(i),pz(i),i=ns+1,ns+5)
!
          nh1= ns+2*np+nZ*(ns/2)
          nh2= nh1 +nZA*(ns/2)
!
          write(11,*) '   '
          write(11,*) '# electrons #'
          write(11,621) (i,ch(i),am(i),ag(i),xg(i),yg(i),zg(i),  &
                          px(i),py(i),pz(i),i=ns+np+1,ns+np+5)
          write(11,621) (i,ch(i),am(i),ag(i),xg(i),yg(i),zg(i),  &
                          px(i),py(i),pz(i),i=ns+2*np+1,ns+2*np+5)
          write(11,621) (i,ch(i),am(i),ag(i),xg(i),yg(i),zg(i),  &
                          px(i),py(i),pz(i),i=nh1+1,nh1+5)
!
  621     format('i=',i6,1p3d10.2,1x,3d10.2,1x,3d8.1)
!
          close(11)
        end if
      end if
!
!
!************************************
!*   Step 1 : Molecular dynamics.   *
!************************************
!*-------------------------------------------------------------
      call moldyn (size,ipar,igrp,ifDebye,ifrefl,            &
                   if_start,ns,np,nq,nCLp)
!*---------------- ++ ++ ++ ++ +++ ----------------------------
!**************************************************************
!* Restart data.
!
      if(ionode) then
        OPEN (unit=12,file=praefixc//'.12'//suffix2,       &
                       status='replace',form='unformatted')
!
        write(12) it,is,ns,np,nq,nCLp,itab,item3
        write(12) xg,yg,zg,px,py,pz,ch,am,ag
        write(12) etx,ety,etz,btx,bty,btz
! 
        write(12) x0,y0,z0
        write(12) fchar,fcharA,heavy,heavyA,massi,nZ,nZA
        write(12) R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
        write(12) pi,tg,dt,dth,pthe,tmax
        write(12) t,phi,tht,iiwrt1,iiwrt2,iiwrt3,iiwrt4
        write(12) ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,        &
                  Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,          &
                  sex,sez,sbx,sbz,time
!
        write(12) iwa,iwb,iwc,iwd
        write(12) r_sp,d_sp,n_sp,ch_ion,wt_ion,rd_cp,rd_hp,  &
                  ch_el,wt_el,rd_el 
        write(12) W_1p,Nele0
        close(12)
!
!!      OPEN (unit=17,file=praefixi//'.17'//suffix2,        &
!!            status='unkown',form='unformatted')
!!                            ++++++
!!      write(17) nipl0,lipl0           ! is the same... 
!       close(17)
! 
!**********************************************************
!
        OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
!   ---------------------------          ++++++
        call lplots
!   ---------------------------
        call plote
        close(77)
      end if
!
      return
      end subroutine RUN_MD
!
!
!----------------------------------------------------------------------
      subroutine FLOPEN (nframe,igrp)
!----------------------------------------------------------------------
      implicit none
!
      include   'param_em3p7_Ca.h'
      integer*4 igrp,nframe
!
!      data suffix/'1','2','3','4','5','6','7','8','9','a'/
!      common/mldata/ suffix(10)
!     suffix2= '1'
!
      if(igrp.gt.10) then
        write(*,*) '*Stop File > 10...'
        stop
      end if
!
!   --------------------------
      if(.not.ionode) return
!   --------------------------
!c
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'igrp=',igrp,' uses ',                &
                        praefixs//'_config.START'//suffix0
        close(11)
!
        OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
                              status='unknown',form='formatted')
!                             ! just starts from here
        call gopen (nframe) 
        close (77)
      end if
!c 
      return
      end subroutine FLOPEN
!
!
!----------------------------------------------------------------------
      subroutine moldyn (size,ipar,igrp,ifDebye,ifrefl,              &
                         if_start,ns,np,nq,nCLp)
!------------------------++ ++ ++ ++ +++ ------------------------------
!*  Double precision.
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include 'fftw3.f03'     ! Generic FFTW3
!     include 'aslfftw3.f03'  ! NEC 2003 style
!
      include 'param_em3p7_Ca.h'
      include 'mpif.h'

      type(C_PTR), save :: plan,pinv
      integer, save ::  n_thread
      real(C_DOUBLE), dimension(mx,my,mz) :: qq
      real(C_DOUBLE), dimension(mx,my,mz), save :: psi
      real(C_DOUBLE), dimension(mx,my,mz) :: qq_c,psi_c 
!+++++
!
      integer*4  size,ipar,igrp,ifrefl,ifDebye,  &
                 ns,np,nq,nCLp,n_twice   ! for restart
      logical    if_start                ! for restart
!
      real*8     xg,yg,zg,px,py,pz,ch,am,ag 
      common/mainda/ xg(npq0),yg(npq0),zg(npq0), &
                     px(npq0),py(npq0),pz(npq0), &
                     ch(npq0),am(npq0),ag(npq0)
      real*8     vx(npq0),vy(npq0),vz(npq0)
!
      real*8     ffr(npq0,3)   !! npq0*3
      real*8     x0,y0,z0
      common/initpos/ x0(npq0),y0(npq0),z0(npq0)
!
      real*8     fchar,fcharA,heavy,heavyA
      integer*4  nZ,nZA
      common/charg/ fchar,fcharA
      common/hev_el/ heavy,heavyA
      common/iroha/ nZ,nZA
!
      integer*4       iiwrt1,iiwrt2,iiwrt3,iiwrt4,iiwrt0
      common/iiwrite/ iiwrt1,iiwrt2,iiwrt3,iiwrt4
!
      integer*4   npio
      parameter  (npio=ns0+np0+nq0)
      real*4    x4(npio),y4(npio),z4(npio),  &
                ch4(npio),am4(npio),ag4(npio)
!
      real*4         t,teql,xp_leng,Rgy0,Rgy1,Rgy2,Rgy3
      CHARACTER*8    label,date_now*10
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ t,xp_leng
!
      real*4        ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,          &
                    Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,            &
                    sex,sez,sbx,sbz,time,ekin20(3000*19)
      common/ehist/ ekin(3000),ppot(3000),ekn1(3000),ekn2(3000), &
                    etot(3000),Rgyi(3000),Rgye(3000),Rgyf(3000), &
                    Rgyc(3000),Rion(3000),Nele(3000),dPot(3000), &
                    ecr(3000),elj(3000),sex(3000),sez(3000),     &
                    sbx(3000),sbz(3000),time(3000)
      equivalence  (ekin(1),ekin20(1))
!
      integer*4     i,j,k,kk,jj,ibox,neigh,it,is,iwa,iwb,iwc,iwd,    &
                    iwrt1,iwrt2,iwrt3,iwrt4,istop,iwrta,iwrtb,iwrtc, &
                    iwrtd,ix,iy,iz,nskip,nsk,ncoe,wrt2

      real*8        pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax, &
                    a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest,     &
                    r2,rcut2,                                    &
                    E_C_s,E_C_PME,E_C_r,E_LJ,E_elas,             &
                    E_C_r1,E_C_r2,E_LJ2,                         &
                    vxav,vyav,vzav,rr1,rr0,ph0,ranff
      real*4        phi,tht,cptot,fchar4,fcharA4,Temp4,rgmax,    & 
                    xmax4,ymax4,zmax4,dx,dy,dz,                  &
                    vm,s0,s1,s2,rcore,rr,r1,ani,ane,             &
                    vmax2,dgaus2
      real*8        cpu0,cpu1,dcpu,walltm0,walltm1
!  
      real*8          ctime1,ctime2,walltime1,walltime2
      common/ps_time/ ctime1,ctime2,walltime1,walltime2 
!
      common/parm1/ it,is
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm4/ phi,tht,rgmax
      common/physc/ a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
!
      common/parm9/ cptot
      common/imemo/ iwa,iwb,iwc,iwd
      common/iotim/ iwrt1,iwrt2,iwrt3,iwrt4
      common/abterm/ istop
      COMMON/ENERGY/ E_C_s,E_C_PME,E_C_r,E_LJ,E_elas
!
      real*8         Temp,epsCLJ,epsLJ,m_gamma,Pot0,W_1p,Nele0,   &
                     gamma,gam1,gam2,gax,gay,gaz,fnml
      real*8         R_sp,D_sp,N_sp,Lambda_D,massi,               &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el, &
                     R_cn1,R_cn2,Z_cn,Z_cn1,Z_cn2,rcut_Clf,rcutlj
!
      real*8          R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
      common/cntubes/ R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
!
      COMMON/ELSTA/  Temp,epsCLJ,epsLJ
      common/energ0/ W_1p,Nele0
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,              &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
      common/cutoffrd/ rcut_Clf,rcutlj
!
      integer*4      ncel,lcel,ibind
      common/srflist/ ncel(nc3),lcel(nbxc,nc3)
      COMMON/boxind/ ibind(27,nc3)
!
      logical        cr_table
      integer*4      itab,nipl0,lipl0
      common/srflst0/ cr_table,itab,nipl0(n00),lipl0(nbxs,n00)
!
      real*8     intens,lambda,E0,ak,ct,omega,Eta,Bta,tp,taup,  &
                 ff,p_xyz,yg0,yg02,xz0,xz02
      common/electr/ intens,lambda
      common/swaves/ E0,ak,ct,omega
!
!+++++
!!    integer*4  mx,my,mz       ! in parameter statement
      integer*4  l,m,n,ll,mm,nn !,mxh,myh,mzh
!     parameter  (mxh=mx,myh=my,mzh=mz)
!     parameter  (mxh=mx/2+1,myh=my/2+1,mzh=mz/2+1)
!
      integer*4  ifeh,item3ab,item3,npar3,ierror,wrt7
      real*8     Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
      real*8     xx0,yy0,zz0,etx,ety,etz,btx,bty,btz,       &
                 EE(mx,my,mz),HH(mx,my,mz),           &
                 xx,yy,zz,etxi,etyi,etzi,btxi,btyi,btzi,    &
                 prho,akx,aky,akz
      common/itre/ item3ab,item3
      common/hx2l/ Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
      common/feld/ etx(mx,my,mz),ety(mx,my,mz),etz(mx,my,mz), &
                   btx(mx,my,mz),bty(mx,my,mz),btz(mx,my,mz)
!
      real*8    etx1(mx,my,mz),ety1(mx,my,mz),etz1(mx,my,mz), &
                cjx(mx,my,mz),cjy(mx,my,mz),cjz(mx,my,mz),    &
                dd(3),ta(3),r,rscCL,rCL,ccj,pp1,pp2,pp3,      &
                chvx,chvy,chvz
!
      real*4    setx,setz,sbtx,sbtz,psi2,delV
      real*4    eetx(mxh,myh,mzh),eety(mxh,myh,mzh),eetz(mxh,myh,mzh), &
                bbtx(mxh,myh,mzh),bbty(mxh,myh,mzh),bbtz(mxh,myh,mzh), &
                qq1(mxh,myh,mzh),psi1(mxh,myh,mzh)  
      real*8    cp1,cp2,cp3,cp4,cp5,cp6,cp7,wall1,wall2,wall3,         &
                wall4,wall5,wall6,wall7
!+++++
      logical    first11,first06,eq_phase,ifeqlib,    &
                 ifskip_e,ifskip_p,first_kk,first_es, &
                 read_10
      data       first11/.true./,first06/.true./,first_kk/.true./, &
                 first_es/.true./,read_10/.true./
!
!
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'moldyn: Round-robin partitioning is used...' 
        write(11,*) ' '
        close(11)
      end if
!
!--------------------------
!*  Initial condition.
!--------------------------
!
      iiwrt0= iiwrt1
!
      if(if_start) then
        tg= 0.
!
        is= 0
        it= 0
        iwa= -1
        iwb= -1
        iwc= -1
        iwd= -1
!
        do i= 1,nCLp
        x0(i)= xg(i)
        y0(i)= yg(i)
        z0(i)= zg(i)
        end do
!
        eq_phase= .true.  ! True during equilibration
        cr_table= .true.  ! create p-table
        n_twice= 0        ! just go
!
        itab= 0   
        item3= 0   
!
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        etx(l,m,n)= 0
        ety(l,m,n)= 0
        etz(l,m,n)= 0
        btx(l,m,n)= 0
        bty(l,m,n)= 0
        btz(l,m,n)= 0
        end do
        end do
        end do
!
        do i= 1,3000*19
        ekin20(i)= 0
        end do
      else
!
!  In restart 
        eq_phase= .false.
        cr_table= .true. 
        n_twice= 0  ! -1 ! restart by two steps is not necessary 
!
        if(ionode) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'write x0,y0,z0...'
          write(11,132) (i,x0(i),y0(i),z0(i),i=1,5)
          write(11,132) (i,x0(i),y0(i),z0(i),i=ns+1,ns+5)
          write(11,132) (i,x0(i),y0(i),z0(i),i=ns+np+1,ns+np+5)
  132     format('i=',i6,2x,1p3d11.4)
!
          close(11)
        end if
      end if
!
!---------------------------------------------------------
!* Control particle list update both at t=0 and restart
!---------------------------------------------------------
!  Lx3 in common
!
      dth= 0.5d0*dt
      rcut2 = rcut_Clf**2  ! Cutoff, cm**2 
!
!   -----------------
      call Labels
!   -----------------
!
      do k= 1,nc3
      ncel(k)= 0             ! Clear the cell register.
      end do
!
!* All particles
!
      do 100 i = 1,nCLp
      ix= isizeX*(x0(i)-xmin3)/Lx3 +1.0001
      iy= isizeY*(y0(i)-ymin3)/Ly3 +1.0001
      iz= isizeZ*(z0(i)-zmin3)/Lz3 +1.0001
!
      if(ix.le.1 .or. ix.ge.isizeX) go to 100
      if(iy.le.1 .or. iy.ge.isizeY) go to 100
      if(iz.le.1 .or. iz.ge.isizeZ) go to 100
!
      ibox = ix + isizeX*(iy-1 + isizeY*(iz-1))
!
      ncel(ibox)= ncel(ibox) +1
      lcel(ncel(ibox),ibox)= i
  100 continue
!
!
      kk= 0
      do 200 i= ipar,nCLp,size
      kk= kk +1
!
      nipl0(kk)= 0
!
      ix= isizeX*(x0(i)-xmin3)/Lx3 +1.0001
      iy= isizeY*(y0(i)-ymin3)/Ly3 +1.0001
      iz= isizeZ*(z0(i)-zmin3)/Lz3 +1.0001
!
      if(ix.le.1 .or. ix.ge.isizeX) go to 200
      if(iy.le.1 .or. iy.ge.isizeY) go to 200
      if(iz.le.1 .or. iz.ge.isizeZ) go to 200
!
      ibox = ix + isizeX*(iy-1 + isizeY*(iz-1))
!
!* Particles other than pellet are counted here
!
      do 230 k = 1,27
      neigh= ibind(k,ibox)
      if(neigh.eq.0) go to 230     ! Far away
!
      do 240 jj= 1,ncel(neigh)      ! Find ions in the boxes around i-th
      j= lcel(jj,neigh)             !  j-th belongs to this box.
      if(j.le.i) go to 240
!        ++++++
!
      dx= x0(i) -x0(j)  ! no folding
      dy= y0(i) -y0(j)
      dz= z0(i) -z0(j)
!
      r2 = dx**2 + dy**2 + dz**2
      if(r2.lt.rcut2) then
        nipl0(kk)= nipl0(kk) +1
        lipl0(nipl0(kk),kk)= j
      end if
  240 continue
  230 continue
  200 continue
! ++++++ initial +++++ at restart ++++++++++++
!
!-------------------------------------------------------
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,*) 'Preparation for write(13)...'
        close(11)
!
        fchar4 = fchar
        fcharA4 = fcharA
        Temp4  = Temp  
        xmax4  = xmax3
        ymax4  = ymax3
        zmax4  = zmax3
!
        OPEN (unit=13,file=praefixc//'.13'//suffix2,     &
              status='unknown',form='unformatted')
!
        write(13) ns,np,nq,fchar4,fcharA4,Temp4,xmax4,ymax4,zmax4
!
        do 7 i= 1,npio
        ch4(i)= ch(i)
        am4(i)= am(i)
        ag4(i)= ag(i)
    7   continue
!
        write(13) ch4,am4,ag4
        close(13)
      end if
!
!***
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Preparation for write(29)...'
        write(11,*) '  mxh,myh,mzh=',mxh,myh,mzh 
        close(11)
!
        OPEN (unit=29,file=praefixc//'.29'//suffix2,     &
              status='unknown',form='unformatted')
!
!  Half size (mxh,myh,mzh)
!
        write(29) mxh,myh,mzh
        write(29) Lx3,Ly3,Lz3
        write(29) xmax3,ymax3,zmax3,xmin3,ymin3,zmin3
        close(29)
      end if
!
!***
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Preparation for write(30)...'
        close(11)
!
        OPEN (unit=30,file=praefixc//'.30'//suffix2,     &
              status='unknown',form='unformatted')
!
        write(30) mxh,myh,mzh
        write(30) Lx3,Ly3,Lz3
        write(30) xmax3,ymax3,zmax3,xmin3,ymin3,zmin3
        close(30)
      end if
!
!* To avoid appending to a stale file
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Preparation for write(23)...'
        close(11)
!
        OPEN (unit=23,file=praefixc//'.23'//suffix2,    &
              status='unknown',form='unformatted')
!
        write(23) ns,np,nq
        write(23) fchar,fcharA,massi
        write(23) am,ch
!
        close(23)
      end if
!
      call clocks (cpu0,walltm0,size)  ! initiate clock (sec)
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'before <moldyn> (sec)=',cpu0 -ctime1
        close(11)
      end if
!  
!-------------------------------------------------------
!  dt= 0.0005e-15 -> EM: dx/dt= 4.d10 > 3.d10  Courant cond. safe!
!        (100,200,200), (5.,10.,10.) 10^2 Ang^3 -> 5.0 Ang
!
 1000 continue
      call clocks (cpu1,walltm1,size)
      dcpu= cpu1 -cpu0 ! sec
!
      if(tg.gt.tmax) go to 2000    
!
!  The mod(it,5)=0 step was executed in the previous run, so 
!  a new step should start at mod(it,5)=1 !
!
      if((dcpu/60.d0).gt.cptot .and. mod(it,5).eq.0) then
        if(ionode) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) '*Normal stop: previous step it=',it
          close(11)
        end if
!
        go to 2000
      end if
!
      if(istop.ge.1) then 
        wrt2= 50 +ipar
        write(wrt2,*) 'Abnormal termination istop=',istop
        write(wrt2,*) '  ipar,t=',ipar,tg
        go to 2000
      end if
!
      it= it + 1
      tg= tg + dt
!
      if(it.eq.1) tg= 0.d0  ! sec
      t= tg                 ! for iiwrt1, write(23)...
!  
!-------------------------------
!  Define write-out intervals.
!-------------------------------
!
!    *********************
      item3= item3 +1        ! starts at item3= 1
      n_twice= n_twice +1    !           n_twice= 1
!    *********************
!
      teql= 0.0005d-15 
      ifeqlib= (if_start).and.(t.lt.teql)
!               ++++++++
!
      if(ifeqlib) then
        iiwrt1= 1000 ! ??
      else
        iiwrt1= iiwrt0
      end if
!
      iwrt1= mod(it,iiwrt1)   ! iiwrt1= 1,2,..,iiwrt1-1,0   ! time 
      iwrt2= mod(it,iiwrt2)   ! iwrtb(t,dtwr2)    ! general field plot
      iwrt3= mod(it,iiwrt3)   ! iwrtc(t,dtwr3)    ! write(29), write(30)
      iwrt4= mod(it,iiwrt4)   ! iwrtd(t,1.0e-15)  ! plots/restart data
!*
      if(iwrt1.eq.0) then
        if(.not.ifeqlib) then
          is= is +1
!
          if(is.ge.3000) then
            call rehist
          end if
        else
          is= 1
        end if
      end if
!
!     is_yes= .false.                                  ! no
!     if(iwrt1.eq.0 .and..not.ifeqlib) is_yes= .true.  ! yes
!
!-------------------------
!*   Velocity update.
!-------------------------
!---------------------------------------------------------
!* Add laser field Ex
!   800nm (0.8 micron) - nearly uniform AC field
!
!  t is in 10^-15 sec (omega in 10^15 unit)
      omega = c1*(2*pi/lambda)    ! 2.36d+15 /sec
      ak    = 2*pi/lambda  ! 7.85d+4/cm <- 0.8 micron
!
      if(it.eq.1) then
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) '  '
        write(11,*) '*Drive of E*H field...'
        write(6,123) omega,ak
  123   format(' omega (1/sec, c*2*pi/lambda)=',1pd14.7,/,  &
               ' ak (1/cm, 2*pi/lambda)=',d14.7)
        close(11)
      end if
      end if
!
!
      E0= 0.d0
      if(.not.ifeqlib) then
!  READ_conf is used
        E0= sqrt(intens/6.d17) *7.0d+7  ! statV/cm
!
!       E0 = sqrt(1.d22/6.d17) *7.d+7 ! 70, 1*10^22 W/cm^2 
!       E0 = sqrt(4.d20/6.d17) *7.d+7 ! 60, 6[4]*10^20 W/cm^2 
!       E0 = sqrt(1.7d19/6.d17)*7.d+7 ! 40, 1.7*10^19 W/cm^2
!       E0 = sqrt(6.d17/6.d17) *7.d+7 ! 20, 5[6]*10^17 W/cm^2
!                                     !   -> 2.1d+12 V/m (eps0: MKSA). 
!                                     !   Then, 1 statV= 300 V, and
!                                     !      one has 7.d+7 statV/cm      
!       E0 = sqrt(3.0d16/6.d17) *7.d+7 ! 10, 3*10^16 W/cm^2 
!
!       lambda = 800.d-7   ! 800nm= 0.8d-4 cm - period 3 fs
!
        tp= (2.6666667d0 *3.d-15)/3.3333333d0
        taup= 5.d-15/3.3333333d0
!
!  Pulse
!       if(t.le.tp) then
!         Eta = E0 * sin(ky*yg -omegat*t) *exp(-((t-tp)/taup)**2)
!
!  Semi-infinite (sine)
!       else
!         Eta = E0 * sin(ky*yg -omega*t)
!       else
!         Eta = 0
!       end if
!
!       do i= 1,nCLp
!       ffr(i,1)= ffr(i,1) +pref_EM*ch(i)*Eta   ! in x-direction
!       ffr(i,3)= ffr(i,3) +pref_EM*ch(i)*Eta   ! now in z-direction
!       end do
      end if
!---------------------------------------------------------
!
!* During initialization, fix CNT C(+Z)
!
      if(ifeqlib) then
        do i= 1,ns          ! C(+Z) and Au of CNT
        ffr(i,1)= 0
        ffr(i,2)= 0
        ffr(i,3)= 0
        end do
!
        do i= ns+1,ns+np    ! H+ of pellet
        ffr(i,1)= 0
        ffr(i,2)= 0
        ffr(i,3)= 0
        end do
!
        do i= ns+np+1,nCLp  ! electron
        ffr(i,1)= 0
        ffr(i,2)= 0
        ffr(i,3)= 0
        end do
      end if
!
!  envelop
      yg0= 3.5d0*0.800d-4
      xz0= 0.1d-4        ! =0.1 micron.is a cutoff distance
      yg02= yg0**2
      xz02= xz0**2
!
      p_xyz= 0.800d-4*(-3.d0 +tg/2.6666667d-15)
!                             ++
!  for initial condition
      if(.not.ifeqlib .and. it.eq.1) then
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        xx= xmin3 +Lx3*(l-1)/mx +hx/2.d0
        yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
        zz= zmin3 +Lz3*(n-1)/mz +hz/2.d0
        ff= exp(-( (yy-p_xyz)**2/yg02 +(xx**2+zz**2)/xz02))
!
        etz(l,m,n)= E0*sin(omega*tg -ak*yy) *ff
        btx(l,m,n)= E0*cos(omega*(tg-dth) -ak*yy) *ff
        end do
        end do
        end do
      end if
!
! Electromagnetic fields
!  B(n+1/2): B(n-1/2), E(n)
!
!      call clocks (cp1,wall1,size)
!
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      xx= xmin3 +Lx3*(l-1)/mx +hx/2.d0
      yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
      zz= zmin3 +Lz3*(n-1)/mz +hz/2.d0
      ff= exp(-( (yy-p_xyz)**2/yg02 +(xx**2+zz**2)/xz02))
!
      btx(l,m,n)= btx(l,m,n) -cdt*((etz(l,m+1,n) -etz(l,m-1,n))/hy2   &
                                  -(ety(l,m,n+1) -ety(l,m,n-1))/hz2)  &
                  -dt*omega*E0*sin(omega*tg -ak*yy)*ff
!
      bty(l,m,n)= bty(l,m,n) -cdt*((etx(l,m,n+1) -etx(l,m,n-1))/hz2   &
                                  -(etz(l+1,m,n) -etz(l-1,m,n))/hx2) 
      btz(l,m,n)= btz(l,m,n) -cdt*((ety(l+1,m,n) -ety(l-1,m,n))/hx2   &
                                  -(etx(l,m+1,n) -etx(l,m-1,n))/hy2)
      end do
      end do
      end do
!
      do n= 1,mz
      do m= 1,my
      do l= 1,mx
      xx= xmin3 +Lx3*(l-1)/mx +hx/2.d0
      yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
      zz= zmin3 +Lz3*(n-1)/mz +hz/2.d0
      ff= exp(-( (yy-p_xyz)**2/yg02 +(xx**2+zz**2)/xz02))
!
      HH(l,m,n)= E0*cos(omega*tg -ak*yy)*ff  ! at tg
      end do
      end do
      end do
!
      ifeh= 2
      call filt3 (btx,bty,btz,EE,HH,ifeh)
!
!       call clocks (cp2,wall2,size)
!
!  Faster by one pass
!  vx(n) <- vx(n+1/2) ??
!
      do n= 1,mz
      do m= 1,my
      do l= 1,mx
      cjx(l,m,n)= 0
      cjy(l,m,n)= 0
      cjz(l,m,n)= 0
      end do
      end do
      end do
!
!     ********************************************************
      do i= 1,nCLp
      m_gamma= sqrt(am(i)**2 +(px(i)**2+py(i)**2+pz(i)**2)/c2) 
      vx(i)= px(i)/m_gamma
      vy(i)= py(i)/m_gamma
      vz(i)= pz(i)/m_gamma
      end do
!     ********************************************************
!
!  (plus,minus) sampling
      do 50 i= 1,nCLp
      l= mx*(xg(i)-xmin3)/Lx3 +1.5000d0
      m= my*(yg(i)-ymin3)/Ly3 +1.5000d0
      n= mz*(zg(i)-zmin3)/Lz3 +1.5000d0
!
      if(l.lt.1 .or. l.ge.mx .or.  &       ! .ge. two-points
         m.lt.1 .or. m.ge.my .or.  &
         n.lt.1 .or. n.ge.mz) go to 50
!
      xx= mx*(xg(i)-xmin3)/Lx3 +1.0001d0
      yy= my*(yg(i)-ymin3)/Ly3 +1.0001d0
      zz= mz*(zg(i)-zmin3)/Lz3 +1.0001d0
      ll= xx  ! ll<xx<ll+1
      mm= yy
      nn= zz
!
      chvx= ch(i)*vx(i)
      chvy= ch(i)*vy(i)
      chvz= ch(i)*vz(i)
!
      cjx(ll+1,m,n)= cjx(ll+1,m,n) +chvx*(xx-ll)   ! vx(n): vx(n+1/2) ??
      cjx(ll,  m,n)= cjx(ll,  m,n) +chvx*(1+ll-xx)
      cjx(l,mm+1,n)= cjx(l,mm+1,n) +chvx*(yy-mm)
      cjx(l,mm  ,n)= cjx(l,mm,  n) +chvx*(1+mm-yy)
      cjx(l,m,nn+1)= cjx(l,m,nn+1) +chvx*(zz-nn)
      cjx(l,m,nn  )= cjx(l,m,nn  ) +chvx*(1+nn-zz) 
!
      cjy(ll+1,m,n)= cjy(ll+1,m,n) +chvy*(xx-ll)
      cjy(ll,  m,n)= cjy(ll,  m,n) +chvy*(1+ll-xx)
      cjy(l,mm+1,n)= cjy(l,mm+1,n) +chvy*(yy-mm)
      cjy(l,mm  ,n)= cjy(l,mm,  n) +chvy*(1+mm-yy)
      cjy(l,m,nn+1)= cjy(l,m,nn+1) +chvy*(zz-nn)
      cjy(l,m,nn  )= cjy(l,m,nn  ) +chvy*(1+nn-zz) 
!
      cjz(ll+1,m,n)= cjz(ll+1,m,n) +chvz*(xx-ll)
      cjz(ll,  m,n)= cjz(ll,  m,n) +chvz*(1+ll-xx)
      cjz(l,mm+1,n)= cjz(l,mm+1,n) +chvz*(yy-mm)
      cjz(l,mm  ,n)= cjz(l,mm,  n) +chvz*(1+mm-yy)
      cjz(l,m,nn+1)= cjz(l,m,nn+1) +chvz*(zz-nn)
      cjz(l,m,nn  )= cjz(l,m,nn  ) +chvz*(1+nn-zz) 
   50 continue 
!
      do n= 1,mz
      do m= 1,my
      do l= 1,mx
      cjx(l,m,n)= cjx(l,m,n)/3.d0  ! just unity for cjx-cjz
      cjy(l,m,n)= cjy(l,m,n)/3.d0
      cjz(l,m,n)= cjz(l,m,n)/3.d0
      end do
      end do
      end do
!
!       call clocks (cp3,wall3,size)
!
!
      if(first_es) then 
        if(n_twice.ge.1) first_es= .false.
!
!  Fortran 
!    Parallel
!       n_thread= 4
!       call dfftw_init_threads (ierror)
!       call dfftw_plan_with_nthreads (n_thread)
!
!         if(ionode) then
!           write(11,*) '*threads with 4 (must be .not. 0) are ',ierror
!         end if
!
!    FFTW3: NEC 2003 style normalized for the forwrd and backward cycle 
!      Generic style of FFTW3 is not normalized (FFTW homepage)
!
        plan= fftw_plan_r2r_3d  &
             (mz,my,mx,qq,qq_c,FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00, &
              FFTW_ESTIMATE) 
!  
        pinv= fftw_plan_r2r_3d  &
             (mz,my,mx,psi_c,psi,FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00, &
              FFTW_ESTIMATE)
      end if 
!
!
      if(mod(it,5).eq.1 .or. read_10) then  ! it=1, or restart
        if(n_twice.ge.1) read_10= .false.   ! n_twice= 1, or restart by 0,1
!
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        qq(l,m,n)= 0
        end do
        end do
        end do
!
        prho= 4*pi/dV
        fnml= 8.d0/(mx*my*mz)  !<- FFTW3 generic
!
        do i= 1,nCLp
        l= mx*(xg(i)-xmin3)/Lx3 +1.5000d0
        m= my*(yg(i)-ymin3)/Ly3 +1.5000d0
        n= mz*(zg(i)-zmin3)/Lz3 +1.5000d0
!
        if(l.lt.1 .or. l.gt.mx .or.       & 
           m.lt.1 .or. m.gt.my .or.       &
           n.lt.1 .or. n.gt.mz) go to 777
!
        qq(l,m,n)= qq(l,m,n) +prho*ch(i)
  777   end do
!
        call fftw_execute_r2r (plan,qq,qq_c)
!
        do n= 1,mz 
        do m= 1,my
        do l= 1,mx
        akx= pi*l/Lx3    ! RODFT00, calculated for l,m,n >= 1
        aky= pi*m/Ly3
        akz= pi*n/Lz3
!
        gam1= 10.d0  ! 7.d0 
        gam2= 20.d0  !14.d0 
        gax= gam1*(l-1)/mx
        gay= gam1*(m-1)/my
        gaz= gam2*(n-1)/mz
        gamma= exp(-(gax**2 +gay**2 +gaz**2))
!
        psi_c(l,m,n)= gamma * qq_c(l,m,n)/(akx**2 +aky**2 +akz**2)
!       psi_c(l,m,n)= fnml*gamma * qq_c(l,m,n)/(akx**2 +aky**2 +akz**2)
        end do                               ! FFTW3 generic
        end do
        end do
!
        call fftw_execute_r2r (pinv,psi_c,psi) 
      end if
!
          if(iwrt1.eq.0 .and. ionode) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
          n= 1
          m= 1
          do l= 1,60,3
          write(11,990) l,m,n,psi_c(l,m,n)
  990     format('l,m,n,psi_c=',3i5,2x,1pd11.4)
          end do
!
          close(11)
          end if
!
!  cjx(n+1/2): vx(n+1/2)
!   Separate the transverse component from the longitudinal one: 
!   J_trans= J -(J*EL)*EL/|ELl^2
!
      do n= 2,mz-1 
      do m= 2,my-1
      do l= 2,mx-1
      pp1= -(psi(l+1,m,n) -psi(l-1,m,n))/hx2
      pp2= -(psi(l,m+1,n) -psi(l,m-1,n))/hy2
      pp3= -(psi(l,m,n+1) -psi(l,m,n-1))/hz2
!
!  Cross product of current and ES field
      ccj= (cjx(l,m,n)*pp1 +cjy(l,m,n)*pp2 +cjz(l,m,n)*pp3)  &
                                    /(pp1**2 +pp2**2 +pp3**2)
!
      cjx(l,m,n)= cjx(l,m,n) -ccj*pp1
      cjy(l,m,n)= cjy(l,m,n) -ccj*pp2
      cjz(l,m,n)= cjz(l,m,n) -ccj*pp3
      end do
      end do
      end do
!
!
      if(iwrt3.eq.0 .and. ionode) then
      if(.not.ifeqlib) then
         if(ionode) then
         OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
               status='unknown',position='append',form='formatted')
!
         write(11,*) 'write(30)... t=',t
         close(11)
         end if
!
         OPEN (unit=30,file=praefixc//'.30'//suffix2,              &
              status='unknown',position='append',form='unformatted')
!
         do n= 1,mzh
         do m= 1,myh
         do l= 1,mxh
         qq1(l,m,n)= qq(l,m,n)
         psi1(l,m,n)= psi(l,m,n)
         end do
         end do
         end do
!
         write(30) t,qq1,psi1  ! real*4
         close(30)
      end if
      end if
!
!       call clocks (cp4,wall4,size)
!
! Save for etx
!
      do n= 1,mz
      do m= 1,my
      do l= 1,mx
      etx1(l,m,n)= etx(l,m,n)
      ety1(l,m,n)= ety(l,m,n)
      etz1(l,m,n)= etz(l,m,n)
      end do
      end do
      end do
!
!  ET(n+1): E(n), B(n+1/2)
!   p4dt= 4*pi*5.d-19
!
      p_xyz= 0.800d-4*(-3.d0 +(tg+dth)/2.6666667d-15)
!                              ******
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      xx= xmin3 +Lx3*(l-1)/mx +hx/2.d0
      yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
      zz= zmin3 +Lz3*(n-1)/mz +hz/2.d0
      ff= exp(-( (yy-p_xyz)**2/yg02 +(xx**2+zz**2)/xz02))
!
      etx(l,m,n)= etx(l,m,n)                               &
                  +cdt*((btz(l,m+1,n) -btz(l,m-1,n))/hy2   &
                       -(bty(l,m,n+1) -bty(l,m,n-1))/hz2)  &
                  -p4dt*cjx(l,m,n)/dV
      ety(l,m,n)= ety(l,m,n)                               &
                  +cdt*((btx(l,m,n+1) -btx(l,m,n-1))/hz2   &
                       -(btz(l+1,m,n) -btz(l-1,m,n))/hx2)  &
                  -p4dt*cjy(l,m,n)/dV
!
      etz(l,m,n)= etz(l,m,n)                               &
                  +cdt*((bty(l+1,m,n) -bty(l-1,m,n))/hx2   &
                       -(btx(l,m+1,n) -btx(l,m-1,n))/hy2)  &
                  -p4dt*cjz(l,m,n)/dV                      &
                  +dt*omega*E0*cos(omega*(tg+dth) -ak*yy)*ff
      end do                      !  at tg+dth
      end do
      end do
!
      do n= 1,mz
      do m= 1,my
      do l= 1,mx
      xx= xmin3 +Lx3*(l-1)/mx +hx/2.d0
      yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
      zz= zmin3 +Lz3*(n-1)/mz +hz/2.d0
      ff= exp(-( (yy-p_xyz)**2/yg02 +(xx**2+zz**2)/xz02))
!
      EE(l,m,n)= E0*sin(omega*(tg+dth) -ak*yy)*ff
      end do
      end do
      end do
!
      ifeh= 1  ! E
      call filt3 (etx,ety,etz,EE,HH,ifeh)
!
!
!  ET(n+1/2): the intermediate time
!
      do n= 1,mz
      do m= 1,my
      do l= 1,mx
      etx1(l,m,n)= (etx(l,m,n) +etx1(l,m,n))/2.d0
      ety1(l,m,n)= (ety(l,m,n) +ety1(l,m,n))/2.d0
      etz1(l,m,n)= (etz(l,m,n) +etz1(l,m,n))/2.d0
      end do
      end do
      end do
!
!       call clocks (cp5,wall5,size)
!
!  EL(n+1/2)
!   All region (infinite)        Syncro by ipar=1,2,...,size ranks
!
!                                    +++ ++++++ ++++++ +++++ 
      call forces (xg,yg,zg,ch,am,ag,ffr,E_C_r1,E_C_r2,E_LJ2,    &
                   Lx3,Ly3,Lz3,rcut_Clf,rcutLJ,prefC_LJ,pref_LJ, &
                   size,ipar,igrp,n_twice,ns,np,nCLp) 
!
      E_C_r = E_C_r2    ! with slow and fast changes
      E_LJ  = E_LJ2
!
!       call clocks (cp6,wall6,size)
!
!  This includes vx(n) ??
!    E(n+1/2), B(n+1/2)
!
      do i= 1,nCLp
!                1/cm      1/sec 
!     ***
      l= mx*(xg(i)-xmin3)/Lx3 +1.5000d0
      m= my*(yg(i)-ymin3)/Ly3 +1.5000d0
      n= mz*(zg(i)-zmin3)/Lz3 +1.5000d0
!
      if(l.lt.1 .or. l.ge.mx .or.  &  ! except L< 1 or L>= mx, two-points
         m.lt.1 .or. m.ge.my .or.  &
         n.lt.1 .or. n.ge.mz) then
!                                          4/01/2017
        ff= exp(-( (yg(i)-p_xyz)**2/yg02 +(xg(i)**2+zg(i)**2)/xz02))
!
        etxi= 0
        etyi= 0
        etzi= E0*sin(omega*(tg+dth) -ak*yg(i))*ff  ! statV/cm
        btxi= E0*cos(omega*(tg+dth) -ak*yg(i))*ff  !  only (etz, btx) extention
        btyi= 0
        btzi= 0
      else
        xx= mx*(xg(i)-xmin3)/Lx3 +1.0001d0
        yy= my*(yg(i)-ymin3)/Ly3 +1.0001d0
        zz= mz*(zg(i)-zmin3)/Lz3 +1.0001d0
        ll= xx  ! ll<xx<ll+1
        mm= yy
        nn= zz
!
!  statV/cm
!    Transverse field
        etxi= ((xx-ll)*etx1(ll+1,m,n) +(1+ll-xx)*etx1(ll,m,n)       & 
              +(yy-mm)*etx1(l,mm+1,n) +(1+mm-yy)*etx1(l,mm,n)       &
              +(zz-nn)*etx1(l,m,nn+1) +(1+nn-zz)*etx1(l,m,nn))/3.d0
        etyi= ((xx-ll)*ety1(ll+1,m,n) +(1+ll-xx)*ety1(ll,m,n)       &
              +(yy-mm)*ety1(l,mm+1,n) +(1+mm-yy)*ety1(l,mm,n)       &
              +(zz-nn)*ety1(l,m,nn+1) +(1+nn-zz)*ety1(l,m,nn))/3.d0
        etzi= ((xx-ll)*etz1(ll+1,m,n) +(1+ll-xx)*etz1(ll,m,n)       &
              +(yy-mm)*etz1(l,mm+1,n) +(1+mm-yy)*etz1(l,mm,n)       &
              +(zz-nn)*etz1(l,m,nn+1) +(1+nn-zz)*etz1(l,m,nn))/3.d0 &
              + E0*sin(omega*(tg+dth) -ak*yg(i))*ff 
!
        btxi= ((xx-ll)*btx(ll+1,m,n) +(1+ll-xx)*btx(ll,m,n)         &
              +(yy-mm)*btx(l,mm+1,n) +(1+mm-yy)*btx(l,mm,n)         &
              +(zz-nn)*btx(l,m,nn+1) +(1+nn-zz)*btx(l,m,nn))/3.d0   &
              + E0*cos(omega*(tg+dth) -ak*yg(i))*ff 
        btyi= ((xx-ll)*bty(ll+1,m,n) +(1+ll-xx)*bty(ll,m,n)         &
              +(yy-mm)*bty(l,mm+1,n) +(1+mm-yy)*bty(l,mm,n)         &
              +(zz-nn)*bty(l,m,nn+1) +(1+nn-zz)*bty(l,m,nn))/3.d0
        btzi= ((xx-ll)*btz(ll+1,m,n) +(1+ll-xx)*btz(ll,m,n)         &
              +(yy-mm)*btz(l,mm+1,n) +(1+mm-yy)*btz(l,mm,n)         &
              +(zz-nn)*btz(l,m,nn+1) +(1+nn-zz)*btz(l,m,nn))/3.d0
      end if
!
!  p(n+1): p(n), ffr(n+1/2), etx(n+1/2), btx(n+1/2), but v(n) ??
!  *fields:  statV/cm
!    ffr: Longitudinal
!    etx: Transverse
!
      px(i)= px(i) +dt*(ffr(i,1)                                   &
                       +ch(i)*(etxi +(vy(i)*btzi-vz(i)*btyi)/c1))
      py(i)= py(i) +dt*(ffr(i,2)                                   &
                       +ch(i)*(etyi +(vz(i)*btxi-vx(i)*btzi)/c1))
      pz(i)= pz(i) +dt*(ffr(i,3)                                   &
                       +ch(i)*(etzi +(vx(i)*btyi-vy(i)*btxi)/c1))
!
!  vx(n+1): from px(n+1)
      m_gamma= sqrt(am(i)**2 +(px(i)**2+py(i)**2+pz(i)**2)/c2)
      vx(i)= px(i)/m_gamma
      vy(i)= py(i)/m_gamma
      vz(i)= pz(i)/m_gamma
!
!  xg(n+3/2): xg(n+1/2)
      xg(i)= xg(i) +dt*vx(i)
      yg(i)= yg(i) +dt*vy(i)
      zg(i)= zg(i) +dt*vz(i)
      end do
!
!
      if(.false.) then
!     if(iwrt2.eq.0 .and. ionode) then
        if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        l= mx/2+1
        n= mz/2+1
        write(11,*) 'tg=',tg
        write(11,*) ' etz,btx; EE,HH (given sin/cos)...'
!
        do m=1,my,20
        write(11,999) m,etz(l,m,n),btx(l,m,n),EE(l,m,n),HH(l,m,n)
  999   format(i4,1p2d13.3,2x,2d13.3)
        end do
!
        close(11)
        end if
      end if
!
!       call clocks (cp7,wall7,size)
!     if(ionode) then
!       write(11,993) cp7-cp1,cp2-cp1,cp3-cp2,cp4-cp3,cp5-cp4, &
!                     cp6-cp5,cp7-cp6
! 993   format('*Ctime ctime=',f7.3,' 1,2,3,4,5,6=',6f7.3)
!     end if
!       
!---------------------
!* Reflection walls
!---------------------
!* Relax electron distributions within the target
!  Exclude ions !
!
      if(ifeqlib) then
        ifeqlib= .false.
!       cccccccccccccccc
!
!  Top and bottom of the cans
        rr1= R_sp  ! 30.d-08 =30 Ang
!a           **** <- READ_CONF
!        
!       if(.false.) then
!       ----------------
        do i= ns+1,ns+np
        rr0= ranff(0.d0)*rr1  ! cm
        ph0= 2*pi*ranff(0.d0)
!
        xg(i)= rr0*sin(ph0)
        yg(i)= rr0*cos(ph0)
        zg(i)= zg(i)
!
        xg(i+np)= xg(i) +0.1d-8*ranff(0.d0) ! elec= H(+): ns+1 -> ns+np+1
        yg(i+np)= yg(i) +0.1d-8*ranff(0.d0)
        zg(i+np)= zg(i)
        end do
!       ----------------
!       end if
!
        if(it.eq.1) then
          if(ionode) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) '    '
          write(11,*) '# Make Coulomb-optimized distributions...'
          write(11,*) ' proton and e (cm)...'
          write(11,330) (i,xg(i),yg(i),zg(i),xg(i+np),yg(i+np), &
                          zg(i+np),i=ns+1,ns+5)
  330     format(' i=',i7,1p3d13.4,2x,3d13.4)
          end if
        end if
!
!
        if_start= .false.  ! make F after equilibration
!       ---------------
        if(eq_phase) then
          eq_phase= .false.  ! equil phase ends here
!
          if(ionode) then
            OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                  status='unknown',position='append',form='formatted')
!
            write(11,*) 'Initial time at restart...'
            write(11,*) ' 1600: ifeqlib= .false.'
            write(11,*) ' 1640: if_start, eq_phase=',if_start,eq_phase
            close(11)
          end if
!
          cr_table= .true.   ! create p-table: restart
          n_twice= 0         ! must be set =0
!
          itab= 0   
          item3= 0
          tg= 0.
           t= tg
!
          is= 0
          it= 0
          iwa= -1
          iwb= -1
          iwc= -1
          iwd= -1
!
          if(ionode) then
            OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                  status='unknown',position='append',form='formatted')
!
            write(11,*) '    '
            write(11,*) '# Use optimized distribution (given T)...'
            close(11)
          end if
!                    ++
          do i= 1,ns+np  ! CNT + H(+)
          px(i)= 0.d0
          py(i)= 0.d0
          pz(i)= 0.d0
          end do
!
          vmax2= 0.7d0*pthe  ! sqrt(Temp_erg), cm/sec
!                    ++
          do i= ns+np+1,nCLp
          px(i)= 0 !am(i)*dgaus2(vmax2) 
          py(i)= 0 !am(i)*dgaus2(vmax2)
          pz(i)= 0 !am(i)*dgaus2(vmax2)
          end do
!
          go to 1000
        end if
      end if
!
!     if(ifrefl.eq.1) then
!       call reflect_box (xg,yg,zg,px,py,pz,ag,nCLp)
!     else if(ifrefl.eq.2) then
!       call reflect_2Dbox (xg,yg,zg,px,py,pz,xmax,size,ipar,nCLp)
!     end if
!
!
!------------------------------
!*  Diagnosis section.
!------------------------------
!
      if(iwrt1.eq.0 .and. ionode) then
      if(.not.ifeqlib) then
          OPEN (unit=13,file=praefixc//'.13'//suffix2,             &
                status='unknown',position='append',form='unformatted')
!
          call nonfdP (xg,yg,zg,x4,y4,z4,npio)
          write(13)  t,x4,y4,z4
          close(13)
      end if
      end if
!
!* dt= 0.5 (every one of 0.5x10^-15) 
!
      if(iwrt3.eq.0 .and. ionode) then   ! 1.d-15 sec 
      if(.not.ifeqlib) then
         if(ionode) then
         OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
               status='unknown',position='append',form='formatted')
!
         write(11,*) 'write(29) t=',t
         close(11)
         end if
!
         OPEN (unit=29,file=praefixc//'.29'//suffix2,             &
               status='unknown',position='append',form='unformatted')
!
!  Field minus EE, HH
!  real*4
         write(29) t
!
         do n= 2,mz,2
         do m= 2,my,2
         do l= 2,mx,2
         ll= l/2 
         mm= m/2 
         nn= n/2 
         xx= xmin3 +Lx3*(l-1)/mx +hx/2.d0
         yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
         zz= zmin3 +Lz3*(n-1)/mz +hz/2.d0
         ff= exp(-( (yy-p_xyz)**2/yg02 +(xx**2+zz**2)/xz02))
!
         eetx(ll,mm,nn)= etx(l,m,n)
         eety(ll,mm,nn)= ety(l,m,n)
         eetz(ll,mm,nn)= etz(l,m,n) -E0*sin(omega*(tg+dth)-ak*yy)*ff
         bbtx(ll,mm,nn)= btx(l,m,n) -E0*cos(omega*(tg+dth)-ak*yy)*ff
         bbty(ll,mm,nn)= bty(l,m,n)
         bbtz(ll,mm,nn)= btz(l,m,n)
         end do
         end do
         end do
! 
         write(29) eetx,eety,eetz,bbtx,bbty,bbtz  ! real*4
!***
         close(29)
      end if
      end if
!---------------------------------------------------------------------
!
! 1. Energy.
!
      if(iwrt1.eq.0) then
      if(ionode) then
!
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
! --------------------------------------- on major nodes --------------
!
        s0= 0.
        s1= 0.
        s2= 0.
!
        do i= 1,ns   ! c2 is omitted
        s0= s0 +(px(i)**2+py(i)**2+pz(i)**2)/(2.d0*am(i))
        end do
!                                           !  subtract rest energy mc^2
        do i= ns+1,ns+np
        s1= s1 +(px(i)**2+py(i)**2+pz(i)**2)/(2.d0*am(i))
        end do
!
        do i= ns+np+1,nCLp   ! all electrons
        s2= s2 +(px(i)**2+py(i)**2+pz(i)**2)/(2.d0*am(i))
        end do
!
        s0= s0/ns    ! C and Au
        s1= s1/np
        s2= s2/(nCLp-(ns+np))
!
!** 10.28.2016
        if(it.eq.1) then
          write(11,*) ' Wkin_E: s0,s1,s2 (per (p^2/2m))'
          write(11,*)  s0,s1,s2
          write(11,*) ' '
        end if
!
!*  em field
        setx= 0
        setz= 0
        sbtx= 0
        sbtz= 0
        psi2= 0
!
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        setx= setx +etx(l,m,n)**2 +ety(l,m,n)**2
        setz= setz +etz(l,m,n)**2
!
        sbtx= sbtx +btx(l,m,n)**2 +bty(l,m,n)**2
        sbtz= sbtz +btz(l,m,n)**2
        end do
        end do
        end do
!
        do n= 2,mz-1
        do m= 2,my-1
        do l= 2,mx-1
        pp1= -(psi(l+1,m,n) -psi(l-1,m,n))/hx2
        pp2= -(psi(l,m+1,n) -psi(l,m-1,n))/hy2
        pp3= -(psi(l,m,n+1) -psi(l,m,n-1))/hz2
        psi2= psi2 +pp1**2 +pp2**2 +pp3**2
        end do
        end do
        end do
!
!             ****** per mx*my*mz
        delV= 1.d0/(8*pi)
        setx= sqrt(delV*setx)/(mx*my*mz)
        setz= sqrt(delV*setz)/(mx*my*mz)
        sbtx= sqrt(delV*sbtx)/(mx*my*mz)
        sbtz= sqrt(delV*sbtz)/(mx*my*mz)
!
!       ***********
        time(is)= t
!
!
        rr= 0
        ani= 0
        rcore= 7*R_sp   ! tentative radius (in micron)
!
        do i= ns+1,ns+np
        r1= sqrt(xg(i)**2 +yg(i)**2 +zg(i)**2)
        rr= amax1(rr, r1)
        if(r1.le.rcore) ani= ani +1.   ! ions within r < rcore
        end do
        Rion(is)= rr  ! radius of ion sphere front
!
        ncoe= 0
        ane = 0
        do i= ns+np+1,ns+np+np
        r1= sqrt(xg(i)**2 +yg(i)**2 +zg(i)**2)
        if(r1.le.Rion(is)) ncoe= ncoe +1   ! electrons inside of ion sphere
        if(r1.le.rcore) ane= ane +1.       ! electrons within r < rcore
        end do
!          **    *
        if(it.le.1) then
          Nele0= ncoe
          write(11,*) '  '
          write(11,*) 'Ne(r<Rion)(0) = ',Nele0
        end if
!
        Nele(is)= ncoe/(Nele0+1.e-5)  ! electrons remaining within the initial sphere
!
        Pot0= 0
        nskip= nCLp/2000
        nsk= nskip/2
!
        do i= 1,nCLp,nskip
        do j= nsk,nCLp,nskip
        dx= xg(i) -xg(j)
        dy= yg(i) -yg(j)
        dz= zg(i) -zg(j)
        Pot0= Pot0 + ch(i)*ch(j)  &
                              /amax1(1.e-8,sqrt(dx**2 +dy**2 +dz**2))
        end do              !        1 Ang
        end do
!
        dPot(is)= Pot0/2000.
!
!---------------------------------------------------------------------
!          **    *
        if(it.le.1) then
!                ++
           W_1p= s2   ! normalize by W_el/nq
           write(11,*) 'W_kin,el(0) = ',W_1p
        end if
!
        ekin(is) = s0 ! C(+Z), by /np
        ekn1(is) = s1 ! H(+)
        ekn2(is) = s2 ! All electrons
!
        ecr (is) = E_C_r
        elj (is) = E_LJ
!*
        ppot(is) = 0
        etot(is) = s0 +s1 +s2 +E_C_r +E_LJ      &
                        +setx +setz +sbtx +sbtz
!
        sex(is)= setx
        sez(is)= setz
        sbx(is)= sbtx
        sbz(is)= sbtz
! 
!* Gyration radius
!
        Rgy0= 0.
        do i= 1,ns
        Rgy0= Rgy0 +xg(i)**2 +yg(i)**2 +zg(i)**2
        end do
!
        Rgy1= 0.
        do i= ns+1,ns+np
        Rgy1= Rgy1 +xg(i)**2 +yg(i)**2 +zg(i)**2
        end do
!
        Rgy2= 0.         ! ***  pellet electrons
        do i= ns+np+1,ns+np+np
        Rgy2= Rgy2 +xg(i)**2 +yg(i)**2 +zg(i)**2
        end do
!
        Rgy3= 0.       ! ***  electrons from C(+Z)
        do i= ns+np+np+1,nCLp
        Rgy3= Rgy3 +xg(i)**2 +yg(i)**2 +zg(i)**2
        end do
!  Average
        Rgyc(is)= sqrt(Rgy0/ns)          ! C+Au
        Rgyi(is)= sqrt(Rgy1/(np+1.e-5))  ! H
        Rgye(is)= sqrt(Rgy2/(np+1.e-5))  ! el1
        Rgyf(is)= sqrt(Rgy3/((nZ+nZA)*(ns/2)))  ! el2
!
!
        if(first11) then
          first11= .false.
!
          if(ipar.eq.1) then
            write(11,700)
  700       format(/,' time(sec):  E_kin0/nC  E_kin1/np  E_kin2/ne ',  &
                  ' E_coulb    E_LJ       setx       setz      ',      &
                  ' sbtx       sbtz       E_tot      Rgyi       ',     &
                  'Rgye       Rgy_C     wall(sec)')
          end if
        end if
!
        write(11,710) time(is),ekin(is),ekn1(is),ekn2(is),ecr(is), &
                      elj(is),sex(is),sez(is),sbx(is),sbz(is),etot(is), &
                      Rgyi(is),Rgye(is),Rgyc(is),dcpu
  710   format('t=',1pd9.2,13d11.2,d11.3) 
!
        close(11)
! --------------------------------------- on major nodes --------------
      end if
      end if
!
! 2. Particle plots.
! ---------------------- 
!*  Total potential
! ---------------------- 
      if(.not.ifeqlib) then
!---
      if(iwrt2.eq.0) then
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,*) '*plots are made... at tg=',tg
        close(11)
!
        OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
!
        call ppl3da (xg,yg,zg,ch,ag,Rgyi(is),Rgye(is),ns,np,nq,  &
                     nCLp,igrp) 
        call ppl3dz (xg,yg,zg,ch,ag,Rgyi(is),Rgye(is),ns,np,nq,  &
                     nCLp,igrp)
!
        call vdistr (xg,yg,zg,vx,vy,vz,ns,np,nCLp) 
        call edistr (vx,vy,vz,am,ns,np,nCLp)
        close (77)
!
        OPEN (unit=23,file=praefixc//'.23'//suffix2,               &
              status='unknown',position='append',form='unformatted')
        write(23) t,xg,yg,zg,vx,vy,vz
        close(23)
      end if
      end if
!---
      end if
!
! --------------------------------------- on major nodes --------------
      GO TO 1000
!
 2000 continue
      if(ionode) then
!
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) ' final: t, tmax (sec)=',tg,tmax
        write(11,*) ' final: cpu1/60., cptot (min)=',dcpu/60.d0,cptot
        close(11)
      end if
!
      return
      end subroutine moldyn
!
!
!-----------------------------------------------------------------------
      subroutine forces (x,y,z,ch,am,ag,ffr,E_C_r1,E_C_r2,E_LJ2,       &
                         Lx3,Ly3,Lz3,rcut_Clf,rcutlj,prefc_lj,pref_lj, &
                         size,ipar,igrp,n_twice,ns,np,nCLp)
!---------------------------------------+++++++-------------------------
      implicit  none
!*
      include    'param_em3p7_Ca.h'
      include    'mpif.h'
!
      integer*4  size,ipar,igrp,ns,np,nCLp,n_twice
!
      real*8     x(npq0),y(npq0),z(npq0),                   &
                 ch(npq0),am(npq0),ag(npq0),ffr(npq0,3),    &
                 ffc(npq0,3),E_C_r1,E_C_r2,E_LJ2,E_LJ,      &
                 Lx3,Ly3,Lz3,alj,dx0,dy0,dz0
!
      real*8          fec,x0,y0,z0
      common/forcepl/ fec(npq0,3)                           ! kokodake
      common/initpos/ x0(npq0),y0(npq0),z0(npq0)            ! moldyn
!
      real*4         time,xp_leng
      common/headr2/ time,xp_leng
!  
      integer*4      i,j,k,l,kk,kc,ll,ncel,lcel,nipl,lipl,liplc,kmax, &
                     ibind,istop,ibox,neigh,ierror,iwrt1,iwrt2,iwrt3, &
                     iwrt4,ix,iy,iz,ix1,iy1,iz1,istp,npar
      common/srflist/ ncel(nc3),lcel(nbxc,nc3)                       ! karita
      common/srflis2/ nipl(n00),lipl(nbxs,n00),liplc(nbx2,n10),kmax  ! kokodake
      common/boxind/ ibind(27,nc3)
      common/abterm/ istop
      common/iotim/  iwrt1,iwrt2,iwrt3,iwrt4
!
      logical        cr_table
      integer*4      itab,nipl0,lipl0,itabs
      common/srflst0/ cr_table,itab,nipl0(n00),lipl0(nbxs,n00)   ! nipl0 wa moldyn
      common/plupdat/ itabs
!
!* afre table
!     integer*4       naf
!     parameter       (naf=(ns0/npq0)*0.3)  ! param_em3p7_Ca.h
      integer*4       nafl,iafl,jafl
      common/afretap/ nafl,iafl(nbxc*naf),jafl(nbxc*naf)
!
      real*8         rcut_Clf,rcutlj,prefC_LJ,pref_LJ,        &
                     dx,dy,dz,tt,r2,r,forceV,ccel,            &
                     driwu,driwu2,addpot,slj,                 &
                     rsi,snt,rcl,rsccl,rsclj,preflj,          &
                     rcut,rcut2,unif1(3),unif2(3)
!-----------
      integer*4     it,is,wrt2
      common/parm1/ it,is
!
      logical    first,if_pc
      data       first/.true./
!*---------------------------------------------------------------
!
      driwu2 = 1.25992104989487316476721060728d0  ! 2**(1/3)
      driwu  = dsqrt(driwu2)                      ! 2**(1/6)
!
      rcut2  = rcut_Clf**2         ! cutoff 
!
      addpot = 1.d0/driwu2**6 -1.d0/driwu2**3
      rsccl = 1.0d-8               ! 1 ang= 10^{-8} cm
!
!*-------------------------------------------
!* update particle table infrequentLy3
!*-------------------------------------------
!
      itab= itab +1
      if(mod(itab,itabs).eq.1 .or. cr_table) then
!
!   Keep it twice for restart n_twice= 0 and 1, 
!   or once for n_twice= 1
!
        if(n_twice.eq.0) cr_table= .true.
        if(n_twice.ge.1) cr_table= .false.     
!
! ++++++ initial/every itabs/restart ++++++++++++
!
!* Step 1: find particles in local boxes
!
        Lx3= xmax3 -xmin3
        Ly3= ymax3 -ymin3
        Lz3= zmax3 -zmin3
!    -----------------
        call Labels
!    -----------------
!
        do k= 1,nc3
        ncel(k)= 0             !  clear the cell register.
        end do
!
!  Step 1: LR-part for 5 steps
!* all particles
!
        do i = 1,nCLp
        ix= isizeX*(x(i)-xmin3)/Lx3 +1.0001
        iy= isizeY*(y(i)-ymin3)/Ly3 +1.0001
        iz= isizeZ*(z(i)-zmin3)/Lz3 +1.0001
!
        if(ix.le.1 .or. ix.ge.isizeX) go to 100
        if(iy.le.1 .or. iy.ge.isizeY) go to 100
        if(iz.le.1 .or. iz.ge.isizeZ) go to 100
!
        ibox = ix + isizeX*(iy-1 + isizeY*(iz-1))
!
        ncel(ibox)= ncel(ibox) +1     ! members of ncel(ibox) of time t
        lcel(ncel(ibox),ibox)= i      !  ibox= 1,...,nc3 cells
  100   end do
!
!* exit
        istop= 0
!
        do k= 1,nc3
        if(ncel(k).gt.nbxs) then
          istop= 1
          go to 140
        end if
        end do
  140   continue
!
        call mpi_allreduce (istop,istp,1,mpi_integer,mpi_sum,  &
                            mpi_comm_world,ierror)
        if(istp.ge.1) then
          istop= istp
          if(ionode) then
            OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                  status='unknown',position='append',form='formatted')
!
            call lp_out (z,nCLp,igrp,istop)
            close(11)
          end if
!
          return
        end if
!
!*  preparation for short-range forces. 
!    register particles around i-th particle
!
        kk= 0
        kc= 0
        nafl= 0
!
        do 200 i= ipar,nCLp,size
        kk= kk +1
!
!* copy initial table components as first entries
!
        if_pc= (i.le.ns)
        if(if_pc) kc= kc +1
!
        nipl(kk)= nipl0(kk)            ! nipl0 in moldyn
!       ++++++++++++++
        if(if_pc) then
          do ll= 1,nipl0(kk)
          liplc(ll,kc)= lipl0(ll,kk)
          end do
        else
          do ll= 1,nipl0(kk)
          lipl(ll,kk)= lipl0(ll,kk)
          end do
        end if
!*
        ix= isizeX*(x(i)-xmin3)/Lx3 +1.0001
        iy= isizeY*(y(i)-ymin3)/Ly3 +1.0001
        iz= isizeZ*(z(i)-zmin3)/Lz3 +1.0001
!
        if(ix.le.1 .or. ix.ge.isizeX) go to 200
        if(iy.le.1 .or. iy.ge.isizeY) go to 200
        if(iz.le.1 .or. iz.ge.isizeZ) go to 200
!
        ibox = ix + isizeX*(iy-1 + isizeY*(iz-1))
!
        do 230 k = 1,27
        neigh= ibind(k,ibox)
        if(neigh.eq.0) go to 230     ! far away
!
        do 240 l= 1,ncel(neigh)      ! find ions in the boxes around i-th
        j= lcel(l,neigh)             !  j-th belongs to this box.
!
        if(j.le.i) go to 240
!          ++++++
        dx= x(i) -x(j)  ! no folding
        dy= y(i) -y(j)
        dz= z(i) -z(j)
!***
        r2 = dx**2 +dy**2 +dz**2
        if(r2.lt.rcut2) then
!
          do ll= 1,nipl0(kk)          ! skip if j are members of t=0
          if(j.eq.lipl0(ll,kk)) go to 240  
          end do
!
          if(if_pc) then
            if(nipl(kk).ge.nbx2) then
              nafl= nafl +1
              iafl(nafl)= i
              jafl(nafl)= j
            else
              nipl(kk)= nipl(kk) +1   ! nipl(kk)=1 and ...
              liplc(nipl(kk),kc)= j
            end if
!
          else
            if(nipl(kk).ge.nbxs) then
              nafl= nafl +1
              iafl(nafl)= i
              jafl(nafl)= j
            else
              nipl(kk)= nipl(kk) +1   ! nipl(kk)=1 and ...
              lipl(nipl(kk),kk)= j
            end if
          end if
!
        end if
  240   continue
  230   continue
  200   continue
!
!* LR-forces: all pairs for 'i-th' particle
!  Keep for 5 time steps
!
        E_C_r1= 0
!
        do i= 1,npq0
        fec(i,1)= 0
        fec(i,2)= 0
        fec(i,3)= 0
        end do
!
!  E_C_r1, fec: reduction, each part for PARALLEL and REDUCTION
!    differences are on ipar value 
!$OMP PARALLEL DEFAULT(NONE)                         &
!$OMP SHARED(ipar,nCLp,size,x,y,z,Lx3,Ly3,Lz3,       &
!$OMP        rcut_Clf,rscCL,ch)                      &
!$OMP PRIVATE(i,j,ix,iy,iz,ix1,iy1,iz1,dx,dy,dz,     &
!$OMP        r,rCL,tt)                               &
!$OMP REDUCTION(+:fec,E_C_r1)
!$OMP DO SCHEDULE(STATIC,1)
!                  +++
        do i= ipar,nCLp,size 
        ix1= isizeX*(x(i)-xmin3)/Lx3 +1.0001
        iy1= isizeY*(y(i)-ymin3)/Ly3 +1.0001
        iz1= isizeZ*(z(i)-zmin3)/Lz3 +1.0001
!
        do j= i+1,nCLp 
        dx = x(i) -x(j)
        dy = y(i) -y(j)
        dz = z(i) -z(j)
!
        ix= isizeX*(x(j)-xmin3)/Lx3 +1.0001
        iy= isizeY*(y(j)-ymin3)/Ly3 +1.0001
        iz= isizeZ*(z(j)-zmin3)/Lz3 +1.0001
!
        r = sqrt(dx**2 + dy**2 + dz**2)
        if(r.gt.rcut_Clf                     &
          .or. ix.le.1 .or. ix.ge.isizeX     &
          .or. iy.le.1 .or. iy.ge.isizeY     &
          .or. iz.le.1 .or. iz.ge.isizeZ     &
          .or. ix1.le.1 .or. ix1.ge.isizeX   &
          .or. iy1.le.1 .or. iy1.ge.isizeY   &
          .or. iz1.le.1 .or. iz1.ge.isizeZ   &
          ) then
!             ++ +++++
          rCL= dmax1(r,rscCL)
          tt = ch(i)*ch(j)/(rCL**2 *r)
!
          fec(i,1) = fec(i,1) + tt*dx
          fec(i,2) = fec(i,2) + tt*dy
          fec(i,3) = fec(i,3) + tt*dz
!
          fec(j,1) = fec(j,1) - tt*dx
          fec(j,2) = fec(j,2) - tt*dy
          fec(j,3) = fec(j,3) - tt*dz
!
          E_C_r1 = E_C_r1 + ch(i)*ch(j)/rCL
        end if
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL
!
!* exit
        istop= 0
        if(nafl.gt.(nbxc*naf)) istop= 1
!
        call mpi_allreduce (istop,istp,1,mpi_integer,mpi_sum,  &
                            mpi_comm_world,ierror)
        if(istp.ge.1) then
          istop= istp
          return
        end if
! ++++++ initial/every itabs/restart ++++++++++++
      end if
!
!
!*****************************************
!*  Short-range Coulomb and LJ forces    *
!*****************************************
!
      kk= 0
      kc= 0
!
      E_C_r2 = E_C_r1      ! in every time step t
      E_LJ2  = 0
!
      do i= 1,npq0 
      ffr(i,1)= fec(i,1)   ! slow changes are here
      ffr(i,2)= fec(i,2)
      ffr(i,3)= fec(i,3) 
      end do
!
! (a) all particles registerd in lipl() and liplc()
!
      do 400 i= ipar,nCLp,size
      kk= kk +1
!
      if_pc= (i.le.ns)
      if(if_pc) kc= kc +1
!
      do 600 k= 1,nipl(kk)
      if(if_pc) then
        j= liplc(k,kc)     ! kc=1,..,ns are members
      else
        j= lipl(k,kk)      ! kk>= ns+1 are counted 
      end if
!
      dx = x(i) -x(j)
      dy = y(i) -y(j)
      dz = z(i) -z(j)
!
      r2 = dx**2 + dy**2 + dz**2
      r  = sqrt(r2)
      rcl= dmax1(r,rsccl)
!
      forceV = ch(i)*ch(j)/(rcl**2 *r)
      E_C_r2 = E_C_r2 + ch(i)*ch(j)/rcl
!
!* short-range forces
!
      if(i.le.ns .and. j.le.ns) then   ! strictly, C or Au pair
        rcut= 3.5d-8
!
        dx0= x0(i) -x0(j)      ! crystal distance = alj
        dy0= y0(i) -y0(j)      !  frozen (ugoka nai)
        dz0= z0(i) -z0(j)
        alj= sqrt(dx0**2 +dy0**2 +dz0**2)
!
        if(alj.lt.rcut) then   ! as r < 3.5 Ang
          preflj= prefc_lj     !   only 1st neighbor
          slj = r/(alj/driwu)  !   note: minimum at alj/driwu
!                  ***
          rsclj= 0.81d0*alj/driwu
        else                   ! larger than equib. radius
          rcut= rcutlj         !   exclusion core (sterlic)
          preflj= pref_lj
!
          slj = r/((ag(i)+ag(j))/driwu)
          rsclj= 0.81d0*(ag(i)+ag(j))/driwu
        end if
!
      else                     ! 3.5 Ang, other than C and Au
        rcut= rcutlj           !   exclusion core (sterlic)
        preflj= pref_lj
!
        slj = r/((ag(i)+ag(j))/driwu)
        rsclj= 0.81d0*(ag(i)+ag(j))/driwu
      end if
!
      ccel= 0
      E_LJ= 0
!        *    ++++
      if(r.le.rcut) then
        rsi = 1.d0/dmax1(slj,0.81d0)**2
        snt = rsi*rsi*rsi
        rsclj= dmax1(rsclj,r)
!                                               ++
        ccel = preflj* snt*(snt-0.5d0)/(rsclj *r)  ! c,au or h,el
        E_LJ = (preflj/12.d0)*(snt*(snt -1.d0) -addpot)
      end if
! 
      ffr(i,1) = ffr(i,1) + (forceV +ccel)*dx
      ffr(i,2) = ffr(i,2) + (forceV +ccel)*dy
      ffr(i,3) = ffr(i,3) + (forceV +ccel)*dz
! 
      ffr(j,1) = ffr(j,1) - (forceV +ccel)*dx
      ffr(j,2) = ffr(j,2) - (forceV +ccel)*dy
      ffr(j,3) = ffr(j,3) - (forceV +ccel)*dz
!
      E_LJ2 = E_LJ2 + E_LJ
  600 continue
  400 continue
!
! (b) "afreta" particles - the same procedures as loop 400
!
      do 700 k= 1,nafl
      i= iafl(k)
      j= jafl(k)
!
      dx = x(i) -x(j)
      dy = y(i) -y(j)
      dz = z(i) -z(j)
!
      r2 = dx**2 + dy**2 + dz**2
      r  = sqrt(r2)
      rcl= dmax1(r,rsccl)
!
      forceV = ch(i)*ch(j)/(rcl**2 *r)
      E_C_r2 = E_C_r2 + ch(i)*ch(j)/rcl
!
      if(i.le.ns .and. j.le.ns) then ! c or au
        rcut= 3.5d-8         ! 1.4ang for cnt 1st neighbor 
!
        dx0= x0(i) -x0(j)    ! crystal distance
        dy0= y0(i) -y0(j)
        dz0= z0(i) -z0(j)
        alj= sqrt(dx0**2 +dy0**2 +dz0**2)
!
        if(alj.lt.rcut) then  ! < 3.5 ang
          preflj= prefc_lj     ! onLy3 1st neighbor
          slj = r/(alj/driwu)  ! note: minimum at alj/driwu
!                  ***
          rsclj= 0.81d0*alj/driwu
        else                  ! larger than equib. radius
          rcut= rcutlj         ! exclusion core (sterlic)
          preflj= pref_lj
!
          slj = r/((ag(i)+ag(j))/driwu)
          rsclj= 0.81d0*(ag(i)+ag(j))/driwu
        end if
!
      else                   ! 3.5 ang, other than c and au
        rcut= rcutlj          ! exclusion core (sterlic)
        preflj= pref_lj
!
        slj = r/((ag(i)+ag(j))/driwu)
        rsclj= 0.81d0*(ag(i)+ag(j))/driwu
      end if
!
      ccel= 0
      E_LJ= 0
!        *    ++++
      if(r.le.rcut) then
        rsi = 1.d0/dmax1(slj,0.81d0)**2
        snt = rsi*rsi*rsi
        rsclj= dmax1(rsclj,r)
!                                               ++
        ccel = preflj* snt*(snt-0.5d0)/(rsclj *r)  ! c,au or h,el
        E_LJ = (preflj/12.d0)*(snt*(snt -1.d0) -addpot)
      end if
! 
      ffr(i,1) = ffr(i,1) + (forceV +ccel)*dx
      ffr(i,2) = ffr(i,2) + (forceV +ccel)*dy
      ffr(i,3) = ffr(i,3) + (forceV +ccel)*dz
! 
      ffr(j,1) = ffr(j,1) - (forceV +ccel)*dx
      ffr(j,2) = ffr(j,2) - (forceV +ccel)*dy
      ffr(j,3) = ffr(j,3) - (forceV +ccel)*dz
!
      E_LJ2 = E_LJ2 + E_LJ
  700 continue
!
!
!  "allreduce" is required for synchronization at a time
!
      npar= npq0*3
      call mpi_allreduce (ffr,ffc,npar,mpi_real8,mpi_sum,  &
                          mpi_comm_world,ierror)
!
      do i= 1,npq0
      ffr(i,1)= ffc(i,1)   ! E_C_r1 +E_C_r2
      ffr(i,2)= ffc(i,2)   ! 
      ffr(i,3)= ffc(i,3)
      end do
!
      unif1(1)= E_C_r2
      unif1(2)= E_LJ2
      call mpi_allreduce (unif1,unif2,2,mpi_real8,mpi_sum,  &
                          mpi_comm_world,ierror)
!
      E_C_r2 = unif2(1)   
      E_LJ2  = unif2(2) 
!
!       if(ionode) then
!       wrt2= 50
!       write(112,990) it,E_C_r1,E_C_r2,E_C_r1+E_C_r2,E_LJ2
! 990   format('it=',i4,' r1, r2, sum, LJ2=',1p3d11.3,2x,d11.3)
!       end if
!
      return
      end subroutine forces
!
!
!----------------------------------------------------------
      subroutine lp_out (z,nCLp,igrp,istop)
!----------------------------------------------------------
      implicit none
      include  'param_em3p7_Ca.h'
!
      real*8     z(npq0)
      integer*4  nCLp,igrp,istop,i,l,k,nskip,ILN,  &
                 ncel,lcel,nipl,lipl,liplc,kmax
      common/srflist/ ncel(nc3),lcel(nbxc,nc3)
      common/srflis2/ nipl(n00),lipl(nbxs,n00),liplc(nbx2,n10),kmax
      real*4     xaxis(2001),yaxis(2001),emax,emin
!
      OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
            status='unknown',position='append',form='formatted')
!
      call symbol (7.0,18.0,0.7,'istop=',0.,6)
      call number (999.,999.,0.7,float(istop),0.,5)
!
      l= 0
      nskip= nc3/2000 +1
!
      do k= 1,nc3,nskip
      l= l +1
      xaxis(l)= k
      yaxis(l)= ncel(k)
      end do
!
      ILN= 1
      call lplmax (yaxis,emax,emin,l)
      call lplot1 (2,2,l,xaxis,yaxis,emax,0.0,ILN,'Ncell(N)',8,  &
                   '   #N   ',8,'        ',8)
!
      l= 0
      nskip= nCLp/2000 +1
!
      do i= 1,nCLp,nskip
      l= l +1
      xaxis(l)= i
      yaxis(l)= z(i)
      end do
!
      call lplmax (yaxis,emax,emin,l)
      call lplot1 (2,3,l,xaxis,yaxis,emax,emin,ILN,' Z(i)   ',8,  &
                   '  #i    ',8,'        ',8)
!
      l= 0
      nskip= kmax/2000 +1
!
      do i= 1,kmax,nskip
      l= l +1
      xaxis(l)= i
      yaxis(l)= nipl(i)
      end do
!
      call lplmax (yaxis,emax,emin,l)
      call lplot1 (3,3,l,xaxis,yaxis,emax,0.0,ILN,'Nipl(i) ',8,  &
                   '  #i    ',8,'        ',8)
!
      call chart
      close (77)
!
      return
      end subroutine lp_out
!
!
!----------------------------------------------------------
      subroutine Labels
!----------------------------------------------------------
!*  Indices of the neighboring (27) sub-boxes.
!   Non-periodic: no folding
!
      implicit none
      include 'param_em3p7_Ca.h'
!
      integer*4   isize2,isize4,icmax,              &
                  i,j,k,n,IPx,IMx,JP,JM,KP,KM,ibind
      common/boxind/ ibind(27,nc3)
!
      isize2 =  isizeX*isizeY
      isize4 =  isize2 + isizeX
!
      icmax= nc3 +1
!
      do i = 1,isizeX
      IPx = i + 1
      IMx = i - 1
      if (i.eq.1) IMx = icmax        ! must be excluded
      if (i.eq.isizeX) IPx = icmax   ! - otherwise, double counted
!
      do j = 1,isizeY
      JP = j + 1
      JM = j - 1
      if (j.eq.1) JM = icmax
      if (j.eq.isizeY) JP = icmax
!
      do k = 1,isizeZ
      KP = k + 1
      KM = k - 1
      if (k.eq.1) KM = icmax
      if (k.eq.isizeZ) KP = icmax
!
      n = i + isizeX*(j-1 + isizeY*(k-1))
      ibind( 1,n) = i   +  j*isizeX +  k*isize2 - isize4
!
      ibind( 2,n) = IPx +  j*isizeX +  k*isize2 - isize4
      ibind( 3,n) = IPx + JM*isizeX +  k*isize2 - isize4 
      ibind( 4,n) = IPx + JP*isizeX +  k*isize2 - isize4
      ibind( 5,n) = i   + JP*isizeX +  k*isize2 - isize4
      ibind( 6,n) = i   +  j*isizeX + KP*isize2 - isize4
      ibind( 7,n) = IMx +  j*isizeX + KP*isize2 - isize4
      ibind( 8,n) = IPx +  j*isizeX + KP*isize2 - isize4
      ibind( 9,n) = i   + JM*isizeX + KP*isize2 - isize4
      ibind(10,n) = IMx + JM*isizeX + KP*isize2 - isize4
      ibind(11,n) = IPx + JM*isizeX + KP*isize2 - isize4
      ibind(12,n) = i   + JP*isizeX + KP*isize2 - isize4
      ibind(13,n) = IMx + JP*isizeX + KP*isize2 - isize4
      ibind(14,n) = IPx + JP*isizeX + KP*isize2 - isize4
!
      ibind(15,n) = i   + JM*isizeX + KM*isize2 - isize4
      ibind(16,n) = IMx + JM*isizeX + KM*isize2 - isize4
      ibind(17,n) = IPx + JM*isizeX + KM*isize2 - isize4
      ibind(18,n) = i   + JP*isizeX + KM*isize2 - isize4 
      ibind(19,n) = IMx + JP*isizeX + KM*isize2 - isize4
      ibind(20,n) = IPx + JP*isizeX + KM*isize2 - isize4
      ibind(21,n) = IPx +  j*isizeX + KM*isize2 - isize4
      ibind(22,n) = IMx +  j*isizeX + KM*isize2 - isize4
      ibind(23,n) = i   +  j*isizeX + KM*isize2 - isize4
      ibind(24,n) = IMx + JP*isizeX +  k*isize2 - isize4
      ibind(25,n) = IMx + JM*isizeX +  k*isize2 - isize4
      ibind(26,n) = IMx +  j*isizeX +  k*isize2 - isize4
      ibind(27,n) = i   + JM*isizeX +  k*isize2 - isize4
      end do
      end do
      end do
!
!* Nullify if ibind() is out of bounds
!
      do n= 1,nc3
      do k= 1,27
      if(ibind(k,n).gt.nc3) ibind(k,n)= 0
      end do
      end do
!
      return
      end subroutine Labels
!
!
!------------------------------------------------------------
      subroutine filt3 (srx,sry,srz,EE,HH,ifeh)
!------------------------------------------------------------
      implicit none
!
      include   'param_em3p7_Ca.h'
      integer*4 l,m,n,kk,ifeh
!
      real*8  srx(mx,my,mz),sry(mx,my,mz),srz(mx,my,mz),  &
              ssx(mx,my,mz),ssy(mx,my,mz),ssz(mx,my,mz),  &
              del,c0,c1,c2,c3,tot,                        &
              EE(mx,my,mz),HH(mx,my,mz)
!
      del= 1/2.d0
      c0= 1.d0                  ! *1 =1
      c1= 1.d0/del              ! *6 =6/d=3.
      c2= 1.d0/(sqrt(2.d0)*del) ! *12=12/sq(2)*d=4.24
      c3= 1.d0/(sqrt(3.d0)*del) ! *8 = 8/sq(3)*d=2.31  total= 10.5
      tot= 1/(1.d0+6/del+12/(sqrt(2.d0)*del)+8/(sqrt(3.d0)*del))
!
! Boundary at 2D cross sections
!   Store the field, although it is NOT small !
!     keep srx -> btx, or srz -> etz
!
      if(ifeh.eq.1) then  ! E
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        srz(l,m,n)= srz(l,m,n) -EE(l,m,n)
        end do
        end do
        end do
!
      else if(ifeh.eq.2) then  ! H
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        srx(l,m,n)= srx(l,m,n) -HH(l,m,n)
        end do
        end do
        end do
      end if
!
      do n= 1,mz
      do l= 1,mx
      srx(l, 1,n)= 0
      srx(l,my,n)= 0
      sry(l, 1,n)= 0
      sry(l,my,n)= 0
      srz(l, 1,n)= 0
      srz(l,my,n)= 0
      end do
      end do
!
      do m= 1,my
      do l= 1,mx
      srx(l,m, 1)= 0
      srx(l,m,mz)= 0 
      sry(l,m, 1)= 0
      sry(l,m,mz)= 0
      srz(l,m, 1)= 0
      srz(l,m,mz)= 0
      end do
      end do
!
      do n= 1,mz
      do m= 1,my
      srx( 1,m,n)= 0
      srx(mx,m,n)= 0
      sry( 1,m,n)= 0
      sry(mx,m,n)= 0
      srz( 1,m,n)= 0
      srz(mx,m,n)= 0
      end do
      end do
!
!      
      kk= 0
  100 kk= kk +1
      if(kk.gt.1) go to 700
!             ***
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      ssx(l,m,n)=                                                     & 
          c3*srx(l+1,m+1,n+1) +c2*srx(l+1,m,n+1) +c3*srx(l+1,m-1,n+1) &
         +c2*srx(l+1,m+1,n  ) +c1*srx(l+1,m,n  ) +c2*srx(l+1,m-1,n  ) & 
         +c3*srx(l+1,m+1,n-1) +c2*srx(l+1,m,n-1) +c3*srx(l+1,m-1,n-1) &
         +c2*srx(l  ,m+1,n+1) +c1*srx(l  ,m,n+1) +c2*srx(l  ,m-1,n+1) &
         +c1*srx(l  ,m+1,n  ) +c0*srx(l  ,m,n  ) +c1*srx(l  ,m-1,n  ) &
         +c2*srx(l  ,m+1,n-1) +c1*srx(l  ,m,n-1) +c2*srx(l  ,m-1,n-1) &
         +c3*srx(l-1,m+1,n+1) +c2*srx(l-1,m,n+1) +c3*srx(l-1,m-1,n+1) &
         +c2*srx(l-1,m+1,n  ) +c1*srx(l-1,m,n  ) +c2*srx(l-1,m-1,n  ) &
         +c3*srx(l-1,m+1,n-1) +c2*srx(l-1,m,n-1) +c3*srx(l-1,m-1,n-1) 
      ssy(l,m,n)=                                                     & 
          c3*sry(l+1,m+1,n+1) +c2*sry(l+1,m,n+1) +c3*sry(l+1,m-1,n+1) &
         +c2*sry(l+1,m+1,n  ) +c1*sry(l+1,m,n  ) +c2*sry(l+1,m-1,n  ) &
         +c3*sry(l+1,m+1,n-1) +c2*sry(l+1,m,n-1) +c3*sry(l+1,m-1,n-1) &
         +c2*sry(l  ,m+1,n+1) +c1*sry(l  ,m,n+1) +c2*sry(l  ,m-1,n+1) &
         +c1*sry(l  ,m+1,n  ) +c0*sry(l  ,m,n  ) +c1*sry(l  ,m-1,n  ) &
         +c2*sry(l  ,m+1,n-1) +c1*sry(l  ,m,n-1) +c2*sry(l  ,m-1,n-1) &
         +c3*sry(l-1,m+1,n+1) +c2*sry(l-1,m,n+1) +c3*sry(l-1,m-1,n+1) &
         +c2*sry(l-1,m+1,n  ) +c1*sry(l-1,m,n  ) +c2*sry(l-1,m-1,n  ) &
         +c3*sry(l-1,m+1,n-1) +c2*sry(l-1,m,n-1) +c3*sry(l-1,m-1,n-1) 
      ssz(l,m,n)=                                                     &
          c3*srz(l+1,m+1,n+1) +c2*srz(l+1,m,n+1) +c3*srz(l+1,m-1,n+1) &
         +c2*srz(l+1,m+1,n  ) +c1*srz(l+1,m,n  ) +c2*srz(l+1,m-1,n  ) &
         +c3*srz(l+1,m+1,n-1) +c2*srz(l+1,m,n-1) +c3*srz(l+1,m-1,n-1) &
         +c2*srz(l  ,m+1,n+1) +c1*srz(l  ,m,n+1) +c2*srz(l  ,m-1,n+1) &
         +c1*srz(l  ,m+1,n  ) +c0*srz(l  ,m,n  ) +c1*srz(l  ,m-1,n  ) &
         +c2*srz(l  ,m+1,n-1) +c1*srz(l  ,m,n-1) +c2*srz(l  ,m-1,n-1) &
         +c3*srz(l-1,m+1,n+1) +c2*srz(l-1,m,n+1) +c3*srz(l-1,m-1,n+1) &
         +c2*srz(l-1,m+1,n  ) +c1*srz(l-1,m,n  ) +c2*srz(l-1,m-1,n  ) &
         +c3*srz(l-1,m+1,n-1) +c2*srz(l-1,m,n-1) +c3*srz(l-1,m-1,n-1) 
      end do
      end do
      end do
!
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      srx(l,m,n)= tot*ssx(l,m,n)
      sry(l,m,n)= tot*ssy(l,m,n)
      srz(l,m,n)= tot*ssz(l,m,n)
      end do
      end do
      end do
      go to 100
!
  700 continue
!
!  Restoration by getting the old field
      if(ifeh.eq.1) then  ! E
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        srz(l,m,n)= srz(l,m,n) +EE(l,m,n)
        end do
        end do
        end do
!
      else if(ifeh.eq.2) then  ! H
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        srx(l,m,n)= srx(l,m,n) +HH(l,m,n)
        end do
        end do
        end do
      end if
!
      return
      end subroutine filt3
!
!
!------------------------------------------------------------
      subroutine filt1 (srx)
!------------------------------------------------------------
      implicit none
!
      include   'param_em3p7_Ca.h'
      integer*4 l,m,n,kk
!
      real*8  srx(mx,my,mz),ssx(mx,my,mz),  &
              del,c0,c1,c2,c3,tot
!
      del= 1/2.d0
      c0= 1.d0                  ! *1 =1
      c1= 1.d0/del              ! *6 =6/d=3.
      c2= 1.d0/(sqrt(2.d0)*del) ! *12=12/sq(2)*d=4.24
      c3= 1.d0/(sqrt(3.d0)*del) ! *8 = 8/sq(3)*d=2.31  total= 10.5
      tot= 1/(1.d0+6/del+12/(sqrt(2.d0)*del)+8/(sqrt(3.d0)*del))
!
! Boundary at 2D cross sections
!
      do n= 1,mz
      do l= 1,mx
      srx(l, 1,n)= 0
      srx(l,my,n)= 0
      end do
      end do
!
      do m= 1,my
      do l= 1,mx
      srx(l,m, 1)= 0
      srx(l,m,mz)= 0 
      end do
      end do
!
      do n= 1,mz
      do m= 1,my
      srx( 1,m,n)= 0
      srx(mx,m,n)= 0
      end do
      end do
!
!      
      kk= 0
  100 kk= kk +1
      if(kk.gt.1) go to 700
!             ***
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      ssx(l,m,n)=                                                     & 
          c3*srx(l+1,m+1,n+1) +c2*srx(l+1,m,n+1) +c3*srx(l+1,m-1,n+1) &
         +c2*srx(l+1,m+1,n  ) +c1*srx(l+1,m,n  ) +c2*srx(l+1,m-1,n  ) & 
         +c3*srx(l+1,m+1,n-1) +c2*srx(l+1,m,n-1) +c3*srx(l+1,m-1,n-1) &
         +c2*srx(l  ,m+1,n+1) +c1*srx(l  ,m,n+1) +c2*srx(l  ,m-1,n+1) &
         +c1*srx(l  ,m+1,n  ) +c0*srx(l  ,m,n  ) +c1*srx(l  ,m-1,n  ) &
         +c2*srx(l  ,m+1,n-1) +c1*srx(l  ,m,n-1) +c2*srx(l  ,m-1,n-1) &
         +c3*srx(l-1,m+1,n+1) +c2*srx(l-1,m,n+1) +c3*srx(l-1,m-1,n+1) &
         +c2*srx(l-1,m+1,n  ) +c1*srx(l-1,m,n  ) +c2*srx(l-1,m-1,n  ) &
         +c3*srx(l-1,m+1,n-1) +c2*srx(l-1,m,n-1) +c3*srx(l-1,m-1,n-1) 
      end do
      end do
      end do
!
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      srx(l,m,n)= tot*ssx(l,m,n)
      end do
      end do
      end do
      go to 100
!
  700 continue
      return
      end subroutine filt1
!
!
!------------------------------------------------------------------------
      subroutine reflect_cylinder (x,y,z,vx,vy,vz,Z_cn,R_cn1,R_cn2,  &
                                   ipar1,ipar2)
!------------------------------------------------------------------------
! Reflection by the cylinder
      implicit none
!
      include    'param_em3p7_Ca.h'
      integer*4   ipar1,ipar2,i
!
      real*8     x(npq0),y(npq0),z(npq0),    &
                 vx(npq0),vy(npq0),vz(npq0), &
                 Z_cn,R_cn1,R_cn2,rr,vpara
!
      do 100 i= ipar1,ipar2
      rr= sqrt(x(i)**2 +y(i)**2)
      if(rr.lt.R_cn1 .or. rr.gt.R_cn2) then
        vpara= (vx(i)*x(i) +vy(i)*y(i))/rr  
        vx(i)= vx(i) -2*vpara*x(i)/rr
        vy(i)= vy(i) -2*vpara*y(i)/rr
      end if
!
      if(abs(z(i)).gt.Z_cn) then
        vz(i)= -vz(i)
      end if
  100 continue
!
      return
      end subroutine reflect_cylinder
!
!
!------------------------------------------------------------------------
      subroutine reflect_cylinder_h (x,y,z,vx,vy,vz,Z_cn1,Z_cn2,  &
                                     R_cn1,R_cn2,ipar1,ipar2)
!------------------------------------------------------------------------
! Reflection by a hollow cylinder: Z_cn1 < z < Z_cn2
      implicit none
!
      include    'param_em3p7_Ca.h'
      integer*4   ipar1,ipar2,i
!
      real*8     x(npq0),y(npq0),z(npq0),     &
                 vx(npq0),vy(npq0),vz(npq0),  &
                 Z_cn1,Z_cn2,R_cn1,R_cn2,rr,vpara
!
      do 100 i= ipar1,ipar2
      rr= sqrt(x(i)**2 +y(i)**2)
      if(rr.lt.R_cn1 .or. rr.gt.R_cn2) then
        vpara= (vx(i)*x(i) +vy(i)*y(i))/rr  
        vx(i)= vx(i) -2*vpara*x(i)/rr
        vy(i)= vy(i) -2*vpara*y(i)/rr
      end if
!
      if(z(i).gt.Z_cn2 .or. z(i).lt.Z_cn1) then
        vz(i)= -vz(i)
      end if
  100 continue
!
      return
      end subroutine reflect_cylinder_h
!
!
!  READ configuration data.
!--------------------------------------------------------------------
      subroutine READ_CONF (np,ifrefl)  ! data
!-------------------------- +++++ ++++++ ++ -------------------------
      implicit none
!
      include    'param_em3p7_Ca.h'
      integer*4  np,ifrefl
!  
      integer*4  nZ,nZA
      real*8     fchar,fcharA
      common/iroha/ nZ,nZA
      common/charg/ fchar,fcharA
!
      real*8   intens,lambda
      common/electr/ intens,lambda
!  
      character  praefix*6,text1*40
!
      REAL*8   rcut_Clf,rcutlj
      REAL*8   Temp,epsCLJ,epsLJ
      common/cutoffrd/ rcut_Clf,rcutlj
      COMMON/ELSTA/  Temp,epsCLJ,epsLJ
!
      integer*4       iiwrt1,iiwrt2,iiwrt3,iiwrt4
      common/iiwrite/ iiwrt1,iiwrt2,iiwrt3,iiwrt4
!
      integer*4   itabs
      common/plupdat/ itabs
!----------------------------------------------------------------
      real*8        a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
      common/physc/ a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
!
      real*8        Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
      common/hx2l/  Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
!
      real*8        pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
!
      real*4        phi,tht,cptot,rgmax
      common/parm9/ cptot
      common/parm4/ phi,tht,rgmax
!
      real*8         R_sp,D_sp,N_sp,Lambda_D,massi,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
!----------------------------------------------------------------
      OPEN (unit=08,file=praefixs//'_config.START'//suffix0,   &
                                               form='formatted')
!
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'READ_CONF: Parameter read... start'
        close(11)
      end if
!*
      read (8,'(a40,a6)') text1, praefix    ! String der simulationserkennung
      read (8,'(a40,i12)')   text1,ifrefl   ! =0,1,2 reflection at open/close/1d open
      read (8,'(a40,f12.0)') text1,cptot    ! Maximum cpu time for each run (min)
      read (8,'(a40,d20.0)') text1,tmax     ! Zeit zum Abbruch, 1.d-15 sec
      read (8,'(a40,d20.0)') text1,dt       ! Zeitschritt, in 1.d-15 sec
!
      read (8,'(a40,i12)')   text1,iiwrt1   ! Write out interval for IWRT1, 1.d-15 sec 
      read (8,'(a40,i12)')   text1,iiwrt2   ! Write out for IWRT2
      read (8,'(a40,i12)')   text1,iiwrt3   ! Write out for IWRT3
      read (8,'(a40,i12)')   text1,itabs    ! Particle table update interval
!
      read (8,'(a40,i12)')   text1,np       ! Number of protons
!
      read (8,'(a40,f12.0)') text1,fchar    ! carbon: charge state
      read (8,'(a40,f12.0)') text1,fcharA   ! Au: charge state: 20 to 80 
!
      read (8,'(a40,i12)')   text1,nZ       ! large lump as nZ=1 
      read (8,'(a40,i12)')   text1,nZA      ! large lump as nZA=4 
      read (8,'(a40,d20.0)') text1,massi    ! mass ratio: mp/me
!
      read (8,'(a40,d20.0)') text1,intens   ! intensity W/cm2; 6.d17
      read (8,'(a40,d20.0)') text1,lambda   ! wavelength 800 nm= 8.d-7 cm
!  
      read (8,'(a40,d20.0)') text1,R_sp     ! in R_sp, cm
      read (8,'(a40,d20.0)') text1,D_sp
!   --------------------
!     R_sp= R_sp/1.d-4
!   --------------------
!
      read (8,'(a40,d20.0)') text1,rcut_Clf ! Coulomb cutoff in Angstrom
      read (8,'(a40,d20.0)') text1,rcutLJ   ! Lennard-Jones cutoff in Angstgrom
      read (8,'(a40,d20.0)') text1,epsCLJ   ! LJ potential between C in eV
      read (8,'(a40,d20.0)') text1,epsLJ    ! LJ potential otherwise in eV
      read (8,'(a40,d20.0)') text1,Temp     ! Electron temperature in eV
      read (8,'(a40,d20.0)') text1,rgmax
!
! -------------------------------
! Convert from Angstrom to micron
!   xmax,... are now read from include file
!
      rcut_Clf= sconv*rcut_Clf  ! Ang to cm
      rcutLJ  = sconv*rcutLJ    ! Ang to cm
!
      rgmax= sconv*rgmax
! -------------------------------
!
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'READ_CONF: Parameter read... end'
        close(11)
      end if
!
      close(08)
!
      return
      end subroutine READ_CONF
!
!
!-----------------------------------------------------------
      subroutine nonfdP (xg,yg,zg,x4,y4,z4,npio)
!-----------------------------------------------------------
      implicit none
      include  'param_em3p7_Ca.h'
!
      integer*4  npio,i
      real*8     xg(npq0),yg(npq0),zg(npq0)
      real*4     x4(npio),y4(npio),z4(npio)
!
      do 100 i= 1,npio
      x4(i)= xg(i)
      y4(i)= yg(i)
      z4(i)= zg(i)
  100 continue
!
      return
      end subroutine nonfdP
!
!
!-------------------------------------------------------------
      subroutine init (x,y,z,px,py,pz,ch,am,ag,           &
                       ipar,istop,if_start,ns,np,nq,nCLp)
!----------------------****--------------- ++ ++ ++ +++ ------
      implicit none
      include  'param_em3p7_Ca.h'
!
      integer*4   ipar,istop,wrt2,ns,np,nq,nZ,nZA,nCLp
      logical     if_start
      real*8      fchar,fcharA,heavy,heavyA
!
      common/iroha/ nZ,nZA
      common/charg/ fchar,fcharA
      common/hev_el/ heavy,heavyA
!
      real*8     x(npq0),y(npq0),z(npq0),                &
                 px(npq0),py(npq0),pz(npq0),             &
                 ch(npq0),am(npq0),ag(npq0),             &
                 vx(npq0),vy(npq0),vz(npq0),             &
                 pi,tg,dt,dth,                           &
                 prefC_LJ,pref_LJ,pthe,tmax,             &
                 a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest, &
                 m_gamma,s1,ranff
      real*4     phi,tht,rgmax,vmax2,dgaus2
      integer*4  i,j,k
!
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm4/ phi,tht,rgmax
      common/physc/ a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
!
      real*8         R_sp,D_sp,N_sp,Lambda_D,massi,               &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,               &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
!
      real*8        rcut_Clf,rcutlj,pi2,cci
      common/cutoffrd/ rcut_Clf,rcutlj
!
      integer*4 sgn(3,6)
      data sgn /1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1/
!
      integer*4   nh1,nh2
      real*8      cci1,cci2,cci3,cci4,cci5,cci6
!
!197.d0*wt_ion
      pi2= 2.d0*pi
!
!* 1. Give charge states
!
      do i= 1,ns/2            ! carbon (tube)
      ch(i)= fchar*e_unit   ! esu
      am(i)= 12.d0*wt_ion
      ag(i)= rd_CP
      end do
!
      do i= ns/2+1,ns         ! Au (coat)
      ch(i)= fcharA*e_unit  ! esu
      am(i)= 197.d0*wt_ion
      ag(i)= rd_CP
      end do
!
      do i= ns+1,ns+np        ! proton (pellet)
      ch(i)= ch_ion
      am(i)= wt_ion
      ag(i)= rd_HP
      end do
!
      do i= ns+np+1,ns+2*np   ! electrons in pellet
      ch(i)= ch_el
      am(i)= wt_el
      ag(i)= rd_el
      end do
!
! Heavy electrons for CNT & Au-associated electrons
!
      nh1= ns+2*np+nZ*(ns/2)
      nh2= nh1 +nZA*(ns/2)
!
      do i= ns+2*np+1,ns+2*np+nZ*(ns/2)
      ch(i)= -heavy*e_unit 
      am(i)=  heavy*m_unit !wt_el
      ag(i)= rd_el
      end do
!
      do i= ns+2*np+nZ*(ns/2)+1,nCLp
      ch(i)= -heavyA*e_unit 
      am(i)=  heavyA*m_unit !*wt_el
      ag(i)= rd_el
      end do
!
        if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' '
        write(11,*) 'ns/2, fchar, am(1)=',ns/2,fchar,12.d0*wt_ion
        write(11,*) 'ns/2, fchar, am(2)=',ns/2,fcharA,197.d0*wt_ion
        write(11,*) 'np(H),fchar, am(3)=',np,ch_ion/e_unit,wt_ion
        write(11,*) 'np(e),fchar, am(4)=',np,ch_el/e_unit,wt_el
        write(11,*) 'nh1,  fchar, am(5)=',nZ*(ns/2),                &
                      heavy,heavy*m_unit
        write(11,*) 'nh2,  fchar, am(6)=',nCLp-(ns+2*np+nZ*(ns/2)), &
                      heavyA,heavyA*m_unit
        write(11,*) ' '
!
        close(11)
        end if
!
!
      cci1= 0
      do i= 1,ns/2
      cci1= cci1 +ch(i)
      end do
!
      cci2= 0
      do i= ns/2+1,ns
      cci2= cci2 +ch(i)
      end do
!
      cci3= 0
      do i= ns+1,ns+np
      cci3= cci3 +ch(i)
      end do
!
      cci4= 0
      do i= ns+np+1,ns+2*np
      cci4= cci4 +ch(i)
      end do
!
      cci5= 0
      do i= ns+2*np+1,ns+2*np+nZ*(ns/2)
      cci5= cci5 +ch(i)
      end do
!
      cci6= 0
      do i= ns+2*np+nZ*(ns/2)+1,nCLp
      cci6= cci6 +ch(i)
      end do
!
      cci= cci1+cci2+cci3+cci4+cci5+cci6
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,*) ' cci not equal or smaller than 1.d-10...',cci
!
        write(11,*) 'Sum of ch(1,2,3,4,5,6)... '
        write(11,*) ' cci1=',cci1
        write(11,*) ' cci2=',cci2
        write(11,*) ' cci3=',cci3
        write(11,*) ' cci4=',cci4
        write(11,*) ' cci5=',cci5
        write(11,*) ' cci6=',cci6
        write(11,*) '  Total cci (< d-13 esu)=',cci
        write(11,*) '      '
!
        close(11)
      end if
!
      if(abs(cci).gt.1.d-10) then
        call exit
      end if
! -----------------------------
      if(.not.if_start) return
! -----------------------------
!
!***************
!*  Velocity.  *
!***************
!
      vmax2= 0.7d0*pthe
!
      do 530 i= 1,ns+np    ! ions are cold
      px(i)= 0.d0
      py(i)= 0.d0
      pz(i)= 0.d0
  530 continue
!
      s1= 0.
      do 540 i= ns+np+1,nCLp  ! electrons, c is omitted
      px(i)= 0 ! am(i)*dgaus2(vmax2)
      py(i)= 0 ! am(i)*dgaus2(vmax2)
      pz(i)= 0 ! am(i)*dgaus2(vmax2)
!
      s1= s1 +(px(i)**2+py(i)**2+pz(i)**2)/(2.d0*am(i))
  540 continue                      !  subtract rest energy mc^2
!
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,541) s1/nq
  541   format(' Kinetic energy of an electron /mc^2=',1pd13.5)
        close(11)
      end if
!
      return
      end subroutine init
!
!
!---------------------------------------------------------------
      subroutine ggauss
!---------------------------------------------------------------
!  (-3.,3.) case
      implicit none
!
      real*4  fv,vv0,dv,vv,s,sdv,fun
      integer*4 j,ns,k2,k,i
      common/gaus1/ fv(51),vv0,dv
!
      fv(1)=0.0
!
      vv0= -3.
      dv= 2.*abs(vv0)/50.0
!
      vv= vv0
      do 100 j=1,50
      s=0.0
      ns=1000
      k2=ns/2
      sdv=dv/float(ns)
!
      do 130 k=1,k2
      vv=vv +2.0*sdv
      s=s +4.0*fun(vv-sdv) +2.0*fun(vv)
  130 continue
      s= (s +4.0*fun(vv+sdv) +fun(vv+2.0*sdv))*sdv/3.0
      fv(j+1)= fv(j)+s
  100 continue
!
      do 200 i=1,51
      fv(i)=fv(i)/fv(51)
  200 continue
!
      return
      end subroutine ggauss
!
!
!---------------------------------------------------------------
      function dgaus2 (vmax)
!---------------------------------------------------------------
!  vmax dimension
      implicit none
!
      real*4  dgaus2,vmax
      real*4  fv,vv0,dv,s,sdv,fun,eps,x2,y1,y2
      real*8  ranff
      integer*4 j,ns,k2,k
!
      common/gaus1/ fv(51),vv0,dv
!
      eps= ranff(0.d0)
      do 100 k=1,51
      k2=k
      if(fv(k).gt.eps) go to 200
  100 continue
!
  200 y1= fv(k2-1)
      y2= fv(k2)
      x2= (eps-y2)/(y2-y1)+k2
      dgaus2= vmax*(vv0 +dv*(x2-1.0))
!
      return
      end function dgaus2
!
!
!---------------------------------------------------------------
      function fun (v)
!---------------------------------------------------------------
!  exp(-v**2/2)
      implicit none
      real*4 fun,v
!
      fun= exp(-v**2/2.)
!
      return
      end function fun
!
!
!-------------------------------------------------------------
      subroutine vdistr (xg,yg,zg,vx,vy,vz,ns,np,nCLp)
!-------------------------------------------------------------
      implicit none
      include  'param_em3p7_Ca.h'
!
      real*8     xg(npq0),yg(npq0),zg(npq0),   &
                 vx(npq0),vy(npq0),vz(npq0),vv
      integer*4  ns,np,nCLp,ILN,i,k,ix,iy,iz
!
      real*4     fvx(101),fvy(101),fvz(101),fvs(101),xsc(101), &
                 vmax2,aiv,fmax1,fmin1
!
      real*8        pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
!
!-------------------------
!* Velocity distribution
!-------------------------
!* protons
!
      if(np.eq.0) go to 1000
! ---
      vv= 0
      do i= ns+1,ns+np
      vv= dmax1(vv,vx(i)**2+vy(i)**2+vz(i)**2)
      end do
!
      vmax2= dmax1(sqrt(vv),1.d-5)  ! cm/sec
      aiv= 50/vmax2
!
      do k= 1,101
      xsc(k)= (k -51)/aiv
      fvx(k)= 0
      fvy(k)= 0
      fvz(k)= 0
      end do
! 
      do i= ns+1,ns+np
      ix= aiv*vx(i) +51.
      iy= aiv*vy(i) +51.
      iz= aiv*vz(i) +51.
!
      if(iabs(ix-51).le.50) fvx(ix)= fvx(ix) +1. 
      if(iabs(iy-51).le.50) fvy(iy)= fvy(iy) +1. 
      if(iabs(iz-51).le.50) fvz(iz)= fvz(iz) +1. 
      end do
!
!* Average over x,y,z
      do k= 1,101
      fvs(k)= (fvx(k) +fvy(k))/2.
      end do
!
      call lplmax (fvs,fmax1,fmin1,101)
      ILN= 1
      call hplot1 (2,4,101,xsc,fvs,fmax1,fmin1,ILN,'FVxy(H+)',8, &
                   '  VR    ',8,'        ',8)
!
      call lplmax (fvz,fmax1,fmin1,101)
      call hplot1 (3,4,101,xsc,fvz,fmax1,fmin1,ILN,'FVz(H+) ',8, &
                   '  VZ    ',8,'        ',8)
 1000 continue
!
!* Electrons
      vv= 0 
      do i= ns+np+1,nCLp
      vv= dmax1(vv,vx(i)**2+vy(i)**2+vz(i)**2)
      end do
!
      vmax2= dmax1(sqrt(vv),1.d-5)  ! cm/sec
      aiv= 50/vmax2 
!
      do k= 1,101
      xsc(k)= (k -51)/aiv
      fvx(k)= 0
      fvy(k)= 0
      fvz(k)= 0
      end do
!
      do i= ns+np+1,nCLp
      ix= aiv*vx(i) +51.
      iy= aiv*vy(i) +51.
      iz= aiv*vz(i) +51.
!
      if(iabs(ix-51).le.50) fvx(ix)= fvx(ix) +1. 
      if(iabs(iy-51).le.50) fvy(iy)= fvy(iy) +1. 
      if(iabs(iz-51).le.50) fvz(iz)= fvz(iz) +1. 
      end do
!
!* Average over x,y
      do k= 1,101
      fvs(k)= (fvx(k) +fvy(k))/2.
      end do
!
      call lplmax (fvs,fmax1,fmin1,101)
      ILN= 1
      call hplot1 (2,5,101,xsc,fvs,fmax1,fmin1,ILN,'FVxy(el)',8, &
                   '  VR    ',8,'        ',8)
!
      call lplmax (fvz,fmax1,fmin1,101)
      call hplot1 (3,5,101,xsc,fvz,fmax1,fmin1,ILN,'FVz(el) ',8, &
                   '  VZ    ',8,'        ',8)
!
!* CNT
      vv= 0 
      do i= 1,ns
      vv= dmax1(vv,vx(i)**2+vy(i)**2+vz(i)**2)
      end do
!
      vmax2= dmax1(sqrt(vv),1.d-5)
      aiv= 50/vmax2  ! no dimension
!
      do k= 1,101
      xsc(k)= (k -51)/aiv
      fvx(k)= 0
      fvy(k)= 0
      fvz(k)= 0
      end do
!
      do i= 1,ns
      ix= aiv*vx(i) +51.
      iy= aiv*vy(i) +51.
      iz= aiv*vz(i) +51.
!
      if(iabs(ix-51).le.50) fvx(ix)= fvx(ix) +1. 
      if(iabs(iy-51).le.50) fvy(iy)= fvy(iy) +1. 
      if(iabs(iz-51).le.50) fvz(iz)= fvz(iz) +1. 
      end do
!
!* Average over x,y
      do k= 1,101
      fvs(k)= (fvx(k) +fvy(k))/2.
      end do
!
      call lplmax (fvs,fmax1,fmin1,101)
      ILN= 1
      call hplot1 (2,6,101,xsc,fvs,fmax1,fmin1,ILN,'FVxy(CA)',8, &
                   '  VR    ',8,'        ',8)
!
      call lplmax (fvz,fmax1,fmin1,101)
      call hplot1 (3,6,101,xsc,fvz,fmax1,fmin1,ILN,'FVz(CAu)',8, &
                   '  VZ    ',8,'        ',8)
!    ------------
      call chart
!    ------------
      return
      end subroutine vdistr
!
!
!-----------------------------------------------------------
      subroutine edistr (vx,vy,vz,am,ns,np,nCLp)
!-----------------------------------------------------------
      implicit none
      include  'param_em3p7_Ca.h'
!
      integer*4  ns,np,nq,nCLp,i,k,ix,ILN
      real*8     vx(npq0),vy(npq0),vz(npq0),am(npq0),v7
      real*4     vv,vv1,vv2,vv3,vv4,vmax1,vmax2,aiv1,aiv2,  &
                 fmax1,fmin1,fmax2,fmin2
!
      real*4     fe(101),fi(101),vsc1(101),vsc2(101)
      real*4     time,xp_leng
      COMMON/HEADR2/ time,xp_leng
!-------------------------
!* Energy distribution
!-------------------------
!
      if(np.eq.0) return
! ---
      v7= 0
      do i= ns+1,ns+np
      v7= dmax1(v7,am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)) ! cm/sec
      end do
      vv= 0.5d0*v7
      vmax1= amax1(vv,1.e-5)
!
      v7= 0
      do i= ns+np+1,nCLp
      v7= dmax1(v7,am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2))
      end do
      vv= 0.5d0*v7
      nq= nCLP -(ns+np)
      vmax2= amax1(vv,1.e-5) ! cm/sec
!
      aiv1= 100/vmax1 
      aiv2= 100/vmax2
!
      do k= 1,101
      vsc1(k)= (k-1)/aiv1
      vsc2(k)= (k-1)/aiv2
      fi(k)= 0
      fe(k)= 0
      end do
!                   ++
      do i= ns+1,ns+np
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)  ! cm/sec
      vv= 0.5d0*v7
      ix= aiv1*vv +1.01
      if(ix.gt.0 .and. ix.le.101) fi(ix)= fi(ix) +1.
      end do
!                ++
      do i= ns+np+1,nCLp
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      vv= 0.5d0*v7
      ix= aiv2*vv +1.01
      if(ix.gt.0 .and. ix.le.101) fe(ix)= fe(ix) +1.
      end do
!
      call lplmax (fi,fmax1,fmin1,101)
      call lplmax (fe,fmax2,fmin2,101)
!
      ILN= 1
      call hplot1 (2,2,101,vsc1,fi,fmax1,fmin1,ILN,'F(H+)   ',8, &
                   ' Lin(E) ',8,'        ',8)
      call hplot1 (2,3,101,vsc2,fe,fmax2,fmin2,ILN,'F(ele)  ',8, &
                   ' Lin(E) ',8,'        ',8)
!
!* Log-Log plot
      vv1= -1000
      vv2=  1000
      vv3= -1000
      vv4=  1000
!
!  Protons
      do i= ns+1,ns+np
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      if(v7.gt.0.d0) then
        vv= alog10(real(0.5d0*v7))
      else
        vv= 0.d0
      end if
      vv1= amax1(vv1,vv)
      vv2= amin1(vv2,vv)
      end do
!
!  Electrons
      do i= ns+np+1,nCLp
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      if(v7.gt.0.d0) then
        vv= alog10(real(0.5d0*v7))
      else
        vv= 0.d0
      end if
      vv3= amax1(vv3,vv)
      vv4= amin1(vv4,vv)
      end do
!
      vv1= vv1 +0.5
      vv3= vv3 +0.5
!
      vmax1= amax1((vv1 -vv2),1.e-5)
      vmax2= amax1((vv3 -vv4),1.e-5)
!
      aiv1= 100/vmax1
      aiv2= 100/vmax2
!
      do k= 1,101
      vsc1(k)= vv2 +(k-1)/aiv1
      vsc2(k)= vv4 +(k-1)/aiv2
!
      fi(k)= 0
      fe(k)= 0
      end do
!
      do i= ns+1,ns+np
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      if(v7.gt.0.d0) then
        vv= alog10(real(0.5d0*v7))
      else
        vv= 0.d0
      end if
      ix= aiv1*(vv -vv2) +1.01
      if(ix.gt.0 .and. ix.le.101) fi(ix)= fi(ix) +1.
      end do
!
      do i= ns+np+1,nCLp
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      if(v7.gt.0.d0) then
        vv= alog10(real(0.5d0*v7))
      else
        vv= 0.d0
      end if
      ix= aiv2*(vv -vv4) +1.01
      if(ix.gt.0 .and. ix.le.101) fe(ix)= fe(ix) +1.
      end do
!
      call lplmax (fi,fmax1,fmin1,101)
      call lplmax (fe,fmax2,fmin2,101)
!
      ILN= 1
      call hplot1 (3,2,101,vsc1,fi,fmax1,fmin1,ILN,'Log.FH+ ',8, &
                   ' Log(E) ',8,'        ',8)
      call hplot1 (3,3,101,vsc2,fe,fmax2,fmin2,ILN,'Log.Fele',8, &
                   ' Log(E) ',8,'        ',8)
!    ------------
      call chart
!    ------------
!
      return
      end subroutine edistr
!
!
!------------------------------------------------------
      subroutine averg1 (q,qav,is)
!------------------------------------------------------
      implicit none
      real*4 q(5000),qav
      integer*4  is,i
!
      qav= 0
!
      do 100 i= is-9,is
      qav= qav +q(i)
  100 continue
!
      qav= qav/10.e0
!
      return
      end subroutine averg1
!
!
!------------------------------------------------------
      subroutine rehist
!------------------------------------------------------
      implicit none
      include  'param_em3p7_Ca.h'
!
      real*4        ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,          &
                    Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,            &
                    sex,sez,sbx,sbz,time
      common/ehist/ ekin(3000),ppot(3000),ekn1(3000),ekn2(3000), &
                    etot(3000),Rgyi(3000),Rgye(3000),Rgyf(3000), &
                    Rgyc(3000),Rion(3000),Nele(3000),dPot(3000), &
                    ecr(3000),elj(3000),sex(3000),sez(3000),     &
                    sbx(3000),sbz(3000),time(3000)
!
      integer*4     it,is,k,k1
      common/parm1/ it,is
!
      integer*4       iiwrt1,iiwrt2,iiwrt3,iiwrt4
      common/iiwrite/ iiwrt1,iiwrt2,iiwrt3,iiwrt4
!
      real*4        phi,tht,rgmax
      common/parm4/ phi,tht,rgmax
!
      do k= 1,is
      if(mod(k,2).eq.0) then
!
      k1= k/2
      ekin(k1)= (ekin(k-1) +ekin(k))/2
      ppot(k1)= (ppot(k-1) +ppot(k))/2
      ekn1(k1)= (ekn1(k-1) +ekn1(k))/2
      ekn2(k1)= (ekn2(k-1) +ekn2(k))/2
      etot(k1)= (etot(k-1) +etot(k))/2
      Rgyi(k1)= (Rgyi(k-1) +Rgyi(k))/2
      Rgye(k1)= (Rgye(k-1) +Rgye(k))/2
!
       Rgyc(k1)= ( Rgyc(k-1) + Rgyc(k))/2
      Rgyf(k1)= (Rgyf(k-1) +Rgyf(k))/2
       ecr(k1)= ( ecr(k-1) + ecr(k))/2
       elj(k1)= ( elj(k-1) + elj(k))/2
      Rion(k1)= (Rion(k-1) +Rion(k))/2
      Nele(k1)= (Nele(k-1) +Nele(k))/2
      dPot(k1)= (dPot(k-1) +dPot(k))/2
!
       sex(k1)= (sex(k-1) + sex(k))/2
       sez(k1)= (sez(k-1) + sez(k))/2
       sbx(k1)= (sbx(k-1) + sbx(k))/2
       sbz(k1)= (sbz(k-1) + sbz(k))/2
      
      time(k1)= (time(k-1) +time(k))/2
      end if
      end do
!
      is  = is/2
      iiwrt1= 2*iiwrt1
!
      if(ionode) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Rehist: new is, iiwrt1=',is,iiwrt1
        close(11)
      end if
!
      return
      end subroutine rehist
!
!
!------------------------------------------------------
      subroutine lplots
!------------------------------------------------------
      implicit none
!
      real*4        ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,          &
                    Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,            &
                    sex,sez,sbx,sbz,time,ekin0,etot1(3000)
      common/ehist/ ekin(3000),ppot(3000),ekn1(3000),ekn2(3000), &
                    etot(3000),Rgyi(3000),Rgye(3000),Rgyf(3000), &
                    Rgyc(3000),Rion(3000),Nele(3000),dPot(3000), &
                    ecr(3000),elj(3000),sex(3000),sez(3000),     &
                    sbx(3000),sbz(3000),time(3000)
!
      character*8    label,date_now*10
      real*4         t,xp_leng
      integer*4      it,is
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ t,xp_leng
      common/parm1/  it,is
!
      real*4       emax0,emin0,emax1,emin1,emax2,emin2,          &
                   pmax,pmin,elmax,elmin,emax,etot0,             &
                   etmax,etmin,Rgyi1,Rgyi2,RTmax,RTmin,rmax,rmin,&
                   e1max,e1min,e2max,e2min,e3max,e3min,e4max,e4min
      integer*4    ILN,ILG,i
!
!     ekin(1)= ekin(3)
!     ekin(2)= ekin(3)
!     etot(1)= etot(3)
!     etot(2)= etot(3)
!
      if(is.lt.5) return  ! too few a case
!
      call lplmax (ekin,emax0,emin0,is)
      call lplmax (ekn1,emax1,emin1,is)
      call lplmax (ekn2,emax2,emin2,is)
      call lplmax ( ecr,pmax,pmin,is)
      call lplmax ( elj,elmax,elmin,is)
      ekin0= (ekin(1)+ekin(2)+ekin(3))/3
      emax= amax1(emax0,emax1,emax2)
!
      ILN= 1
      ILG= 2
!
      call lplot1 (2,4,is,time,ekin,emax0,0.0,ILN,'KIN/C+Au',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,5,is,time,ekn1,emax1,0.0,ILN,'KIN/H(+)',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,ekn2,emax2,0.0,ILN,'KIN/ele ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,4,is,time, ecr,pmax,pmin,ILN,'E_CL(in)',8, &
                 '  time  ',8,'        ',8)
!
      call lplot1 (3,5,is,time, elj,elmax,elmin,ILN,'E_LJ(in)',8, &
                 '        ',8,'        ',8)
!
      do i= 1,is
      etot1(i)= ekin(i)+ekn1(i)+ekn2(i)+ecr(i)+elj(i) 
      end do
!
      etot0= etot1(1)
      do i= 1,is
      etot1(i)= etot1(i) -etot0
      end do
!
      call lplmax (etot1,etmax,etmin,is)
      call lplot1 (3,6,is,time,etot1,etmax,etmin,ILN,'dE.TOTAL',8, &
                 '  time  ',8,'        ',8)
!------------------------
      CALL CHART
!------------------------
!
! page 2
      call lplmax (sex,e1max,e1min,is)
      call lplmax (sez,e2max,e2min,is)
      call lplmax (sbx,e3max,e3min,is)
      call lplmax (sbz,e4max,e4min,is)
      emax1= amax1(e1max,e2max)
      emax2= amax1(e3max,e4max)
      emin1= amin1(e1min,e2min)
      emin2= amin1(e3min,e4min)
!
      ILN= 1
      ILG= 2
!
      call lplot1 (2,4,is,time,sex,emax1,emin1,ILN,'sex+sez ',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,5,is,time,sez,emax1,emin1,ILN,'sez     ',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,sbx,emax2,emin2,ILN,'sbx+sbz ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,4,is,time,sbz,emax2,emin2,ILN,'sez     ',8, &
                 '        ',8,'        ',8)
!
!------------------------
      CALL CHART
!------------------------
!
! page 3
      CALL SYMBOL (3.0,3.0,0.7,'Wion=', 0.,5)
      call number (6.0,3.0,0.7,ekin0,0.,5)
      CALL SYMBOL (3.0,2.0,0.7,'W0=', 0.,3)
      call number (6.0,2.0,0.7,etot0,0.,5)
!
      call lplmax (Rgyi,Rgyi1,Rgyi2,is)
      call lplot1 (2,4,is,time,Rgyi,Rgyi1,0.0,ILN,'Rgy(H+) ',8, &
                   '        ',8,'        ',8)
!
      call lplmax (Rgyc,RTmax,RTmin,is)
      call lplot1 (2,5,is,time,Rgyc,RTmax,0.0,ILN,'Rgy(CAu)',8, &
                   '        ',8,'        ',8)
!
      call lplmax (Rgyf,rmax,rmin,is)
      call lplot1 (2,6,is,time,Rgyf,rmax,0.0,ILN,'Rgy(e)  ',8, &
                   '  time  ',8,'        ',8)
!
      call lplmax (Rion,rmax,rmin,is)
      call lplot1 (3,4,is,time,Rion,rmax,0.0,ILN,'Rmax(H+)',8, &
                   '        ',8,'        ',8)
!
      call lplmax (Nele,rmax,rmin,is)
      call lplot1 (3,5,is,time,Nele,rmax,0.0,ILN,'Ne(r<Ri)',8, &
                   '        ',8,'        ',8)
!
      call lplmax (dPot,rmax,rmin,is)
      call lplot1 (3,6,is,time,dPot,rmax,rmin,ILN,'d.Pot   ',8, &
                   '  time  ',8,'        ',8)
!------------------------
      CALL CHART
!------------------------
!
      return
      end subroutine lplots
!
!
!------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
!------------------------------------------------------
      implicit none
!
      real*4 f(is),fmax,fmin
      integer*4  is,i
!
      fmax= -1.d10
      fmin=  1.d10
!
      do 100 i= 1,is
      fmax= amax1(fmax,f(i))
      fmin= amin1(fmin,f(i))
  100 continue
!
      return
      end subroutine lplmax
!
!
!------------------------------------------------------------------
      subroutine ppl3da (xg,yg,zg,ch,ag,Rgy1,Rgy2,ns,np,nq,  &
                         nCLp,igrp) !,ifskip_e,ifskip_p)
!------------------------------------------------------------------
      include    'param_em3p7_Ca.h'
!
      integer*4  ns,np,nq,nCLp,igrp,nZ,nZA
      common/iroha/ nZ,nZA
!  
      real*8     xg(npq0),yg(npq0),zg(npq0),ch(npq0),ag(npq0)
      real*4     Rgy1,Rgy2,rgmax,rmax1,ypp,zpp,amass
!   
      logical    ifskip_e,ifskip_p
!
      real*8        pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax,  &
                    Temp,epsCLJ,epsLJ
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm4/ phi,tht,rgmax
      COMMON/ELSTA/ Temp,epsCLJ,epsLJ
!
      real*8         R_sp,D_sp,N_sp,Lambda_D,massi,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
!
      CHARACTER*8    label,date_now*10,cax*1
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ t,xp_leng
!
      HH= 0.7
!     CALL SYMBOL (0.5,18.0,HH,LABEL,0.,8)
!     CALL SYMBOL (4.5,18.0,HH,praefixc//'.77'//suffix(igrp),0.,28)
!
      CALL SYMBOL ( 0.5,1.0,HH,date_now, 0.,10)
      CALL SYMBOL (17.0,1.0,HH,'t=', 0.,2)
      call number (20.0,1.0,HH,t,0.,7)
!
      Temp4 = Temp
      Dens4 = D_sp
      Rsp4  = R_sp  ! in micron
      Nsp4  = N_sp
!     CALL SYMBOL (16.4,16.0,HH,'Temp=', 0.,5)
!     call number (18.5,16.0,HH,Temp4,0.,5)
!     CALL SYMBOL (16.4,15.2,HH,'Dens=', 0.,5)
!     call number (18.5,15.2,HH,Dens4,0.,5)
!     CALL SYMBOL (16.4,14.4,HH,'Rsp=', 0.,4)
!     call number (18.5,14.4,HH,Rsp4,0.,5)
!
!     CALL SYMBOL (16.4,13.4,HH,'Nsp=', 0.,4)
!     call number (18.5,13.4,HH,float(Nsp4),0.,5)
!     CALL SYMBOL (16.4,12.6,HH,'np=', 0.,3)
!     call number (18.5,12.6,HH,float(np),0.,5)
!
!     CALL SYMBOL (16.4,10.8,HH,'Rgi=', 0.,4)
!     call number (18.5,10.8,HH,Rgy1,0.,5)
!     CALL SYMBOL (16.4,10.0,HH,'Rge=', 0.,4)
!     call number (18.5,10.0,HH,Rgy2,0.,5)
!
      FSIZE= 13.
      HL=  10.
      VD=   8.
!
      pi=  4.*atan(1.0)
      pha= pi*phi/180.
      tha= pi*tht/180.
!
      cph= cos(pha)
      sph= sin(pha)
      cth= cos(tha)
      sth= sin(tha)
!
      rmax1= 0.
!
!     do 100 i= 1,nCLp
      do 100 i= 1,ns+np   ! for all atoms
!     if(abs(x(i)).gt.xmax) go to 100
!     if(abs(y(i)).gt.ymax) go to 100
!     if(abs(z(i)).gt.zmax) go to 100
!
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      rmax1= amax1(rmax1, sqrt(ypp**2+zpp**2))
!     rmax1= 2.8d-2   ! fixed in time
  100 continue
!
!c    rmax1= dmax1(1.d0, rmax1) ! dmax1(rmax1,xmax)
      ps= 0.5*fsize/rmax1
!
!**********************
!*  Draw Axes.        *
!**********************
!
      amass= massi
      CALL SYMBOL (HL+6.5,VD-3.5,HH,'Rmax=', 0.,5)
      call number (HL+9.0,VD-3.5,HH,rmax1,0.,5)
      CALL SYMBOL (HL+6.5,VD-4.3,HH,'massp=', 0.,6)
      call number (HL+9.2,VD-4.3,HH,amass,0.,5)
!
      do 200 i= 1,3
      if(i.eq.1) then
         x1= 0.7*rmax1
         y1= 0.
         z1= 0.
         cax='X'
      else if(i.eq.2) then
         x1= 0.
         y1= 0.7*rmax1
         z1= 0.
         cax='Y'
      else if(i.eq.3) then
         x1= 0.
         y1= 0.
         z1= 0.7*rmax1
         cax='Z'
      end if
!
      xp= x1*cph -y1*sph
      yp= x1*sph +y1*cph
      zp= z1
!
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp  +HL
      yy= ps*zpp  +VD
      CALL PLOT (HL,VD,3)
      CALL PLOT (xx,yy,2)
!
      CALL SYMBOL (xx-0.7,yy-0.5,HH,cax,0.,1)
  200 continue
!
!
!* Protons
!
!     if(.not.ifskip_p) then
!---
      nskip= max(1,np/1000)
!---                  ****
!
      do 300 i= ns+1,ns+np,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 300
      if(yy.lt.0.1 .or. yy.gt.18.) go to 300
!
      dd= 7*ps*ag(i)
      call newcolor (3,0.,0.,1.)  ! blue
      call circle (xx-0.12,yy-0.12,dd,2)
  300 continue
!---
!     end if
!
!* carbon ions
!
      nskip= (ns/2)/1000  ! plot only a few
!
      do 400 i= 1,ns/2,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 400
      if(yy.lt.0.1 .or. yy.gt.18.) go to 400
!
      dd= 7*ps*ag(i)
      call newcolor (3,0.,1.,0.)  ! green
      call circle (xx-0.12,yy-0.12,dd,2)
  400 continue
!
!* Au ions
!
      nskip= (ns/2)/1000
!---                ****
!
      do 430 i= ns/2+1,ns,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 430
      if(yy.lt.0.1 .or. yy.gt.18.) go to 430
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.0,0.7,0.0)  ! gold
      call circle (xx-0.12,yy-0.12,dd,2)
  430 continue
!
!* Electrons in pellet
!
!     if(.not.ifskip_e) then
!---
      nskip= max(1,np/1000)   ! -e
!---                  ****
!
      do 500 i= ns+np+1,ns+2*np,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 500
      if(yy.lt.0.1 .or. yy.gt.18.) go to 500
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.,0.,0.)  ! red
      call circle (xx-0.12,yy-0.12,dd,1)
  500 continue
!
!* Electrons around CNT - heavy
!
!     nskip= 6*ns/1000
      nskip= (nZ+nZA)*(ns/2)/1000
!---                       ****
!
      do 600 i= ns+2*np+1,nCLp,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 600
      if(yy.lt.0.1 .or. yy.gt.18.) go to 600
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.,0.,0.)  ! red
      call circle (xx-0.12,yy-0.12,dd,1)
  600 continue
!---
!     end if
!
      call newcolor (0,1.,0.,0.)  ! reset
!---------------------
      CALL CHART
!---------------------
!
      return
      end
!
!
!------------------------------------------------------------------
      subroutine ppl3dz (xg,yg,zg,ch,ag,Rgy1,Rgy2,ns,np,nq,  &
                         nCLp,igrp)
!------------------------------------------------------------------
!  View from top (x,y)
!
      include    'param_em3p7_Ca.h'
!
      integer*4  ns,np,nq,nCLp,igrp,nZ,nZA
      common/iroha/ nZ,nZA
!  
      real*8     xg(npq0),yg(npq0),zg(npq0),ch(npq0),ag(npq0)
      real*4     Rgy1,Rgy2,rgmax,rmax1,ypp,zpp
!  
!
      real*8        pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax,  &
                    Temp,epsCLJ,epsLJ
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm4/ phi,tht,rgmax
      COMMON/ELSTA/ Temp,epsCLJ,epsLJ
!
      real*8         R_sp,D_sp,N_sp,Lambda_D,massi,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
!
      CHARACTER*8    label,date_now*10,cax*1
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ t,xp_leng
!
      HH= 0.7
!     CALL SYMBOL (0.5,18.0,HH,LABEL,0.,8)
!     CALL SYMBOL (4.5,18.0,HH,praefixc//'.77'//suffix(igrp),0.,28)
!
      CALL SYMBOL ( 0.5,1.0,HH,date_now, 0.,10)
      CALL SYMBOL (17.0,1.0,HH,'t=', 0.,2)
      call number (20.0,1.0,HH,t,0.,7)
!
      Temp4 = Temp
      Dens4 = D_sp
      Rsp4  = R_sp  ! in micron
      Nsp4  = N_sp
!     CALL SYMBOL (16.4,16.0,HH,'Temp=', 0.,5)
!     call number (18.5,16.0,HH,Temp4,0.,5)
!     CALL SYMBOL (16.4,15.2,HH,'Dens=', 0.,5)
!     call number (18.5,15.2,HH,Dens4,0.,5)
!     CALL SYMBOL (16.4,14.4,HH,'Rsp=', 0.,4)
!     call number (18.5,14.4,HH,Rsp4,0.,5)
!
!     CALL SYMBOL (16.4,13.4,HH,'Nsp=', 0.,4)
!     call number (18.5,13.4,HH,float(Nsp4),0.,5)
!     CALL SYMBOL (16.4,12.6,HH,'np=', 0.,3)
!     call number (18.5,12.6,HH,float(np),0.,5)
!
!     CALL SYMBOL (16.4,10.8,HH,'Rgi=', 0.,4)
!     call number (18.5,10.8,HH,Rgy1,0.,5)
!     CALL SYMBOL (16.4,10.0,HH,'Rge=', 0.,4)
!     call number (18.5,10.0,HH,Rgy2,0.,5)
!
      FSIZE= 13.
      HL=  10.
      VD=   8.
!
      pi=  4.*atan(1.0)
      pha= pi*phi/180.
      tha= pi*tht/180.
!
      cph= cos(pha)
      sph= sin(pha)
      cth= cos(tha)
      sth= sin(tha)
!
      rmax1= 0.
!
!     do 100 i= 1,nCLp
      do i= 1,ns+np  ! for all atoms
      xx = - xg(i)     ! arrange the x-direction with ppl3da plot
      yy = yg(i)
!
!     if(abs(x(i)).gt.xmax) go to 100
!     if(abs(y(i)).gt.ymax) go to 100
!
      rmax1= amax1(rmax1, sqrt(xx**2 +yy**2))
!     rmax1= 2.11d-2     ! fixed in time
      end do
!
      ps= 0.5*fsize/rmax1
!
!**********************
!*  Draw Axes.        *
!**********************
!
      rmax4= rmax1
      CALL SYMBOL (HL+6.5,VD-3.5,HH,'Rmax=', 0.,5)
      call number (HL+9.0,VD-3.5,HH,rmax4,0.,5)
!
      do i= 1,2
      if(i.eq.1) then
         x1= - 0.7*rmax1
         y1= 0.
         cax='X'
      else if(i.eq.2) then
         x1= 0.
         y1= 0.7*rmax1
         cax='Y'
      end if
!
      xx= ps*x1 +HL
      yy= ps*y1 +VD
      CALL PLOT (HL,VD,3)
      CALL PLOT (xx,yy,2)
!
      CALL SYMBOL (xx-0.7,yy+0.2,HH,cax,0.,1)
      end do
!
!
!* Protons
!
      nskip= max(1,np/1000)
!
      do i= ns+1,ns+np,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 300
      if(yy.lt.0. .or. yy.gt.18.) go to 300
!
      dd= 7*ps*ag(i)
      call newcolor (3,0.,0.,1.)  ! blue
      call circle (xx-0.12,yy-0.12,dd,2)
  300 end do
!
!* carbon ions
!
      nskip= (ns/2)/1000  ! plot only a few
!
      do i= 1,ns/2,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 400
      if(yy.lt.0. .or. yy.gt.18.) go to 400
!
      dd= 7*ps*ag(i)
      call newcolor (3,0.0,1.0,0.0)  ! green
      call circle (xx-0.12,yy-0.12,dd,2)
  400 end do
!
!* Au ions
!
      nskip= (ns/2)/1000  ! plot only a few
!
      do i= ns/2+1,ns,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 430
      if(yy.lt.0. .or. yy.gt.18.) go to 430
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.0,0.7,0.0)  ! gold
      call circle (xx-0.12,yy-0.12,dd,2)
  430 end do
!
!* Electrons in pellet
!
      nskip= max(1,np/1000)
!
      do 500 i= ns+np+1,ns+np+np,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 500
      if(yy.lt.0. .or. yy.gt.18.) go to 500
!
      dd= 7*ps*ag(i)
!     call newcolor (3,1.,1.,0.)  ! orange
      call newcolor (3,1.,0.,0.)  ! red
      call circle (xx-0.12,yy-0.12,dd,1)
  500 end do
!
!* Electrons around CNT
!
!     nskip= 6*ns/1000
      nskip= (nZ+nZA)*(ns/2)/1000   ! (ns/2)*(6+20)
!  
!---                       ****
!
      do 600 i= ns+np+np+1,nCLp,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 600
      if(yy.lt.0. .or. yy.gt.18.) go to 600
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.,0.,0.)  ! red
      call circle (xx-0.12,yy-0.12,dd,1)
  600 end do
!
      call newcolor (0,1.,0.,0.)  ! reset
!---------------------
      CALL CHART
!---------------------
!
      return
      end
!
!
!---------------------------------------
       subroutine newcolor (ic,r,g,b)
!---------------------------------------
!  ic= 3 tri-color
!  ic= 0 gray scale, r= 0. for black
!
       write(77,*) 'stroke'
!
       if(ic.eq.0) then
         write(77,10) 1.-r  ! 0. for black
   10    format(f4.1,' setgray')
       end if
!
       if(ic.eq.3) then
         write(77,30) r,g,b
   30    format(3f4.1,' setrgbcolor')
       end if
!
       return
       end
!
!
!-------------------------------------------------
       subroutine circle (x,y,d,ic)
!-------------------------------------------------
!*  Open circle centered at (x,y) /or outer edge.
!
      write(77,*) " 3.0 setlinewidth"
!
      pi= 3.1415927
      nc= 13
      dth= 2.*pi/nc
      a= d/2.
!
      x0= x +a
      y0= y
      call plot (x0,y0,3)
!
      do 100 j= 1,nc
      th= dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      call plot (x1,y1,2)
  100 continue
!
      call plot (x1,y1,3)
      write(77,*) " 1.0 setlinewidth"
!
      if(ic.eq.1) return
!------------------------------------
!*  Filled circle centered at (x,y).
!------------------------------------
!
      write(77,*) " 3.0 setlinewidth"
!
      nc= 5
      dth= pi/(2*nc +1)
!
      do 300 j= -nc,nc
      th= 0.5*pi +dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      x2= x1
      y2= 2.*y -y1
!
      call plot (x1,y1,3)
      call plot (x2,y2,2)
  300 continue
!
      call plot (x2,y2,3)
      write(77,*) " 1.0 setlinewidth"
!
      return
      end
!
!
!------------------------------------------------
      integer function iwrta (t,twr)
!------------------------------------------------
      common/imemo/ iwa,iwb,iwc,iwd
!
      iw= t/twr 
      if(iw.gt.iwa) then
        iwa= iw
        iwrta= 0
      else
        iwrta= 1
      end if
!
      return
      end
!
!
!------------------------------------------------
      integer function iwrtb (t,twr)
!------------------------------------------------
      common/imemo/ iwa,iwb,iwc,iwd
!
      iw= t/twr 
      if(iw.gt.iwb) then
        iwb= iw
        iwrtb= 0
      else
        iwrtb= 1
      end if
!
      return
      end
!
!
!------------------------------------------------
      integer function iwrtc (t,twr)
!------------------------------------------------
      common/imemo/ iwa,iwb,iwc,iwd
!
      iw= t/twr 
      if(iw.gt.iwc) then
        iwc= iw
        iwrtc= 0
      else
        iwrtc= 1
      end if
!
      return
      end
!
!
!------------------------------------------------
      integer function iwrtd (t,twr)
!------------------------------------------------
      common/imemo/ iwa,iwb,iwc,iwd
!
      iw= t/twr 
      if(iw.gt.iwd) then
        iwd= iw
        iwrtd= 0
      else
        iwrtd= 1
      end if
!
      return
      end
!
!
!------------------------------------------------
      block data
!------------------------------------------------
      common/ranfff/ ir,iq
      data  ir/3021/,iq/7331/    ! original
!     data  ir/17331/,iq/37711/
      end
!
!
!------------------------------------------------
      function ranff (x)
!------------------------------------------------
!*  ranf= (0,1)
!
      common/ranfff/ ir,iq
!
      REAL*8     ranff,x,INVM
      PARAMETER  (MASK=2**30+(2**30-1),INVM= 0.5D0**31)
      PARAMETER  (LAMBDA=48828125)
!
      IR= IAND( LAMBDA*IR, MASK)
      ranff= IR*INVM
!
!     ask= 371597.
!     ambda= sqrt(ask)
!     qq= 0.3713*ask
!
!     ir= amod( ambda*ir +qq, ask)
!     ranff= ir/ask
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine lplot1 (ix,iy,npt1,x,y,ymax,ymin,il,lab1,n1,lab2,n2,  &
                         lab3,n3)
!-----------------------------------------------------------------------
!  <<warning>>  order and number of arguments /lplot/ have been changed.
!               also, x (time) is defined for all range.
!               date: 5/18/96 at mit.
!***********************************************************************
!   il=1................ linear plot of (x,y)
!   il=2................ log10 plot of (x,log y)
!***********************************************************************
!
      dimension  x(npt1),y(npt1),u(npt1),v(npt1)
      dimension  xcm(6),ycm(6),pl(6),pr(6),ql(6),qr(6)
!
      character*8    lab1,lab2,lab3
      character*8    label,date_now*10,cax*1
      common/headr1/ label,date_now
      common/headr2/ time,xp_leng
      common/pplcom/ nfine,pl1(10),pr1(10),ql1(10),qr1(10),  &
                     xmin1(10),xmax1(10),ymin1(10),ymax1(10)
!
!   for fujitsu.
!     data  xcm/18.46,2*9.867,3*6.18/,
!    *      ycm/16.85,2*7.435,3*4.381/,
!    *      pl/2*2.00,15.132,2.00,8.00,18.20/,
!    *      ql/1.95,10.885,1.95,13.832,7.891,1.95/
!
!   for nec.
      data  xcm/21.0, 2*10.00, 3*7.00/,  &
            ycm/15.0, 2*6.80, 3*3.90/,   &
            pl/2.0,  2.0,14.0, 1.0,9.0,17.0/, &
            ql/2.3, 10.5,2.3, 12.9,7.6,2.3/
      logical  lab_skip
!
      iplot=1
      go to 1
!
!-----------------------------------------------------------------------
      entry hplot1 (ix,iy,npt1,x,y,ymax,ymin,il,lab1,n1,lab2,n2,lab3,n3)
!-----------------------------------------------------------------------
      iplot=2
!
    1 npt= npt1
      isc= 1
!
      do i=1,6
      pr(i)= pl(i) +xcm(i)
      end do
!
      do j=1,6
      qr(j)= ql(j) +ycm(j)
      end do
!
      lab_skip= .false.
      if(il.eq.7) lab_skip= .true.
!
!                 ******************************************************
!*                **  make a copy before the top-left frame is drawn. **
!                 ******************************************************
      hh = 0.70
      hhs= 0.60
!
      i1= iabs(ix)
      j1= iabs(iy)
      if(i1.ge.3) go to 10
      if(j1.eq.3.or.j1.ge.5) go to 10
!                                              ************************
!                                              ** label of the page. **
!                                              ************************
      call symbol (0.1,18.0,hh,label,0.,8)
      call symbol (3.1,18.0,hh,date_now, 0.,10)
      call symbol (15.9,0.1,hh,'t =',0.,3)
      call number (999.0,999.0,hh,time,0.,5)
!
   10 continue
!
      do i=1,npt
      u(i)= x(i)
      end do
!
      xmax= u(npt)
      xmin= u(1)
!                             ************************************
!                             ** three-point average if il > 0  **
!                             ************************************
      if(il.gt.0) then
        v(1)=   y(1)
        v(npt)= y(npt)
        do i=2,npt-1
        v(i)= y(i)
!       v(i)= 0.33333*(y(i-1)+y(i)+y(i+1))
        end do
      else
        do i=1,npt
        v(i)= y(i)
        end do
      end if
!                                                *****************
!                                                **  log. scale **
!                                                *****************
      if(iabs(il).eq.2) then
         do i=1,npt
         if(v(i).gt.0.) then
            v(i)= alog10(v(i))
         else
            v(i)= -10.
         end if
         end do
      end if
!                                **************************************
!                                ** set a new scale and draw a frame.**
!                                **************************************
      if(iplot.eq.2) then
         ymax= -1.e10
         ymin=  1.e10
!
         do i= 1,npt
         ymax= amax1(ymax,v(i))
         ymin= amin1(ymin,v(i))
         end do
!
         if(ymin.ge.0.) then
           ymax= 1.1*ymax
           ymin= 0.
         else
           ymax= amax1(0.,ymax)
           ymin= 1.1*ymin
         end if
      end if
!
      if(ymax.le.ymin) ymax= ymin+1.0
      if(iabs(il).eq.2) then
         if(ymax.gt.0.0) ymax= ymax+1.0
      end if
!
      dx= (xmax-xmin)/xcm(i1)
      dy= (ymax-ymin)/ycm(j1)
      x0= xmin
      y0= ymin
!
      call scalex (pl(i1),ql(j1),x0,y0,dx,dy,isc)
!
      pl1(isc)= pl(i1)
      pr1(isc)= pr(i1)
      ql1(isc)= ql(j1)
      qr1(isc)= qr(j1)
      xmin1(isc)= xmin
      xmax1(isc)= xmax
      ymax1(isc)= ymax
      ymin1(isc)= ymin
!                                                      *************
!                                                      **  frame. **
!                                                      *************
      call plot (pl(i1),ql(j1),3)
      call plot (pl(i1),qr(j1),2)
      call plot (pr(i1),qr(j1),2)
      call plot (pr(i1),ql(j1),2)
      call plot (pl(i1),ql(j1),2)
!                                                    ******************
!                                                    **  tick marks. **
!                                                    ******************
      scx= xcm(i1)/5.0
      scy= ycm(j1)/4.0
!
      x0= pl(i1)
      y1= ql(j1)
      y4= qr(j1)
      y2= y1 +0.25
      y3= y4 -0.25
!
      do 62 k=1,4
      x0= x0 +scx
      call plot (x0,y1,3)
      call plot (x0,y2,2)
      call plot (x0,y3,3)
      call plot (x0,y4,2)
   62 continue
!
      y0= ql(j1)
      x1= pl(i1)
      x4= pr(i1)
      x2= x1 +0.25
      x3= x4 -0.25
!
      do 63 k=1,3
      y0= y0 +scy
      call plot (x1,y0,3)
      call plot (x2,y0,2)
      call plot (x3,y0,3)
      call plot (x4,y0,2)
   63 continue
!                                                     **************
!                                                     ** numbers. **
!                                                     **************
!
      if(.not.lab_skip) then
        call number (pl(i1)-0.5,ql(j1)-0.45,hhs,xmin,0.,101)
        call number (pr(i1)-1.5,ql(j1)-0.45,hhs,xmax,0.,101)
!
        call number (pl(i1)-2.0,ql(j1)     ,hhs,ymin,0.,101)
        call number (pl(i1)-2.0,qr(j1)-0.30,hhs,ymax,0.,101)
      end if
!
!                                                     **************
!                                                     **  labels. **
!                                                     **************
      xc= 0.5*(pl(i1)+pr(i1))
      xu= xc -1.60
      xd= xc -0.20*n2/2
!
      yr= qr(j1)+0.15
      yl= ql(j1)-0.70
!
      call symbol (xu,yr,hh,lab1,0.,n1)
      call symbol (xd,yl,hh,lab2,0.,n2)
!
      xl= pl(i1)-1.50
      yc= 0.5*(ql(j1)+qr(j1))
      call symbol (xl,yc,hh,lab3,0.,n3)
!                                     **********************************
!                                     **  no plot is made if npt1 < 0 **
!                                     **********************************
   70 if(npt1.lt.0) return
!
      call plotl (u(1),v(1),isc,3)
!**
      if(iplot.eq.1) then
         do 100 i=1,npt
         call plotl (u(i),v(i),isc,2)
  100    continue
      else
         do 120 i=1,npt-1
         call plotl (u(i+1),v(i)  ,isc,2)
         call plotl (u(i+1),v(i+1),isc,2)
  120    continue
      end if
!**
      call plotl (u(npt),v(npt),isc,3)
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine scalex (xcm,ycm,x00,y00,dx,dy,isc)
!-----------------------------------------------------------------------
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      x0(isc)= x00
      y0(isc)= y00
      dxi(isc)= 1./dx
      dyi(isc)= 1./dy
!
      xl(isc)= xcm
      yl(isc)= ycm
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine plotl (x,y,isc,ipl)
!-----------------------------------------------------------------------
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      xcm= xl(isc) +dxi(isc)*(x -x0(isc))
      ycm= yl(isc) +dyi(isc)*(y -y0(isc))
!
      call plot (xcm,ycm,ipl)
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine values (x,y,height,val,theta,ifmat)
!-----------------------------------------------------------------------
!  << values >>
!     1. function
!        (1) to draw variable
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) x,y               (i)       absolute coordinate value
!        (2) height            (i)       draw out size on paper
!        (3) val               (i)       variable
!        (4) theta             (i)       angle
!        (5) ifmat             (i)       format type
!     3. called by
!             (** nothing **)
!     4. calls
!             (** number **)
!             (** symbol **)
!-----------------------------------------------------------------------
!        ifmat = (n100)*100 + keta
!        n100 = 0 : integer format
!        n100 = 1 : f format ::  number(x,y,height,val,theta,keta)
!        n100 = 2 : e format ::
!        n100 = 3 : power of ten format
!        n100 = othewise : not write out
!-----------------------------------------------------------------------
!
      real*4 val
      character chr13*13,chr12*12,chr3*3
      character*1 minus,zero,blank
      parameter(ratio = 6./7. )
      data minus/'-'/,zero/'0'/,blank/' '/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (ifmat.lt.0) return
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      n100 = ifmat/100
      keta = ifmat - n100*100
!
      if (n100.eq.0) then
        call number(x,y,height,val,theta,ifmat)
!       call number(x,y,height,val,theta,-1)
      else if (n100.eq.1) then
        call number(x,y,height,val,theta,ifmat)
!       call number(x,y,height,val,theta,keta)
      else if (n100.eq.2) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pe13.6)') val
          chr12(1:4) = chr13(1:3)//'e'
          numsym = 4
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else if (val.eq.0) then
            chrval = val
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+2) = chr13(2:keta+2)//'e'
            numsym = keta + 2
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or.   &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
        call symbol(999.,999.,height,chr3,theta,numsy1)
      else if (n100.eq.3) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pe13.6)') val
          chr12(1:6) = chr13(1:3)//'x10'
          numsym = 6
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+5) = chr13(1:keta+2)//'x10'
            numsym = keta + 5
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+4) = chr13(2:keta+2)//'x10'
            numsym = keta + 4
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        sint = sin(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
!
!                                             *******************
!                                             ** exponent part **
!                                             *******************
!
        h2 = height * 5./7.
        x1 = (numsym+1)* height * ratio
        y1 = height * 4./7.
        if (abs(theta).lt.1e-04) then
          x1 = x + x1
          y1 = y + y1
        else
          x2 =     x1 * cost - y1 * sint
          y1 = y + x1 * sint + y1 * cost + h2*cost
          x1 = x + x2                    - h2*sint
        end if
        call symbol(x1,y1,h2,chr3,theta,numsy1)
      end if
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine wdash (x1,y1,x2,y2,ipen )
!-----------------------------------------------------------------------
!  << wdash  >>                      ver 2.00   16.mar.1990
!
!     1. function
!        (1) to draw line from (x1,y1) to (x2,y2) by wdash
!                            in absolute coordinate
!     2. arguments            (i/o)     (meaning)
!        (1) x1,x2,y1,y2       (i)       absolute coordinate value
!        (2) ipen              (i)       pen type of 'wdash'
!     3. called by
!             (** eqcntr  **)
!             (** wdashl  **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
!       ipen : meaning           - : 0.05 (cm)
!        1   :       line     -------------------
!        2   :  dash line     --- --- --- --- ---
!        3   :  dash line     -- -- -- -- -- -- --
!        4   :  dash line     - - - - - - - - - -
!        5   :  1 point dash  ---- - ---- - ---- -
!        6   :  2 point dash  --2.0-- - - --2.0--
!   otherwise:  line          ---------------------
!-----------------------------------------------------------------------
!
      h1  =  0.05
      h2  =  2.0 * h1
      h3  =  3.0 * h1
      h4  =  4.0 * h1
      h20 = 20.0 * h1
      call plot ( x1 , y1 , 3 )
      k = - 1
      if(ipen.lt.2) then
        go to 999
      else if(ipen.eq.2) then
        hh1 = h3
        hh2 = h1
      else if (ipen.eq.3) then
        hh1 = h2
        hh2 = h1
      else if (ipen.eq.4) then
        hh1 = h1
        hh2 = h1
      else if (ipen.eq.5) then
        hh1 = h4
        hh2 = h1
        hh3 = h1
        hh4 = h1
      else if (ipen.eq.6) then
        hh1 = h20
        hh2 = h1
        hh3 = h1
        hh4 = h1
        hh5 = h1
        hh6 = h1
      end if
      if(ipen.lt.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hhh
  200   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hhh
          d=d+hh1
          goto 200
        end if
      else if (ipen.eq.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hhh
  500   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hhh
          d=d+hh1
          goto 500
        end if
      else if (ipen.eq.6) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hh5
        hh5 = hh6
        hh6 = hhh
  600   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hh5
          hh5 = hh6
          hh6 = hhh
          d=d+hh1
          goto 600
        end if
      end if
  999 call plot ( x2 , y2 , ( 5 + k ) / 2 )
      call plot ( x2 , y2 , 3)
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine daisho(x  ,nx,xmin1,xmax1)
!-----------------------------------------------------------------------
      dimension x(1)
!
      xmax1= x(1)
      xmin1= x(1)
      do 100 i=2,nx
      xmax1= amax1(xmax1,x(i) )
      xmin1= amin1(xmin1,x(i) )
  100 continue
      return
      end
!
!
!***************************************************************
!*     This program package generates a UNIX postscript        *
!*     graphic file when called by calcomp-compatible          *
!*     /plot23.f/.                                             *
!***************************************************************
!----------------------------------------------------------
!      PostScript header by fortran
!        T. Ogino (Nagoya University) February 27, 1992
!      Modified to conform GSIPP commands
!        Motohiko Tanaka (NIFS)       November 23, 1993
!
!----------------------------------------------- 5/27/96 -------
!     This PS-Adobe-2.0 header allows us full paging features in
!     the Ghostview.  To scroll up the page (backward), click the 
!     page number and press two buttons of mouse simultaneously.
!
!     Consult: A.Saitou (Kyoto U.)  The definition of /@eop  
!    needs stroke for line drawings (not in the TeX header).
!---------------------------------------------------------------
       subroutine gopen (nframe)
!----------------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
!*  This is an Adobe-2.0 postscript file.
!
       write(77,10)
   10  format('%!PS-Adobe-2.0',/       &
              '%%Pages: (atend)',/     &
              '%%PageOrder: Ascend',/  &
              '%%EndComments',/        &
              '%%BeginDocument')
!
!%%%%%%%%%%%%%%%%%%% Procedure Defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     write(77,11) 
!  11 format('%%BoundingBox: 150. 400. 550. 600.')
!
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/ &
             '/m {moveto} bind def  % x y m -- move to position')
!
      write(77,23) 
   23 format('/tr {/Times-Roman findfont} bind def',/ &
             '/sf {scalefont} bind def',/             &
             '/se {setfont} bind def',/               &
             '/ro {rotate}  bind def',/               & 
             '/tl {translate} bind def',/             &
             '/sc {scale} bind def')
!
      write(77,24) 
   24 format('/@bop          % @bop -- begin the a new page',/ &
             '{erasepage newpath initgraphics',/               &
             '/SaveImage save def',/                           &
             '} bind def')
!
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/ &
             '{stroke showpage',/                    &
             ' SaveImage restore',/                  &
             '} bind def')
!
      write(77,26) 
   26 format('/@end          % @end -- done the whole shebang',/ &
             ' /end load def')
!
      write(77,27) 
   27 format('/dir 0 def')
!
      write(77,29) 
   29 format('/s             % string s -- show the string',/ &
             '{dir 1 eq',/                                    &
             ' {gsave currentpoint translate 90 rotate 0 0 moveto',/ &
             ' show grestore}',/                              &
             ' {show} ifelse',/                               &
             '} bind def')
!
      write(77,31)
   31 format('%%EndDocument',/  &
             '%%EndProlog',/    &
             '%%BeginSetup',/   &
             '/Resolution 300 def',/ &
             '/#copies 1 def',/ &
             '%%EndSetup')
!
!%%%%%%%%%%%%%%%%%%% End of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*  initiate the page one.
!
       nfrm = nframe
!
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%Page:',1x,i2,1x,i2)
!
       write(77,30) 
   30  format('%%BeginPageSetup',/ &
              '%%EndPageSetup',/   &
              '@bop')
!
!
!*  Set magnifying factor (GSIPP to Sun coordinate).
!   Rotate and translate to output on A4-L paper.
!      Left corner ...... (  0.,  0.)
!      Right corner ..... (600.,780.)
!
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
!
       write(77,*) '90.0 ro'
       write(77,*) '50.0 -550.0 tl'
!
!*  If nfrm=4, four frames in a page (top-left frame).
!
       if(nfrm.eq.1) then
          write(77,*) '1.00 1.00 sc'
       else
          write(77,*) '0.50 0.50 sc'
          write(77,*) '0.0 550.0 tl'
       end if
!
       return
       end
!
!
!-----------------------------
       subroutine gclose
!-----------------------------
       call plote
       return
       end
!
!
!-----------------------------
       subroutine plote
!-----------------------------
       write(77,10) 
   10  format('@eop')
       return
       end
!
!
!-----------------------------------------
       subroutine chart
!-----------------------------------------
!*     Four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
!
!
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
!
!*  Frame 1: open a new page.
!
       if(loc.eq.0) then
          call plote
!
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
!
          write(77,10) 
   10     format('%%PageTrailer    % Need for the page count')
!
          write(77,20) lpage,lpage
   20     format('%%Page:',1x,i2,1x,i2)
!
          write(77,30) 
   30     format('%%BeginPageSetup',/ &
                 '%%EndPageSetup',/   &
                 '@bop')
!
          write(77,*) '90.0 ro'
          write(77,*) '50.0 -550.0 tl'
!
          if(nfrm.eq.1) then
             write(77,*) '1.00 1.00 sc'
          else
             write(77,*) '0.50 0.50 sc'
             write(77,*) '0.0  550.0 tl'
          end if
!
          return
       end if
!
!
!-----------------------------------------------------
!      First cancel the previous translation, then
!      make a new translation (scale factor alive).
!-----------------------------------------------------
!*   Frames 2-4:
!
       if(loc.eq.1) then
          write(77,*) '  0.0 -550.0 tl'
          write(77,*) '700.0  550.0 tl'
       end if
!
       if(loc.eq.2) then
          write(77,*) '-700.0 -550.0 tl'
          write(77,*) '   0.0    0.0 tl'
       end if
!
       if(loc.eq.3) then
          write(77,*) '  0.0 0.0 tl'
          write(77,*) '700.0 0.0 tl'
       end if
!
       return
       end
!
!
!------------------------------------
       subroutine factor(fct)
!------------------------------------
       write(77,10) fct,fct
   10  format(f6.2,1x,f6.2,' sc')
       return
       end
!
!
!------------------------------------
       subroutine newpen (ip)
!------------------------------------
       i1=(ip-1)/2
       i2=ip-2*i1
       write(77,*) 'sn'
       pi1=0.40*float(i1-1)
       write(77,30) pi1
   30  format(f3.1,' sl')
       if(i2.ne.1) then
       write(77,*) '[2 2] 0 sd'
       endif
       return
       end
!
!
!-----------------------------
       subroutine linee
!-----------------------------
       write(77,*) 'st'
       return
       end
!
!
!------------------------------------
       subroutine plot (x0,y0,ip)
!------------------------------------
!
       x= x0
       y= y0
       h= 0.
       n= 777
       call sunscl (x,y,h,n)
!
       if(ip.eq.3)  write(77,10) x,y
       if(ip.eq.2)  write(77,20) x,y
       if(ip.eq.-3) write(77,30) x,y
       if(ip.eq.-2) write(77,40) x,y,x,y
   10  format(f5.1,1x,f5.1,' m')
   20  format(f5.1,1x,f5.1,' l')
   30  format(f5.1,1x,f5.1,' tl')
   40  format(f5.1,1x,f5.1,' l sn',1x,f5.1,1x,f5.1,' tl')
!       write(77,*) 'st'
       return
       end
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       character    ica*80,ich(80)*1
       character(*) isymb   !!! NEC: must be necessary !
       equivalence (ica,ich(1))
!
       x= x0
       y= y0
       h= h0
       n= n0
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!*
       ica= isymb
       write(77,*) '(',(ich(i),i=1,n),') s'
!
       return
       end
!
!
!-----------------------------------------------
       subroutine number (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       character  isymb*9
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
!
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!
       if(abs(anu).gt.1.e1 .or.  &
          abs(anu).lt.1.e-1) then
        write(isymb,31) anu
   31   format(1pe9.2)
       else
        write(isymb,32) anu
   32   format(f7.2)
       end if
!
       if(.true.) go to 300
       if(abs(anu).lt.10000.) then  ! 5 digits
         if(abs(anu).gt.0.1) then
           write(isymb,40) anu
   40      format(f6.1)
         else
           if(abs(anu).gt.0.001) then  ! f6.3
             write(isymb,41) anu
   41        format(f6.3)
           else
             if(abs(anu).gt.0.001) then  ! f6.3
               write(isymb,42) anu   ! e9.2
   42          format(1pe9.2)
             else
               write(isymb,40) anu   ! 0.0
             end if
           end if
         end if
!
       else
         if(abs(anu).lt.100000.) then
           write(isymb,51) anu     ! f7.1
   51      format(f7.1)
         else
           write(isymb,52) anu     ! e9.2
   52      format(1pe9.2)
         end if
       end if
  300  continue
!
       write(77,*) '(',isymb,') s'
!
       return
       end
!
!
!-----------------------------------------------
       subroutine number2 (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       character  isymb*9
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
!
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!
      if(n0.eq.1) write(isymb,41) anu
      if(n0.eq.2) write(isymb,42) anu
   41  format(f6.1)
   42  format(f6.2)
!
       write(77,*) '(',isymb,') s'
!
       return
       end
!
!
!
!---------------------------------------------------
       subroutine sunscl (x,y,h,n)
!---------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
!
       if(x.eq.999.) then
         x= x0 +iabs(n0)*h0
       else
         x= fmag*x
         x0= x
       end if
!
       if(y.eq.999.) then
         y= y0
       else
         y= fmag*y
         y0= y
       end if
!
       h= fmag*h
       h0= h
       if(n.ne.777) n0= n
!
       return
       end
