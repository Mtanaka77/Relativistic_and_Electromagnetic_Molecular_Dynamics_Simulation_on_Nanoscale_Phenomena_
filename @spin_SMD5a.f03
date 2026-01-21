!*****************************************************************
!*                                                               *
!*   ### Dissipative Spin Molecular Dynamics Simulation ###      * 
!*                                                               *
!*      Numerical simulation code by:                            *
!*        Author: Motohiko Tanaka, PhD, Chikusa,Nagoya,Japan.    *
!*        https://github.com/Mtanaka77/                          *
!*                                                               *
!*        @spin_SMD5a.f03: Simulation code                       *
!*        param-spinRL5.h: Basic parameters                      *
!*        SAI105_config.START1: Configuration file               *
!*        magnetite8.xyz:  Initial xyz file                      * 
!*                                                               *
!*   Spin dynamics under the microwave H field for electrons,    *
!*   while nuclei feel forces from spins, Neighboring nuclei     * 
!*   and restoring forces on their lattice points.               *
!*                                                               *
!*   Ref: Selective heating mechanism of magnetite metal oxides  *
!*      by a microwave magnetic field, M. Tanaka, H. Kono, and   *
!*      K. Maruyama, Phys. Rev. B, 79, 104420 (2009).            *
!*                                                               *
!** First: 2007/7/07 *********************** Update: 2026/1/16 ***
!
!  Large number of cells... 
!  1) As convergence is poor, a small system is first organized. 
!   Then, a large system is equilibrated with the organized 
!   seed  cell. 
!                                                   11/07/2010
!
!  2) Computation time becomes huge for the MC step, thus the sum 
!   is divided on processors, and they are mpi_allreduced.
!                                                   11/12/2010
!
!  3) Array arguments to /forces/ are passed by common statements,
!   as execution fails on mol205-208 (DL385G2).     11/27/2010
!
!
!  On SR16000
!  1. Comment out - etime
!  2. In param.h for date and time, they are cdate*10, ctime*8 
!  3. File directory: /home/mtanaka/MPI_spin/... (char*29)
!  4. Use n_smp=1 - thread version does not work!  Comment out
!     call fftw_threads...
!  5. Use the script "f90sr3" with my /fftw3/lib
!
!  A) On SMP machines, call fftw_plan only on the major node,
!     Different nodes have different plan numbers.
!
!*****************************************************************
!*                                                               *
!*   Step 1: MC for equilibration at given T and B= 0.           *
!*   Step 2: MD (spin and nuclei dynamics)                       *
!*                                                               *
!*   In kstart= 0 of MC procedures, it is followed by MD steps.  *
!*      kstart > 0 using FT12 data                               *
!*---------------------------------------------------------------*
!*  >>NOTE for Fe and O:                                         *
!*                                                               *
!*       1  2 ... np1 ... np= np1+np2                            *
!*    1  Fe Fe Fe Fe  O   O                                      *
!*    2  Fe       Fe  O   O     rcut .... spins                  *
!*  ...  Fe       Fe  O   O     rcutC ... Coulomb                *
!*  np1  Fe Fe Fe Fe  O   O                                      *
!*  ...  O  O     O   O   O                                      *
!*   np  O  O     O   O   O                                      *
!*                                                               *
!*  Distances                                                    *
!*   Fe(A)-O  1.891A, 3.495A                                     *
!*   Fe(B)-O  2.058A, 3.563A, 3.661A                             *
!*                                                               *
!*   Fe(A)-Fe(A)  3.636 Ang                                      *
!*   Fe(B)-Fe(B)  2.969 Ang                                      *
!*   Fe(A)-Fe(B)  3.481 Ang                                      *
!*                                                               *
!*   O-O   2.85 Ang, 2.97 Ang, 3.088 Ang                         *  
!*                                                               *
!*    ... Definition of i1(),i2() in spin_dynamics               *           
!*    ... multi write(11) must be avoided                        *
!*                                                               *
!* 1. Potential field derived to fit ab initio is used, in which *
!*   short-range force is A*exp(-Br) -C/r**6 +D/r**12            *  
!* 2. Fe(3+) and Fe(2+) are balanced (2+2) on each sub-lattice.  *
!*                                                               *
!************************************************* 2009/02/11 ****
!*                                                               *
!*  a) Give three axes of the crystal in /init/                  *
!*  b) If isolated crystal - do not fold in /spin_dynamics/      *
!*                                                               *
!*    Position of atoms: x,y,z (Ang)                             *
!*    Spin of atoms: spx,spy,spz (multiple of 1 and 2)           *
!*    Site index: site= 1,2 (prime number)                       *
!*     Positions of Fe(2+) andFe(3+) at B-site are shuffled.     *
!*                                                               *
!*    .............................                              *
!*               site=1   site=2                                 *
!*    .............................                              *
!*      spec=1   (A)3+    (B)3+   3+                             *
!*      spec=2            (B)2+   2+                             *
!*      spec=0             O 2-                                  *
!*    .............................                              *
!*                                                               *
!*   1. Energy minimization                                      *
!*     Apply the Metropolis criterion when Usys increases        *
!*     U= - Sum J_ij s_i*s_j +Sum_i g*mu_b B*s_i                 *
!*                                                               *
!*   2. Spin dynamics                                            *
!*     - should start from an energy-minimized state             *
!*                                                               *
!*    ds_i          2*Jij             g*mue_B                    *
!*   ------ = Sum_j ----- s_i x s_j + ------- s_i x Bw )         *
!*     dt            hbar              hbar                      *
!*                                                               *
!*   Original: s_i: (hbar/2), J_ij: cm^2 gram/sec^3              *
!*   For exchange interactions Jaa, Jbb and Jab, first and       *
!*   second neighbors ( r < 4 A) are non-zero.                   *
!*                                                               *
!*   Units:                                                      *  
!*    t_unit= 1.d-12         ! ps                                *
!*    a_unit= 1.d-8          ! Ang                               *
!*    e_unit= 4.80325d-10    ! e                                 * 
!*    m_unit= 1.67261d-24    ! m_H                               *
!*                                                               *
!*****************************************************************
!* >> To get a free format of f90/f03, convert f77 into:         *
!*                                                               *
!*   Inside the vi editor: :%s/^c/!/, :%s/^*/!/, and :wq         *
!*   Outside the vi editor:                                      *
!*     tr 'a-z' 'a-z' <@spin_SMD5a.f >@spin_SMD5a.f03            *
!*                                                               *
!*   On Fortran 2003, "use, intrinsic :: iso_c_binding"          *
!*   Also,format statement on Fortran 2003 is DIFFERENT.         *
!*****************************************************************
!* >> Parallel Fortran 2003/2008                                 *
!*                                                               *
!* % mpif90 -mcmodel=medium -fpic -o ax.out @spin_SMD5a.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3 &> log
!* % mpiexec -n 5 ax.out &                                       *
!*                                                               *
!*  Debian-13 OS may need "-fallow-argument-mismatch"            *
!*****************************************************************
!
      program spin37
      use, intrinsic :: iso_c_binding
!
      include    'mpif.h'
      include    'param-spinRL5.h'
!
      integer(C_INT) size,rank,key,ierror,kstart, &
                     mx,my,mz,ir0,igrp,n_MCsteps,nt_p3m
      real(C_DOUBLE) cputime,wtime,wtime1,wtime2
      character*8    label,cdate*10,ctime*8
!
      common/parm4/  mx,my,mz
      common/parm5/  n_MCsteps,nt_p3m
      common/headr1/ label,cdate
      common/ranfff/ ir0  !<- ir0 rotate by ranff(x)
!     namelist/inp1/ kstart,praefix,suffix2*2,suffix1*2,suffix0*1,ir0
!
!  kstart = 0 ... New run with MC (t=0), followed by MD
!         = 1 ... Restart
!  
!  -------------------------------------------------
      kstart =  0   ! 1
      praefix= '105'
      suffix2= '1a' ! '1b' 
      suffix1= '1a' ! '1a' 
      suffix0= '1'  ! '1' 
!
      label  = 'sp3_nucl'
      nt_p3m =  5   ! in every 5 steps
!  -------------------------------------------------
      mx= meshx     !<<- param-spinRL5.h
      my= meshy
      mz= meshz
!  -------------------------------------------------
!
      call mpi_init (ierror)
      call mpi_comm_rank (mpi_comm_world,rank,ierror)
      call mpi_comm_size (mpi_comm_world,size,ierror)
!
!*******************************************
!*  A group for intra-task communication.  *
!*******************************************
! rank= 0... basic I/O
!
!     igrp = 1     ! group #
!     key= rank
!     call mpi_comm_split (mpi_comm_world,igrp,key,comm_local,ierror)
!
      praefixs = '/home/mtanaka/MPI_spin/SAI'//praefix 
      praefixi = '/home/mtanaka/MPI_spin/sai'//praefix 
      praefixc = '/home/mtanaka/MPI_spin/sai'//praefix 
!
! -------------------------------------------------------------
!* Input by keyboard or input file <aaa
!   Initial FT11 file=praefixc, form='formatted'
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2,form='formatted')
!
        write(11,*) 'This run uses',size,' processors...'
!       write(11,*) '         with',nthreads,' threads.'
!
        write(11,*) 'type kstart (0,1) ...'
!       read(05,'(i1)') kstart 
!
        write(11,*) 'Give the file name by praefix and suffix2'
        write(11,*) 'type praefix (a3)...'
!       read(05,'(a3)') praefix
!
        write(11,*) 'type suffix2 (a2)...'
!       read(05,'(a1)') suffix2
!
        write(11,*) 'type ir0 (i6)...'
!       read(05,'(i5)') ir0
!
        write(11,*) 'kstart =',kstart
        write(11,*) 'praefix=',praefix
        write(11,*) 'suffix2=',suffix2
        write(11,*) 'suffix1=',suffix1
        write(11,*) 'suffix0=',suffix0
        write(11,*) ' random seed: ir0=',ir0
!
        close (11)
      end if
!
!     call mpi_bcast ( kstart,1,mpi_integer,  0,mpi_comm_world,ierror)
!     call mpi_bcast (praefix,3,mpi_character,0,mpi_comm_world,ierror)
!     call mpi_bcast (suffix2,2,mpi_character,0,mpi_comm_world,ierror)
!     call mpi_bcast (suffix1,2,mpi_character,0,mpi_comm_world,ierror)
!     call mpi_bcast (suffix0,2,mpi_character,0,mpi_comm_world,ierror)
!     call mpi_bcast ( ir0,1,mpi_integer,     0,mpi_comm_world,ierror)
!                                             + ++++++++++++++
! -------------------------------------------------------------
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &              
                  status='unknown',position='append',form='formatted')
!
        call date_and_time7 (cdate,ctime)
!
        write(11,*) ' This run uses ',size,' processors...'
        write(11,*) '         '
        write(11,*) '   p3m is called in every',nt_p3m,' steps'
        write(11,*) '         '
        write(11,*) '   kstart =',kstart
        write(11,*) '         '
        write(11,*) '   today = ',cdate
        write(11,*) '   time  = ',ctime
!
        close (11)
!
!* Initial FT77 file=praefixc, form='formatted'
        open (unit=77,file=praefixc//'.77'//suffix2,form='formatted')
!
        nframe= 4
        call gopen (nframe)
        close (77)
      end if
!
      call clocks (cputime,wtime1)
! ------------------------------------------------------------------
      call spin_dynamics (rank,size,igrp,kstart)
! ------------------------------------------------------------------
      call clocks (cputime,wtime2)
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &              
                  status='unknown',position='append',form='formatted')
!
        wtime = wtime2 -wtime1
        write(11,'(" Final cpu, wtime (sec)=",2f10.3)') cputime,wtime
!
        close (11)
      end if
!*
  900 continue
!     call mpi_comm_free (comm_local,ierror)
      call mpi_finalize  (ierror)
!
      stop
      end program spin37
!
!
!------------------------------------------------------------------------
      subroutine date_and_time7 (date_now,time_now)
!------------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer, dimension(8) :: ipresent_time
      character(len=10) :: date_now
      character(len=8)  :: time_now

      call date_and_time (values=ipresent_time)

      write(time_now,'(i2,":",i2,":",i2)') ipresent_time(5:7)
      write(date_now,'(i4,"/",i2,"/",i2)') &
               ipresent_time(1),ipresent_time(2),ipresent_time(3)
!
      return
      end subroutine date_and_time7
!
!
!------------------------------------------------------------------------
      subroutine spin_dynamics (rank,size,igrp,kstart)
!------------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'mpif.h'
      include    'param-spinRL5.h'
!
      integer(C_INT) rank,size,igrp,kstart,np1,np2,            &
                 i1(0:num_proc),i2(0:num_proc),i3(0:num_proc), &
                 i4(0:num_proc),                               &
                 cnt_recv(0:num_proc),disp_recv(0:num_proc),   &
                 cnt_recvC(0:num_proc),disp_recvC(0:num_proc), &
                 cnt_send,ierror,mpierror,i_MC,ifcmp,kk
!
      real(C_DOUBLE) x,y,z,spx,spy,spz,ch,fx,fy,fz,fxC,fyC,fzC, &
                     Jint,r_ij,rintC,fc1,fc2,fcLJ
      common/partcl1/ &
                 x(np0),y(np0),z(np0),spx(np0),spy(np0),spz(np0), &
                 ch(np0)
      common/partcl2/ &
                 fx(np0),fy(np0),fz(np0),fxC(np0),fyC(np0), &
                 fzC(np0),Jint(nbx,np0),r_ij(nbx,np0),rintC(nbx,np0)
!
      real(C_DOUBLE) vx(np0),vy(np0),vz(np0),mass(np0),ag(np0),ep(np0),&
                 fkx(np0),fky(np0),fkz(np0),sp2(np0),spx0(np0),     &
                 spy0(np0),spz0(np0),spx1(np0),spy1(np0),spz1(np0), &
                 wsp1(3,np0),wsp2(3,np0),Jint_1,           &
                 u,u1,u2,spx00(np0),spy00(np0),spz00(np0), &
                 u_1,u_2,aspx(np0),aspy(np0),aspz(np0)
!
      integer(C_INT) nintS,lintS,nintC,lintC,if_LJ,spec,site
      common/partcl3/ & 
                 nintS(np0),lintS(nbx,np0),nintC(np0),lintC(nbx,np0), &
                 if_LJ(np0),spec(np0),site(np0)
!
      integer(C_INT) nlist(np0),lmax,k0,iac1,iac2,irej,modes,kwrite,  &
                 n_of_j,np10,np20,np000,np100,n_MCsteps,nt_p3m,       &
                 ia,ja,ka,i1x,i2x,i1y,i2y,i1z,i2z
!
      real(C_DOUBLE) t8,dt,dts,dt0,dth,Bex,Bey,Bez,spin2,spin3, &
                 Jaa,Jbb,Jab,J00,B00,Bapx,Bapy,Bapz,Bap,        &
                 tau_b,tau_diss,t_adv,fw,fw00,                  &
                 Temp,TCurie,g,mue_B,hbar,KJoule,Kcal,mol,kT,eV,&
                 omg_b,rnd,                                     &
                 f1,g1,b1,prb,qsx,qsy,qsz,rsx,rsy,rsz,hh1,hh2,  &
                 qq,rqq,dthg,dthqq,rdiv,spmin,                  &
                 alp,xx,yy,zz,rr,rr0,toler,deps,Usys8(nhs),     &
                 ss,smax,sav,t_unit,e_unit,a_unit,m_unit,m_Fe,m_O, &
                 pi,pi2,th,ph,tmax,tmax0,cptot,unif1(2),unif2(2),  &
                 xleng0,yleng0,zleng0,sss,rlist(np0),           &
                 buffer1(4),buffer2(4),del_en,wtime,            &
                 fc3,J_ki,wfdt,vth0,vth_O,vth_Fe,vth,           &
                 svx,svy,svz,sqrt2,vmax1,dgaus2,vsq1,vsq2,      &
                 vx1,vy1,vz1,mas,e_sp,e_c_r,e_LJ,e_Coulomb_p3m
!
      real(C_DOUBLE) x0(np0),y0(np0),z0(np0),vx0(np0),vy0(np0),vz0(np0)
      real(C_DOUBLE) x_0(np0),y_0(np0),z_0(np0)
      real(C_DOUBLE) dtwr,dtwr2,cputime
!
      real(C_DOUBLE) rad_Fe,rad_O,elj_Fe,elj_O,rcut,rcutC
      common/atoms/  rad_Fe,rad_O,elj_Fe,elj_O,rcut,rcutC
!
      integer(C_INT) it,is,iw,iwa,iwb,istop,                  &
                    iwrt1,iwrt2,iwrta,iwrtb,iter,itermax,     &
                    i00,i,j,k,l,nsite,notconv,nframe,mx,my,mz,&
                    ir0,lsite(num_proc,np0),ic,nstep_MC,nstep_MCt, &
                    ifdt,kmax,if_vres(np0),ifv
      common/parm1/ it,is
      common/parm2/ dtwr,dtwr2
      common/parm3/ t8,pi,dt,tmax,cptot
      common/parm4/ mx,my,mz
      common/parm5/ n_MCsteps,nt_p3m
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b, &
                    tau_diss,Temp,TCurie
      common/imemo/ iwa,iwb
      common/ranfff/ ir0
      common/itera/  toler,itermax
!
      real(C_float) spinx,spinz,spin7,Bextx,Bextz,magx,magy,magz, &
                    Usys,conv,aitr,psdt,tfix,uss,usb,tsx,tsy,tsz,sum_mb,&
                    U_Fe,U_O,ds_Fe,ds_O,fdt4,vdt4,idt4,timeh,dtrhs
      real(C_float) sx1(3),sy1(3),sz1(3),sx2(3),sy2(3),sz2(3),sn1(3), &
                    csx,csy,csz,axis(100),freq(100),t4,Bex4,Bez4, &
                    spx4(np0),spy4(np0),spz4(np0),    &
                    ch4(np0),x4(np0),y4(np0),z4(np0), &
                    vx4(np0),vy4(np0),vz4(np0)
      real(C_DOUBLE) wx,wy,wz,wn,wsx,wsy,wsz,wp, &
                    wx1,wy1,wz1,wn1,uav,wt1,wx7,wy7,wz7,wn7,uav7,wt7,  &
                    ssx,ssy,ssz,tht,phi,psi,bbw,tsz0,atsz,av_tsz(nhs), &
                    ub,um,ub1,um1,uu1(3),uu2(3),ranff,fdt8,vdt8,ds1,ds2
      common/ehist/ spinx(nhs),spinz(nhs),spin7(nhs),                   &
                    Bextx(nhs),Bextz(nhs),magx(nhs),magy(nhs),magz(nhs),&
                    Usys(nhs),conv(nhs),aitr(nhs),psdt(nhs),tfix(nhs),  &
                    uss(nhs),usb(nhs),tsx(nhs,3),tsy(nhs,3),tsz(nhs,3), &
                    sum_mb(nhs),U_Fe(nhs),U_O(nhs),ds_Fe(nhs),ds_O(nhs),&
                    fdt4(nhs),vdt4(nhs),idt4(nhs),timeh(nhs)
      logical       MC_first,ft06_start
      data          MC_first/.true./,ft06_start/.true./
!
      real(C_DOUBLE) alpha,xleng,yleng,zleng
      integer(C_INT) PP
      common/ewald1/ alpha,fc2
      common/ewald2/ PP
      common/ewald3/ xleng,yleng,zleng
!
      character*8    label,cdate*10,ctime*8,plot_ch*8,char*2
      real(C_float)  time
      common/headr1/  label,cdate
      common/headr2/  time
!
      istop = 0
!
      if(rank.eq.0) then 
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,*)
        if(kstart.eq.0) then
          write(11,'("## MC step with b_mw=0, followed by MD.....",/)') 
        else
          write(11,'("## MD step by restart.....",/)') 
        end if
!
        close (11)
      end if
!      
!--------------------------
!*  Cold start: MC and MD
!--------------------------
      if(kstart.eq.0) then
        t8= 0.d0
!
        it= 0
        is= 0   ! write out history
!
        iwa= -1
        iwb= -1
!
! Parameters are read
!
        call READ_CONF (xleng0,yleng0,zleng0,rank)
!
! xleng is expanded and redeined in "init"
!
        call init (x,y,z,ch,spx,spy,spz,sp2,spec,site,     &
                   xleng0,yleng0,zleng0,np1,np2,np10,np20, &
                   np000,i1x,i2x,i1y,i2y,i1z,i2z,rank,size)
!
        do i= 1,np1+np2
        x_0(i)= x(i)
        y_0(i)= y(i)
        z_0(i)= z(i)
        end do
!
        if(np1+np2.gt.np0) then
          write(06,*) ' # Stop: np > np0...',np1,np2,np0
          return
        end if
!
!--------------------------
!*  Restart from FT12
!--------------------------
      else
        if(rank.eq.0) then
          open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
          write(11,*) '# Read: ',praefixi//'.12'//suffix1
          close (11)
        end if
!
!
        open (unit=12,file=praefixi//'.12'//suffix1, & ! Only suffix1
                                   status='old',form='unformatted')
!
!                   +++++ +++++ +++++ mx*xleng0
        read(12) t8,xleng,yleng,zleng,rcut,rcutC,Temp,TCurie, &
                 tmax,dt,cptot
!                aaaa aa aaaaa        aaaa aaaaa aaaa aaaaaa    
        read(12) x,y,z,vx,vy,vz,ch,mass,ag,ep
        read(12) x_0,y_0,z_0,rintC
        read(12) spx,spy,spz,sp2,spx00,spy00,spz00,r_ij
        read(12) aspx,aspy,aspz
        read(12) Jint,u,spin2,spin3,n_MCsteps,tau_diss
!                                   ccccccccc cccccccc
        read(12) Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,toler,itermax
!                aaa aaa aaa aaaa aaaa aaaa aaaaa aaaaa aaaaaaa
        read(12) dtwr,dtwr2,fw00,atsz,tsz0,av_tsz
!                aaaa aaaaa 
        read(12) it,is,iwa,iwb,ir0
        read(12) np1,np2,spec,site,nintS,lintS,nintC,lintC,if_LJ
        read(12) i1,i2,i3,i4,disp_recv,cnt_recv,disp_recvC,cnt_recvC
        read(12) spinx,spinz,spin7,Bextx,Bextz,magx,magy,magz,       &
                 Usys,conv,aitr,psdt,tfix,uss,usb,tsx,tsy,tsz,sum_mb,&
                 U_Fe,U_O,fdt4,vdt4,idt4,timeh
        close(12)
!
!    # Change READ_CONF #
!          tmax, cptot, dtwr, dtwr2 
!
        call READ_CONF (xleng0,yleng0,zleng0,rank)
      end if
!
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        if(kstart.eq.1) write(11,*) '# kstart=1: Spin dynamics ........'
!
        write(11,*) '>> Input parameters...'
        write(11,'(" domains: mx, my, mz=",3i3)') mx,my,mz
        write(11,'(" kstart=",i3,/)') kstart
        write(11,933) Bapx,Bapy,Bapz,Temp, &
                      Jaa,Jbb,Jab,rcut,rcutC 
        write(11,'(" dt  =",1pd11.3,/," tolerance =",d11.3,/)') dt,toler
  933   format(' Applied H field(/100gauss) = ',3f8.1,/,&
               ' Temp(k) = ',0pf7.1,/,   &
               ' Jaa(meV) = ',1pd11.3,/, &
               ' Jbb(meV) = ',d11.3,/,   &
               ' Jab(meV) = ',d11.3,/,   &
               ' rcut(Ang) = ',d11.3,/,  &
               ' rcutC(Ang)= ',d11.3)
! 
        close (11)
      end if
!
!--------------------------
!* Define constants
!--------------------------
      pi = 4.d0*atan(1.d0)
      pi2= 2*pi
!
      t_unit= 1.d-12         ! ps
      a_unit= 1.d-8          ! Ang
      e_unit= 4.80325d-10    ! e
      m_unit= 1.67261d-24    ! m_H
!
      hbar = 6.62620d-27/(2.0d0*pi)  ! hbar, erg*sec
      mue_B= 9.27410d-21         ! e*hbar/2mc
      KJoule  = 1.d10            ! erg
      Kcal= 4.1868d0 *KJoule     ! 4.18 J/cal
      mol = 6.0220d23
!
      g  = 2.0023d0
      kT = 1.38062d-16*Temp
      eV = e_unit/300.d0
!
      m_Fe= 56.d0
      m_O = 16.d0
!
!      kT/m_unit= 2*1.38d-16*300/1.67d-24= 4.75d10
!      sqrt( kT/m_unit )= 2.22d5      
!        t_unit/a_unit= 1.d-4
      vth0= (t_unit/a_unit) * sqrt(2*kT/m_unit)  ! for hydrogen atom
      vth_O = vth0 /sqrt(m_O)
      vth_Fe= vth0 /sqrt(m_Fe)
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,'(" vth_O,vth_Fe=",1p2d11.3)') vth_O,vth_Fe
!
        close (11)
      end if
!
!* Exchange: J= kT_c/0.3z   cf. Kittel eq.(15-7)
!
!     znei= 6
!     asp = sp2(1)
!     J00= 3.0d0* (3/(2*znei*asp*(asp+1))) *1.38062d-16*TCurie
      J00 = 1.0d-2    ! 0.01 *1 eV = 10 meV
!       if(rank.eq.0) write(11,*) ' J00(eV)=',J00
!
      J00 = J00 * eV  ! unit in e_unit
      B00 = 100.d0    ! unit: 100 gauss
!
      f1 = J00 *t_unit/hbar
      g1 = g*mue_B*B00 *t_unit/hbar
!
      Bap= sqrt(Bapx**2 +Bapy**2 +Bapz**2)
      b1= g*mue_B*B00 *Bap  ! mue*B energy
!
!     tau_r= 0.d0
!     if(Bap.ne.0) tau_r= 1.d12* pi2/((g*mue_B*(B00*Bap))/hbar) ! in ps
      omg_b= pi2/tau_b
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,'(" TCurie=",f7.1,/," J00       =",1pd12.5)') TCurie,J00
        write(11,102) (J00*Jab)*spin2**2,b1*spin2,kT,omg_b,tau_b
  102   format(' Jab*s*s   = ',1pd12.5,/, &
               ' g*mub*s*b = ',  d12.5,/, &
               ' kT(erg)   = ',  d12.5,/, &
               ' omg_b (/ps), tau_b (ps)= ',2d12.5,/)
!
        dtrhs= f1*Jab*sp2(1)**2*dt
        if(abs(dtrhs).gt.0.4d0) then
          write(06,*) ' ## dt*rhs > 0.4d0 ..... stop'
          return
        end if
!
        write(11,103) dtrhs,tau_b,tau_diss
  103   format(' >> Increment of RHS for dt (negative Jab)=',f10.4, &
               /,'   tau_mw  = ',f10.2,  &
               /,'   tau_diss= ',f10.2,/)
!
        write(11,104) f1*Jab*sp2(1)**2,g1*Bap*sp2(1)
  104   format('## RHS of spin equation (per t_unit/hbar)...',/, &
               '   J*si*sj = ',f10.6,/, &
               ' g*mue*b*s = ',f10.6,/)
!
        write(11,105) vth_O
  105   format('## Thermal velocity of oxygen (a_unit/t_unit)...',/, &
               '   v_th(O)=',1pd12.3,/)
!
        close (11)
      end if
!
! ----------------------------
!*  Force constants
! ----------------------------
! ############################
!  unit: Coulomb force ->(e*tau/a)^2/m
!
      fc1  = J00*(t_unit/a_unit)**2/m_unit          ! spin force
      fc2  = (e_unit*t_unit)**2/(m_unit*a_unit**3)  ! Coulomb force
      fcLJ = 48.d0*t_unit**2*kT/(m_unit*a_unit**2)
!
      fc3  = m_unit*(a_unit/t_unit)**2    ! for kinetic energy
!
!* Order of j(r) as (1/r)**n
!
      n_of_j= 1
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,107) fc1,fc2,fc3,n_of_j
  107   format('## Force constants on nuclei ###',/, &
               '  fc_spin=',1pd12.4,/ &
               '  fc_Coul=',d12.4,/   &
               '  fc_kin =',d12.4,/   &
               '  n_of_j(r)=',i3,/)
!
        close (11)
      end if
! ############################
!
!  >> p3m_init ...
!* Define parameters of /ewald1-3/
!  -------------------------------
      pi = 4.d0*atan(1.d0)
!     alpha  = 2*pi/zleng   ! 0.26d0
      alpha  = 2*pi/rcutC   ! near if r < rcutC 
      PP     = p_max 
!  -------------------------------
!
!   fc2 -> prefactor= (t_unit*e_unit)**2/(w_unit*a_unit**3)
!     pref2    = 48.d0*t_unit**2*kbt/(m_unit*a_unit**2)
!
      call interpol_charge_assign_function (rank)
      call calculate_meshift (rank)
      call calculate_differential_operator (rank)
      call calculate_influence_function (rank)
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,601) rcut,rcutC,mx,my,mz 
  601   format('## p3m parameters ###',/, &
               '  rcut  =',f7.2,/ &
               '  rcutC =',f7.2,/ &
               '  [mx, my, mz] domains =',3i3,/)
!
        write(11,*) ' p3m successfully initialized !'
        close (11)
      end if
!
!--------------------------
!* Parallelization index
!--------------------------
!                  **
      do k= 0, size-1
      i1(k)=     k*(np1/size) +1
      i2(k)= (k+1)*(np1/size)
      if(k.eq.size-1) i2(k)= np1
!
      disp_recv(k)= 3*(i1(k) -1)
      cnt_recv(k) = 3*(i2(k) -i1(k) +1)
!
      if(kstart.eq.0 .and. rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
              status='unknown',position='append',form='formatted')
!
        if(k.eq.0) then
          write(11,*) '  '
          write(11,*) 'k, i1,i2, disp_recv, cnt_recv...'
        end if
        write(11,'("   k=",i3,4i8)') k,i1(k),i2(k),disp_recv(k),cnt_recv(k)
!
        close (11)
      end if
      end do
!
!                  **
      do k= 0, size-1
      i3(k)=     k*((np1+np2)/size) +1
      i4(k)= (k+1)*((np1+np2)/size)
      if(k.eq.size-1) i4(k)= np1 +np2
!
      disp_recvC(k)= 3*(i3(k) -1)
      cnt_recvC(k) = 3*(i4(k) -i3(k) +1)
!
      if(kstart.eq.0 .and. rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        if(k.eq.0) then
          write(11,*) '  '
          write(11,*) 'k, i3,i4, disp_recvC, cnt_recvC...'
        end if
        write(11,'("   k=",i3,4i8)') k,i3(k),i4(k),disp_recvC(k),cnt_recvC(k)
!
        close (11)
      end if
      end do
!
!------------------------
!* Interaction table (1)
!------------------------
! Pick up only nearest neighbors for A-A, B-B, A-B
! depending on my spin
!   site(i)= 1   ! A site Fe (8)
!   site(i)= 2   ! B site Fe (16)
!
      do 30 i= 1,np1
      l= 0
!
      do 32 j= 1,np1
      if(j.eq.i) go to 32   ! avoid j=i
!
      xx= x_0(i) -x_0(j)
      yy= y_0(i) -y_0(j)
      zz= z_0(i) -z_0(j)
!
!* Do not fold - if isolated crystal
!  #################################
!
      xx= xx -nint(xx/xleng)*xleng
      yy= yy -nint(yy/yleng)*yleng
      zz= zz -nint(zz/zleng)*zleng
      rr= sqrt(xx**2 +yy**2 +zz**2)
!
      if(rr.gt.rcut) go to 32
      l= l +1
      rlist(l)= rr
      nlist(l)= j
   32 continue
!
      lmax= l
!
      l= 0
   40 l= l +1
      if(l.gt.lmax) go to 30
!
      rr0= 1000
      if(l.gt.nbx) then
        write(06,*) 'l > nbx ... stop at 40 l=l+1...'
        return
      end if
!
! Sort the indices in ascending order
      do 50 k= 1,lmax
      if(rlist(k).le.rr0) then
        rr0= rlist(k)
        j  = nlist(k)
        k0 = k
      end if
   50 continue
!
      nintS(i)= l
      lintS(l,i)= j
      r_ij(l,i)= rr0
!
!  Exchange interaction is small for distant pairs: rr > 3 Ang
!
      nsite= site(i)*site(j)
      lsite(l,i)= nsite
      if(nsite.eq.1) Jint(l,i)= Jaa
      if(nsite.eq.2) Jint(l,i)= Jab
      if(nsite.eq.4) Jint(l,i)= Jbb
!
      rlist(k0)= 1000       ! Exclude in the next search
      go to 40
   30 continue
!
!------------------------
!* Interaction table (2)
!------------------------
!  Coulomb force
!
!    Fe(b)  2.058 x 6 - O  *
!           2.969 x 6 - Fe *
!           3.481 x 6 - Fe
!
!    Fe(a)  1.891 x 4  - O *
!           3.481 x 12 - Fe
!           3.495 x 12 - O
!
!      O    1.891 x 1 - Fe *
!           2.058 x 3 - Fe *
!           2.850 x 3 - O  *
!           2.970 x 6 - O  *
!           3.088 x 3 - O  *
!           3.495 x 3 - Fe    
!
!  Give the LJ-pair limit in sai100 data file
!    r_ij(t=0) < rcutC
!
!* Divided upon processors for a large system... (11/12/2010)
!
      do 60 i= i3(rank),i4(rank)
      l= 0
!
      do 62 j= 1,np1+np2
      if(j.eq.i) go to 62   ! avoid itself
!
      xx= x_0(i) -x_0(j)
      yy= y_0(i) -y_0(j)
      zz= z_0(i) -z_0(j)
!
      xx= xx -nint(xx/xleng)*xleng
      yy= yy -nint(yy/yleng)*yleng
      zz= zz -nint(zz/zleng)*zleng
      rr= sqrt(xx**2 +yy**2 +zz**2)
!
      if(rr.gt.rcutC) go to 62
      l= l +1
      rlist(l)= rr
      nlist(l)= j
   62 continue
!
      lmax= l
!
      l= 0
   70 l= l +1
      if(l.gt.lmax) go to 60
!
      if(l.gt.nbx) then
        write(06,*) 'l > nbx ... Stop at 70, nbx=',nbx
        return
      end if
      rr0= 1000
!
! Sort the indices in ascending order
      do 80 k= 1,lmax
      if(rlist(k).le.rr0) then
        rr0= rlist(k)
        j  = nlist(k)
        k0 = k
      end if
   80 continue
!
      nintC(i)= l
      lintC(l,i)= j
      rintC(l,i)= rr0
!
!  Use the same list, but limit the LJ pair 
!     for r_ij(0) < 3.7 ang (almost 2nd neighbors) -> 7.0 ang for 5x5x5 cells
!                                    9/18/2010                10/27/2010
      if_LJ(i)= 0
      if(rr0.lt.rcutC) if_LJ(i)= 1
!
      rlist(k0)= 1000       ! Exclude in the next search
      go to 70
   60 continue
!
!
!%%%% Only at t=0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
      if(kstart.eq.0) then
!
! ------------------------------------------------
!*  The mass, radius, LJ, and initial velocity
! ------------------------------------------------
!
      do i= 1,np1+np2
      if(i.le.np1) then
        mass(i)= m_Fe
        ag(i)= rad_Fe
        ep(i)= elj_Fe*(Kcal/mol)/kT
!
      else
        mass(i)= m_O
        ag(i)= rad_O
        ep(i)= elj_O*(Kcal/mol)/kT
      end if
      end do
!
!***************
!*  Velocity.  *
!***************
!
      sqrt2= dsqrt(2.d0)
      call ggauss 
!
!* For irons
!
      svx= 0.
      svy= 0.
      svz= 0.
!
      do i= 1,np1
      vmax1= vth0 /sqrt(mass(i)) ! sqrt2*vth0/sqrt(mass(i))
      vx(i)= dgaus2(vmax1)
      vy(i)= dgaus2(vmax1)
      vz(i)= dgaus2(vmax1)
!
      svx= svx +vx(i)
      svy= svy +vy(i)
      svz= svz +vz(i)
      end do
!
      do i= 1,np1
      vx(i)= vx(i) -svx/np1
      vy(i)= vy(i) -svy/np1
      vz(i)= vz(i) -svz/np1
      end do
!
!
      svx= 0.
      svy= 0.
      svz= 0.
!
      do i= np1+1,np1+np2
      vmax1= sqrt2*vth0/sqrt(mass(i))
      vx(i)= dgaus2(vmax1)
      vy(i)= dgaus2(vmax1)
      vz(i)= dgaus2(vmax1)
!
      svx= svx +vx(i)
      svy= svy +vy(i)
      svz= svz +vz(i)
      end do
!
      do i= np1+1,np1+np2
      vx(i)= vx(i) -svx/np2
      vy(i)= vy(i) -svy/np2
      vz(i)= vz(i) -svz/np2
      end do
      end if  
!
!%
!%%%% Only at t=0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------
      if(rank.eq.0) then
        open (unit=21,file=praefixc//'.21'//suffix2,form='formatted')
!
        write(21,*) '### Spin tables ###'
!
        do 901 i= 1,np1
        write(21,936) i,sp2(i),site(i),nintS(i)
  936   format('i=',i5,/,'  spin,site,nint=',f5.1,2i3)
!
        write(21,'(" lintS=",15i5)') (lintS(k,i),k=1,nintS(i))
        write(21,'(" lsite=",15i5)') (lsite(k,i),k=1,nintS(i))
        write(21,'(" r_ij=",10f6.2)') (r_ij(k,i),k=1,nintS(i))
  901   continue
!
        write(21,*) '  '
        write(21,*) '### coulomb force pairs...'
        do 910 i= 1,np1+np2
        if(i.eq.1) write(21,*)     '>> For iron ...'
        if(i.eq.np1+1) write(21,*) '>> For oxygen ...'
!
        if(i.le.np1) then
          if(site(i).eq.1) char= ' A'
          if(site(i).eq.2) char= ' B'
        else
          char= '  '
        end if
!
        write(21,940) i,nintC(i),char,site(i),spec(i)
        write(21,941) (k,lintC(k,i),rintC(k,i),k=1,nintC(i))
  940   format('i=',i5,' nintC=',i3,' site(a,b)=',a2,' site, spec=',2i2)
  941   format(' k,j,r_ij= ',i3,i5,f10.3)
  910   continue
        close (21)
!
        open (unit=13,file=praefixc//'.13'//suffix2,form='unformatted')
!
        write(13) xleng,yleng,zleng,rcut,rcutC,dt
        write(13) Jaa,Jbb,Jab,Bapx,Bapy,Bapz,Jint
        write(13) np1,np2,spec,site,nintS,lintS
        write(13) nintC,lintC,if_LJ
        close(13)
      end if
!------------------------------------------------------------------
!******************************************************************
!*  Time Loop or Minimization Loop begins                         *
!******************************************************************
!
      if(kstart.eq.1) del_en= 0
!
!                     !<<-- go to 1000: 250 lines below
 1000 continue
      dth= 0.5d0*dt
!
      it= it + 1
      t8= t8 + dt
!    ++++++++++++
!
!* All nodes must share the same info
!  time synchronization is achieVed here...
!
      call clocks (cputime,wtime)
      if(t8.gt.tmax) istop= 1
      if(wtime/60.d0.gt.cptot) istop= 1
!
      call mpi_bcast (istop,1,mpi_integer,0,mpi_comm_world,ierror)
      if(istop.eq.1) go to 7000
!
!     if(kstart.ne.0) then
!       dtwr =  tau_b/80.d0
!       dtwr2 = tau_b/4.d0
!     end if
!
      iwrt1= iwrta(t8,dtwr)
      iwrt2= iwrtb(t8,dtwr2)
      iter = 0
!
      if(it.eq.1) then
        wx1= 0
        wy1= 0
        wz1= 0
        wn1= 0
        uav= 0
        wt1= 0
      end if
!
!*****************************************************************
!* MC step is executed by all nodes if kstart= 0                 * 
!*  Arrays spx00...spz00() are used in MC - do not use spx-spz ! *
!*****************************************************************
!
 2100 continue
      if(kstart.eq.0) then
!     ++++++++++++++++++++
!
!  kstart= 0
        if(MC_first) then   ! Organize a seed one cell
          np100= np10       ! One cell - spx-spz
          nstep_MC = 50001
          kwrite   = 10000
!
        else                ! Initialize a large system with the organized seed cell
          np100= np1
          nstep_MC = n_MCsteps  ! 1000001 
          kwrite   =   50000    !   50000
        end if
!
        if(rank.eq.0) then
          open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
          if(MC_first) then
            write(11,*)
            write(11,*) 'Step 1: MC for a small (one cell) system...'
          else
            write(11,*)
            write(11,*) 'Step 2: Equilibration of a large system...'
            write(11,*) '     maximum iteration count=',n_MCsteps
          end if
!
          close (11)
        end if
!
        do i= 1,np1
        spx00(i)= spx(i)
        spy00(i)= spy(i)
        spz00(i)= spz(i)
        end do
!                     !<<-- 2000: 24 lines below
        go to 2000    !<<-- kstart= 0
      else
!
!  kstart= 1
!  +++++++++
        if(it.eq.1) then
          call magnetiz (spx,spy,spz,g,wx7,wy7,wz7,wn7,u1,uav7,wt7,np1,-1)
!                                              +++                    +-1
!* Sign of magnetic field is defined here
!
          if(wz7.gt.0) then
            fw00= 1.d0    ! Choose so that -m*b < 0
          else
            fw00= -1.d0
          end if
        end if
!
        go to 3000    !!<-- 3000: 211 lines below 
      end if
!
!**************************************
!* Energy minimization by MC process  *
!**************************************
!  Flip a spin
!
 2000 i_MC= 0  !<-- kstart= 0
!
      Bex= 0
      Bey= 0
      Bez= 0
!
      u= 1.d10
      iac1= 0
      iac2= 0
      irej= 0
!
 2300 i_MC= i_MC +1
      j= mod(int(np100*ranff(0.d0)+0.00001),np100) +1
!                  +++                      +++
      th=  pi*ranff(0.d0)
      ph= pi2*ranff(0.d0)
!
      spx0(j)= spx00(j)   ! Save
      spy0(j)= spy00(j)
      spz0(j)= spz00(j)
!
      spx00(j)= sp2(j)*sin(th)*cos(ph)
      spy00(j)= sp2(j)*sin(th)*sin(ph)
      spz00(j)= sp2(j)*cos(th)
!
      u1= 0.d0
!
!* Calculating pairs both ways - 0.5*
!
      if(MC_first) then
!
        do 200 i= 1,np100
        u1= u1 + b1 * (Bex*spx00(i) +Bey*spy00(i) +Bez*spz00(i))
!
        do 210 l= 1,np000   !<<-- np000 
        if(l.eq.i) go to 210
!
        xx= x(i)-x(l)
        yy= y(i)-y(l)
        zz= z(i)-z(l)
!                            **       + Small system
        xx= xx -nint(xx/xleng0)*xleng0
        yy= yy -nint(yy/yleng0)*yleng0
        zz= zz -nint(zz/zleng0)*zleng0
        rr= sqrt(xx**2 +yy**2 +zz**2)
!
        if(rr.gt.rcut) go to 210
        nsite= site(i)*site(l)
!
        if(nsite.eq.1) Jint_1= Jaa
        if(nsite.eq.2) Jint_1= Jab
        if(nsite.eq.4) Jint_1= Jbb
        u1= u1 - 0.5d0*J00* Jint_1*(spx00(i)*spx00(l) &
                              +spy00(i)*spy00(l) +spz00(i)*spz00(l))
  210   continue
  200   continue
!
      else
!
!* Divided upon processors for a large system... (11/12/2010)
!
        do 220 i= i1(rank),i2(rank)
        u1= u1 + b1 * (Bex*spx00(i) +Bey*spy00(i) +Bez*spz00(i))
!
        do 220 k= 1,nintS(i)
        l= lintS(k,i)
        u1= u1 - 0.5d0*J00* Jint(k,i)*(spx00(i)*spx00(l) &
                              +spy00(i)*spy00(l) +spz00(i)*spz00(l))
  220   continue
!
!
        u_1= u1
        call mpi_allreduce (u_1,u_2,1,mpi_real8,mpi_sum,mpi_comm_world, &
                            mpierror)
        u1= u_2
      end if
!
!*****************************************************
!* Apply the Metropolis criterion when u increases   *
!*****************************************************
!
      if(u1.lt.u) then            ! Accept the flip
        u= u1                     ! keep lowest energy
        iac1= iac1 +1
      else
!
        prb= exp(-(u1-u)/kT)      ! (u1-u) for one excitation
        if(prb.gt.ranff(0.d0)) then
          iac2= iac2 +1
          u= u1                   ! <-- keep this energy
!
        else
          irej= irej +1           ! Reject the trial flip and
          spx00(j)= spx0(j)       ! restore the spin
          spy00(j)= spy0(j)
          spz00(j)= spz0(j)
        end if
      end if
!                          

      if(rank.eq.0) then  !                                          Add
        call magnetiz (spx00,spy00,spz00,g,wx7,wy7,wz7,wn7,u1,uav7,wt7,np100,1)
!                                                                      +++++ 1
        if(mod(i_MC,kwrite).eq.1) then  !                              Average                
          call magnetiz (spx00,spy00,spz00,g,wx7,wy7,wz7,wn7,u1,uav7,wt7,np100,7)
!                                                                        +++++ 7
          bbw= B00*Bez
!
          if(rank.eq.0) then
            open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
            write(11,201) i_MC,iac1,iac2,irej,u/np100,u1/np100, &
                          bbw,wx7,wy7,wz7     !    +++      +++
  201       format(' i_MC,iac1,iac2,irej=',i7,3i7,' <u>,<u1>=', &
                    1p2d11.3,' bz=',d9.2,' <m>=',0p3f6.2)
            close (11)
          end if
!                                                                      Reset
          call magnetiz (spx00,spy00,spz00,g,wx7,wy7,wz7,wn7,u1,uav7,wt7,np100,0)
!                                                                        +++++ 0
        end if
      end if
!                <- 2300: 114 lines above          
      if(i_MC.lt.nstep_MC) go to 2300
!     -------------------------------
!
      if(MC_first) then
        MC_first= .false.
!      ++++++++++++++++++
!
        i= 0
        do 230 ia= i1x,i2x
        do 230 ja= i1y,i2y
        do 230 ka= i1z,i2z
!
! To enable convergence to equilibrium...  11/13/2010
!
        if(mod(ia,3).eq.0 .and. &
           mod(ja,3).eq.0 .and. &
           mod(ka,3).eq.0) then
!
! Seed cell
          do 240 l= 1,np000  !<<--- np000 
          i= i +1
          spx(i)= spx00(l)  ! spx() are used in MD
          spy(i)= spy00(l)
          spz(i)= spz00(l)
  240     continue
!
        else
!!        do 243 l= 1,np10
          do 243 l= 1,np000
          i= i +1
!
! To avoid always being s=0... 11/23/2010
!  Note: specifying spx(i)= sp2(i)*ranff(0) is not good
          th=  pi*ranff(0.d0)
          ph= pi2*ranff(0.d0)
          spx(i)= sp2(i)*sin(th)*cos(ph)
          spy(i)= sp2(i)*sin(th)*sin(ph)
          spz(i)= sp2(i)*cos(th)
  243     continue
        end if
  230   continue
!
        if(rank.eq.0) then
          open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
          write(11,*) 'np1, np000*cells=',np1,i
          close (11)
        end if
!
        go to 2100      !<-- 1000: 247 lines above
!
      else
        do 250 i= 1,np1
        spx(i)= spx00(i)  ! spx() are used in MD
        spy(i)= spy00(i)
        spz(i)= spz00(i)
!
        aspx(i)= spx00(i)
        aspy(i)= spy00(i)
        aspz(i)= spz00(i)
  250   continue
      end if
!
      kstart= 1
      it= 0   !<<-- 
      t8= 0
      go to 1000     !<-- 1000: 326 lines above
!
!-----------------------------------------------------------
!* Spin dynamics with dissipation
!    start from an equilirated state ...
!
!    ds_i                 2*jij       g*mue_B
!   ------ = s_i x (sum_j ----- s_j - ------- be)
!     dt                   hbar        hbar
!
!     f1= J00 *t_unit/hbar        J00= 10 meV
!     g1= g*mue_B*B00 *t_unit/hbar   B00= 100 gauss
!-----------------------------------------------------------
!              !<-- ksrart= 1
 3000 continue
!                                           *****
      fw= fw00*(1.d0 -exp(-t8/tau_b)) * sin(omg_b*t8)  ! -m*b < 0
      Bex= Bapx*fw  ! Bapx * sin(omg_b*t)
      Bey= Bapy*fw
      Bez= Bapz*fw
!
      iter= 0
  400 iter= iter +1
!
      do 301 i= 1,np1
      if(iter.eq.1) then
        spx1(i)= spx(i)
        spy1(i)= spy(i)
        spz1(i)= spz(i)
      end if
!
      spx0(i)= spx1(i)
      spy0(i)= spy1(i)
      spz0(i)= spz1(i)
  301 continue
!
      do 300 i= i1(rank),i2(rank)
      qsx= 0
      qsy= 0
      qsz= 0
!
      do 320 k= 1,nintS(i)
      j= lintS(k,i)
!
!* Fold back - corrected 11/23/2010
!
      xx= x(i) - x(j)
      yy= y(i) - y(j)
      zz= z(i) - z(j)
      xx= xx -nint(xx/xleng)*xleng
      yy= yy -nint(yy/yleng)*yleng
      zz= zz -nint(zz/zleng)*zleng
!
      rr = sqrt(xx**2 +yy**2 +zz**2)
      J_ki= Jint(k,i)*(r_ij(k,i)/rr)**n_of_j
!
      qsx= qsx +f1*J_ki*(spx1(j) +spx(j))
      qsy= qsy +f1*J_ki*(spy1(j) +spy(j))
      qsz= qsz +f1*J_ki*(spz1(j) +spz(j))
  320 continue !   f1*J_ki * (s_1 +s)
!
!*********************************************************
!* Equation: ds_i/dt= s_i*qq/hbar -(s_i-s_i0)/tau,       *
!*               where qq= sum_j (2*Jij*s_j) -g*mue_B*B  *
!*********************************************************
!
      qsx= qsx -g1*Bex
      qsy= qsy -g1*Bey
      qsz= qsz -g1*Bez
!
!* RHS
      hh2= dth/tau_diss
      rsx= spx(i) +dth*(spy(i)*qsz -spz(i)*qsy) 
      rsy= spy(i) +dth*(spz(i)*qsx -spx(i)*qsz) 
      rsz= spz(i) +dth*(spx(i)*qsy -spy(i)*qsx) &
                            -hh2*(spz(i) -2.d0*aspz(i))
!
! RHS_para
      qq  = qsx**2 +qsy**2 +qsz**2
      rqq = (rsx*qsx +rsy*qsy +rsz*qsz)/qq
!
      dthqq= dth**2 * qq
      rdiv= 1.d0/(1.d0 +dthqq)
!
      spx1(i)= (rsx +dth*(rsy*qsz -rsz*qsy) + dthqq*rqq*qsx)*rdiv
      spy1(i)= (rsy +dth*(rsz*qsx -rsx*qsz) + dthqq*rqq*qsy)*rdiv
      spz1(i)= (rsz +dth*(rsx*qsy -rsy*qsx) + dthqq*rqq*qsz)*rdiv
!
!-------------------------------
!* Mix value of two iterations
!-------------------------------
!
      alp= 0.7d0
      spx1(i)= alp*spx1(i) +(1.d0 -alp)*spx0(i)
      spy1(i)= alp*spy1(i) +(1.d0 -alp)*spy0(i)
      spz1(i)= alp*spz1(i) +(1.d0 -alp)*spz0(i)
  300 continue
!
!----------------------
!* Unify the data
!----------------------
      i00= i1(rank) -1
!
      do 360 i= i1(rank),i2(rank)
      wsp1(1,i-i00)= spx1(i) 
      wsp1(2,i-i00)= spy1(i) 
      wsp1(3,i-i00)= spz1(i) 
  360 continue
!
      cnt_send= 3*(i2(rank) -i1(rank) +1)
      call mpi_allgatherv (wsp1,cnt_send,          mpi_real8, &
                           wsp2,cnt_recv,disp_recv,mpi_real8, &
                           mpi_comm_world,ierror)
!
      do 370 i= 1,np1
      spx1(i)= wsp2(1,i)
      spy1(i)= wsp2(2,i)
      spz1(i)= wsp2(3,i)
  370 continue
!
!* Tolerance= 1.d-10 and dt= 2.d-3 (ps)
!   -> iter= 6, error of 0.29% at t= 1000 ps
!
      toler= 1.d-10
      smax= 0
      sav= 0
      notconv= 0
!
      do 390 i= 1,np1
      ss= ((spx1(i) -spx0(i))**2 +(spy1(i) -spy0(i))**2  &
                                 +(spz1(i) -spz0(i))**2) &
                 /(spx0(i)**2 +spy0(i)**2 +spz0(i)**2)
      smax= dmax1(ss,smax)
      sav = sav +ss
!
      if(abs(ss).gt.toler) then
        notconv= notconv +1
      end if
  390 continue
!                                 <- Return to 116 lines above
      if(notconv.ne.0 .and. iter.lt.itermax) go to 400  !!<<-- 
!     ------------------------------------------------
!
!* Decrease in spz(i) should be converted to s_perp(i), as sp2(i)= const.
!
      do 500 i= 1,np1
      ss= sp2(i)**2 -spz1(i)**2  ! Could be negative by numerics
!
      spmin= 0.05d0*sp2(i)
      if(ss.gt.spmin**2) then
        qq= sqrt(ss/(spx1(i)**2 +spy1(i)**2))
      else
        qq= spmin   ! phase information must be retained
      end if
!
      spx(i)= qq*spx1(i)
      spy(i)= qq*spy1(i)
      spz(i)= spz1(i)
  500 continue
!
! ---------------------------------------------------------------
!    Mass unit: m_unit is included in fc1, fc2
! ---------------------------------------------------------------
!* Use a sub-time mesh, such that  v*dts < 1
!
      ifdt= 0
!
      do i= 1,np1+np2
      x0(i)= x(i)     !<- x0(i) these 43 lines below, x(i)
      y0(i)= y(i)
      z0(i)= z(i)
!
      vx0(i)= vx(i)   !<-- vx(i)
      vy0(i)= vy(i)
      vz0(i)= vz(i)
      end do
!
!  p3m - Slow change 
!
      if(mod(it,nt_p3m).eq.1) then
        call p3m_perform (x,y,z,ch,fkx,fky,fkz,np1+np2,iwrt1, &
                          rank,e_Coulomb_p3m)
      end if
!
!
 4000 if(ifdt.eq.0) then
        kmax= 1
        dts= dt
!
      else
        kmax= 2*kmax  !<-- Small steps
        dts= dt/kmax
!
        do i= 1,np1+np2
        x(i)= x0(i)   !<- x0(i) 43 lines above
        y(i)= y0(i)
        z(i)= z0(i)
!
        vx(i)= vx0(i)
        vy(i)= vy0(i)
        vz(i)= vz0(i)
        end do
      end if
!
!
      do 5000 kk= 1,kmax
      call forces (fc1,fc2,fcLJ,alpha,e_sp,e_c_r,e_LJ, &
                   i1,i2,i3,i4,rank,np1,np2,it)
!
!----------------------
!* Unify the forces
!----------------------
! 1) Spin force
      i00= i1(rank) -1
!
      do i= i1(rank),i2(rank)
      wsp1(1,i-i00)= fx(i) 
      wsp1(2,i-i00)= fy(i) 
      wsp1(3,i-i00)= fz(i) 
      end do
!
      cnt_send= 3*(i2(rank) -i1(rank) +1)
      call mpi_allgatherv (wsp1,cnt_send,          mpi_real8, &
                           wsp2,cnt_recv,disp_recv,mpi_real8, &
                           mpi_comm_world,ierror)
!
      do i= 1,np1
      fx(i)= wsp2(1,i)
      fy(i)= wsp2(2,i)
      fz(i)= wsp2(3,i)
      end do
!
! 2) Coulomb force
!
      i00= i3(rank) -1
!
      do i= i3(rank),i4(rank)
      wsp1(1,i-i00)= fxC(i) 
      wsp1(2,i-i00)= fyC(i) 
      wsp1(3,i-i00)= fzC(i) 
      end do
!
      cnt_send= 3*(i4(rank) -i3(rank) +1)
      call mpi_allgatherv (wsp1,cnt_send,            mpi_real8, &
                           wsp2,cnt_recvC,disp_recvC,mpi_real8, &
                           mpi_comm_world,ierror)
!
      do i= 1,np1+np2
      fxC(i)= wsp2(1,i)
      fyC(i)= wsp2(2,i)
      fzC(i)= wsp2(3,i)
      end do
!
! 3) Energies
!
      uu1(1)= e_sp
      uu1(2)= e_c_r
      uu1(3)= e_LJ
      call mpi_allreduce (uu1,uu2,3,mpi_real8,mpi_sum,mpi_comm_world, &
                          mpierror)
!
      e_sp = uu2(1)
      e_c_r= uu2(2)
      e_LJ = uu2(3)
!
!
!* Full forces should be used from t=0 to avoid atom thermal diffusion
!
      do i= 1,np1+np2
      vx(i)= vx(i) +dts*(fx(i) +fxC(i) +fkx(i))/mass(i)
      vy(i)= vy(i) +dts*(fy(i) +fyC(i) +fky(i))/mass(i)
      vz(i)= vz(i) +dts*(fz(i) +fzC(i) +fkz(i))/mass(i)
!
      x(i)= x(i) +dts*vx(i) !<- x(i)
      y(i)= y(i) +dts*vy(i)
      z(i)= z(i) +dts*vz(i)
      end do
!
!  Check all atoms separately
!
      fdt8= 0
      vdt8= 0
      ifv= 0
!
      do i= 1,np1+np2
      fdt8= dts* dmax1(abs(fx(i)+fxC(i)+fkx(i)), &
                       abs(fy(i)+fyC(i)+fky(i)), &
                       abs(fz(i)+fzC(i)+fkz(i)) )/mass(i)
      vdt8= dts* dmax1(abs(vx(i)),abs(vy(i)),abs(vz(i)))
!
      vth= vth0/sqrt(mass(i))
      if(fdt8.gt.0.2*vth .or. vdt8.gt.0.3d0) then
!
        if(.false.) then
!       if(ifdt.ge.5 .and. rank.eq.0) then
          open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
           write(11,570) it,t8,i,ifdt,dts,fdt8/vth,vdt8
           write(11,571) x(i),y(i),z(i),vx(i),vy(i),vz(i), &
                          dts*(fx(i)+fxC(i)+fkx(i))/(mass(i)*vth), &
                          dts*(fy(i)+fyC(i)+fky(i))/(mass(i)*vth), &
                          dts*(fz(i)+fzC(i)+fkz(i))/(mass(i)*vth)
  570      format('it,t8=',i8,f8.1,'  i,ifdt=',i4,i2, &
                  '   dts,fdt8/vth,vdt8=',1pd10.2,0p2f10.5)
  571      format(5x,'x:',1p3d11.3,3x,'vx:',3d11.3,/,5x,'fdt:',3d11.3)
!
           close (11)
        end if
!
        ifv= ifv +1
        if(ifv.ge.1) then
          ifdt= ifdt +1
          go to 4000
        end if
      end if
      end do
 5000 continue 
!
! -----------------------------------------------
!* Calculate the system energy (in all steps)
! -----------------------------------------------
!
 5300 continue 
      ub1= 0
      um1= 0
!
      do 590 i= i1(rank),i2(rank)
      ub1= ub1 + b1 * (Bex*spx(i) +Bey*spy(i) +Bez*spz(i))
!
      do 590 k= 1,nintS(i)
      j= lintS(k,i)
!
!* Fold back
      xx= x(i) - x(j)
      yy= y(i) - y(j)
      zz= z(i) - z(j)
      xx= xx -nint(xx/xleng)*xleng
      yy= yy -nint(yy/yleng)*yleng
      zz= zz -nint(zz/zleng)*zleng
!
      rr = sqrt(xx**2 +yy**2 +zz**2)
      J_ki= Jint(k,i)*(r_ij(k,i)/rr)**n_of_j
!
      um1= um1 - 0.5d0*J00* J_ki*(spx(i)*spx(j) +spy(i)*spy(j) +spz(i)*spz(j))
  590 continue
!
      uu1(1)= ub1
      uu1(2)= um1
      uu1(3)= 0
      call mpi_allreduce (uu1,uu2,3,mpi_real8,mpi_sum,mpi_comm_world, &
                          mpierror)
      ub= uu2(1)
      um= uu2(2)
      u1= ub + um
!
!
!******************************************************************
! ---- Diagnosis ----------------- On the major node ------------ *
!******************************************************************
!
      if(rank.eq.0) then
!
        open (unit=11,file=praefixc//'.11'//suffix2, &
                  status='unknown',position='append',form='formatted')
!
!* Accumulate in every step                                       Add
        call magnetiz (spx,spy,spz,g,wx1,wy1,wz1,wn1,u1,uav,wt1,np1,1)
!                                                                  +1
!* Energy absorption due to phase lag: 
!       - <m(t)*b(t)/b_0> = - integral m*b/b_0 dt /t 
!
      if(iwrt1.eq.0) then
        ss= 0
        do 735 i= 1,np1
        ss= ss +spx(i)*Bex +spy(i)*Bey +spz(i)*Bez  ! Bez= bz/B00
  735   continue
!
        del_en= del_en - g*ss/np1
!
! ------------------------------
!*  For FT13 output            *
! ------------------------------
        open (unit=13,file=praefixc//'.13'//suffix2, &
                  status='unknown',position='append',form='unformatted')
!
        do i= 1,np1
        spx4(i)= spx(i)
        spy4(i)= spy(i)
        spz4(i)= spz(i)
        end do
!
        do i= 1,np1+np2
        ch4(i)= ch(i)
        x4(i)= x(i)
        y4(i)= y(i)
        z4(i)= z(i)
        vx4(i)= vx(i)
        vy4(i)= vy(i)
        vz4(i)= vz(i)
        end do
!
        write(13) t8
        write(13) (spx4(i),spy4(i),spz4(i),i=1,np1)
        write(13) (ch4(i),x4(i),y4(i),z4(i),i=1,np1+np2)
        write(13) (vx4(i),vy4(i),vz4(i),i=1,np1+np2)
        close(13)
!
!**********************************
!*  Make xyz file for DS Viewer   *
!**********************************
!
        open (unit=23,file=praefixc//suffix2//'mg.xyz', &
                  status='unknown',position='append',form='formatted')
!                         +++++++
!
        write(23,'(i6)') np1+np2
        write(23,'(a30)') 'All atoms in the entire system'
!
        do i= 1,np1
        write(23,123) 'Fe',x(i),y(i),z(i)
  123   format(a2,3f12.6)
        end do
!
        do i= np1+1,np1+np2
        write(23,123) 'O ',x(i),y(i),z(i)
        end do
!
        close (23)
!
!**************************
!*  History plots         *
!**************************
!
        is= is +1
        if(is.gt.nhs) call rehist (rank)
!
        timeh(is)= t8
        do 737 i= 1,np1
        if(spec(i).eq.2) then  ! get one of Fe(2+)
          spinx(is)= spx(i)
          spinz(is)= spz(i)
          spin7(is)= sqrt(spx(i)**2 +spy(i)**2 +spz(i)**2)
          go to 738
        end if
  737   continue
  738   continue
!
        Bextx(is)= B00*Bex
        Bextz(is)= B00*Bez
!                                                   ic= 7: Get average
        call magnetiz (spx,spy,spz,g,wx1,wy1,wz1,wn1,u1,uav,wt1,np1,7)
!                                                                  +7
        magx(is)= wx1  ! m= -g*mue_B*s
        magy(is)= wy1  ! u= -m dot b = g*mue_B s*b 
        magz(is)= wz1
!
!
        U_Fe(is)= 0
        U_O(is) = 0
        ds1= 0
        ds2= 0
!
        do i= 1,np1
        U_Fe(is)= U_Fe(is) +0.5*fc3*m_Fe*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        ds1     = ds1 +sqrt((x(i)-x_0(i))**2 +(y(i)-y_0(i))**2 &
                                                     +(z(i)-z_0(i))**2)
        end do
!
        do i= np1+1,np1+np2
        U_O(is)= U_O(is) +0.5*fc3*m_O*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        ds2    = ds2 +sqrt((x(i)-x_0(i))**2 +(y(i)-y_0(i))**2 &
                                                     +(z(i)-z_0(i))**2)
        end do
!
        ds_Fe(is) = ds1/np1 
        ds_O(is)  = ds2/np2
!
        u1 = u1 + U_Fe(is) + U_O(is)
!
!% mpif90 -mcmodel=medium -fpic -o ay.out @spin_SMD5h.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3 &> log
!                                                                *
        Usys8(is)= u1/(np1+np2)  ! <Usys>= <U_O +U_Fe>
        Usys(is)=  u1/(np1+np2)
!
        uss(is) =  um/np1
        usb(is) =  ub/np1
!
        U_Fe(is)= U_Fe(is)/np1
        U_O(is) = U_O(is) /np2
!
        e_sp  = e_sp /(np1+np2)
        e_c_r = e_c_r/(np1+np2)
        e_LJ  = e_LJ /(np1+np2)
!
        t4 = t8
        Bex4= B00*Bex
        Bez4= B00*Bez
!
!--------------------
!* Spin Temperature
!--------------------
        do 741 k= 1,3
        sn1(k)= 0
        sx1(k)= 0
        sy1(k)= 0
        sz1(k)= 0
        sx2(k)= 0
        sy2(k)= 0
        sz2(k)= 0
  741   continue
!
        do 743 i= 1,np1
        csx= spx(i)/sqrt(spx(i)**2+spy(i)**2+spz(i)**2)
        csy= spy(i)/sqrt(spx(i)**2+spy(i)**2+spz(i)**2)
        csz= spz(i)/sqrt(spx(i)**2+spy(i)**2+spz(i)**2)
!
        if(spec(i).eq.1) then    ! Fe(3+)
          if(site(i).eq.1) k= 1  !  at (A)
          if(site(i).eq.2) k= 2  !  at (B)
        else
          k= 3                   ! Fe(+2) at (B)
        end if
!
        sn1(k)= sn1(k) +1
        sx1(k)= sx1(k) +csx
        sy1(k)= sy1(k) +csy
        sz1(k)= sz1(k) +csz
        sx2(k)= sx2(k) +csx**2
        sy2(k)= sy2(k) +csy**2
        sz2(k)= sz2(k) +csz**2
  743   continue
!
        do 747 k= 1,3
        tsx(is,k)= sqrt(sx2(k)/sn1(k) -(sx1(k)/sn1(k))**2)
        tsy(is,k)= sqrt(sy2(k)/sn1(k) -(sy1(k)/sn1(k))**2)
        tsz(is,k)= sqrt(sz2(k)/sn1(k) -(sz1(k)/sn1(k))**2)
  747   continue
!
! ------------------------------------------------------------
!* av_tsz() is used in the MC/MD step
!
        ss= 0.
        do 748 k= 1,3
        ss= ss +tsz(is,k)
  748   continue
        av_tsz(is)= ss/3.d0
!
        if(is.eq.5) tsz0= (av_tsz(3)+av_tsz(4)+av_tsz(5))/3.
! ------------------------------------------------------------
        open (unit=11,file=praefixc//'.11'//suffix2, &
                  status='unknown',position='append',form='formatted')
!
        conv(is)= sav/np1
        aitr(is)= iter
        psdt(is)= dtrhs
        sum_mb(is)= del_en
!
!       if(is.eq.1) then
        if(ft06_start) then
          ft06_start= .false.
!
          write(11,770)
  770     format(//,'# MD Run is performed #',/, &
            '      t8   it   is     itr  Usys      U_Fe      U_O     ',&
            '  conv      f*dts/m_i    v*dt      e_sp     e_c_r    ',   &
            ' e_LJ     magz   deV_x(Fe)   deV_x(O)   cpu(min)')
          write(11,*)
        end if
!
        bbw= B00 * Bez
        fdt4(is)= fdt8/vth_O
        vdt4(is)= vdt8
        idt4(is)= ifdt
! 
        write(11,771) t8,it,is,iter,Usys(is),U_Fe(is),U_O(is),conv(is), &
                 fdt8,vdt8,e_sp,e_c_r,e_LJ,magz(is),ds_Fe(is),ds_O(is), &
                 wtime/60.d0
  771   format('t=',f7.1,i8,i5,i4,1p6d10.2,3d10.2,0pf10.5,2f10.3,f8.2)
!
!                                                   ic= 0: Reset wx-wn
        call magnetiz (spx,spy,spz,g,wx1,wy1,wz1,wn1,u1,uav,wt1,np1,0)
!                                                                  +0
        close (11)
      end if
!
!
      if(iwrt2.eq.0) then
        open (unit=77,file=praefixc//'.77'//suffix2, &
                    status='unknown',position='append',form='formatted')
        time= t8  !<- real*4
!
        call plot_spin (x,y,z,spx,spy,spz,spec,site,np1)
        call distr_spin (spx,spy,spz,spec,site,np1)
!
        close (77)
      end if
!
!
      if(iwrt1.eq.0) then
        open (unit=12,file=praefixc//'.12'//suffix2, & ! usual= 2
                                  status='unknown',form='unformatted')
!
        write(12) t8,xleng,yleng,zleng,rcut,rcutC,Temp,TCurie, &
                  tmax,dt,cptot
        write(12) x,y,z,vx,vy,vz,ch,mass,ag,ep
        write(12) x_0,y_0,z_0,rintC
        write(12) spx,spy,spz,sp2,spx00,spy00,spz00,r_ij
        write(12) aspx,aspy,aspz
        write(12) Jint,u,spin2,spin3,n_MCsteps,tau_diss
        write(12) Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,toler,itermax
        write(12) dtwr,dtwr2,fw00,atsz,tsz0,av_tsz
        write(12) it,is,iwa,iwb,ir0
        write(12) np1,np2,spec,site,nintS,lintS,nintC,lintC,if_LJ
        write(12) i1,i2,i3,i4,disp_recv,cnt_recv,disp_recvC,cnt_recvC
        write(12) spinx,spinz,spin7,Bextx,Bextz,magx,magy,magz, &
                  Usys,conv,aitr,psdt,tfix,uss,usb,tsx,tsy,tsz,sum_mb, &
                  U_Fe,U_O,fdt4,vdt4,idt4,timeh
        close(12)
      end if
!
      close (11)
      end if
!
!******************************************************************
! ---- Diagnosis Ended ------------- On the major node ---------- *
!******************************************************************
!
        if(istop.eq.1) go to 7000
      go to 1000
! ************* End of the Loop ***********************************
!
 7000 continue
!
      if(rank.eq.0) then
!
        call date_and_time7 (cdate,ctime)
!
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,*) ' Job is finished  t8, cptot=',t8,cptot
        write(11,*) ' today = ',cdate
        write(11,*) ' time  = ',ctime

        write(11,*) ' Final: t8, tmax=',t8,tmax
        close (11)
!
        open (unit=77,file=praefixc//'.77'//suffix2, &
                    status='unknown',position='append',form='formatted')
        call lplots 
        call gclose
!
        close (77)
      end if
!
      return
      end subroutine spin_dynamics 
!
!
!-----------------------------------------------------------------------
      subroutine forces (fc1,fc2,fcLJ,alpha,e_sp,e_c_r,e_LJ, &
                         i1,i2,i3,i4,rank,np1,np2,it)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'param-spinRL5.h'
!
      real(C_DOUBLE) x,y,z,spx,spy,spz,ch,fx,fy,fz,fxc,fyc,fzc, &
                 Jint,r_ij,rintC,fc1,fc2,fcLJ,alpha,xleng,yleng,zleng,&
                 e_sp,e_c_r,e_LJ,djdr,r2,r,xx,yy,zz,            &
                 ff,pi,sqrtpi,erfc,forceV,rsi,snt,snt2,rLJ,ccel
      common/ewald3/ xleng,yleng,zleng

      integer(C_INT) n_of_j,i1(0:num_proc),i2(0:num_proc),i3(0:num_proc),&
                 i4(0:num_proc),rank,nintS,lintS,nintC,lintC, &
                 if_LJ,spec,site,np1,np2,i,j,k,l,it
      common/partcl1/ &
                 x(np0),y(np0),z(np0),spx(np0),spy(np0),spz(np0),ch(np0)
      common/partcl2/ &
                 fx(np0),fy(np0),fz(np0),fxC(np0),fyC(np0), &
                 fzC(np0),Jint(nbx,np0),r_ij(nbx,np0),rintC(nbx,np0)
      common/partcl3/ & 
                 nintS(np0),lintS(nbx,np0),nintC(np0),lintC(nbx,np0), &
                 if_LJ(np0),spec(np0),site(np0)
!
      real(C_DOUBLE) a1,a2,a3,a4,a5,pp,ar,tt
      parameter ( a1 = 0.254829592d0, a2 = -0.284496736d0 )
      parameter ( a3 = 1.421413741d0, a4 = -1.453152027d0 )
      parameter ( a5 = 1.061405429d0, pp =  0.3275911d0   )
!
      real(C_DOUBLE) a,b,c,d,bb
      parameter ( a =  1.8277098d3, &   ! e^2/a
                  b =  4.9254591d0, &   ! 1/a
                  c = -2.1362173d0, &   ! e^2/a^5
                  d =  7.4679860d1  )   ! e^2/a^11
!*---------------------------------------------------------------
!
      pi= 4.d0*atan(1.d0)
      sqrtpi= sqrt(pi)
!
      do i= 1,np1
      fx(i)= 0   ! note: not defined for i > np1
      fy(i)= 0
      fz(i)= 0
      end do
!
      do i= 1,np1+np2
      fxC(i)= 0
      fyC(i)= 0
      fzC(i)= 0
      end do
!
      e_sp = 0
      e_LJ = 0
      e_c_r= 0
!
!  Spin forces
! ####################################
!
      do i= i1(rank),i2(rank) 
      do k= 1,nintS(i) 
      j= lintS(k,i)
!
! 2*0.5 is omitted in djdr
!     fc1= J00*(t_unit/a_unit)**2/m_unit          ! spin force
!     fc2= (e_unit*t_unit)**2/(m_unit*a_unit**3)  ! coulomb force
!
!* Fold back
      xx= x(i) - x(j)
      yy= y(i) - y(j)
      zz= z(i) - z(j)
      xx= xx -nint(xx/xleng)*xleng
      yy= yy -nint(yy/yleng)*yleng
      zz= zz -nint(zz/zleng)*zleng
!
      r = sqrt(xx**2 +yy**2 +zz**2)
      djdr= -fc1 *Jint(k,i)*(r_ij(k,i)/r**2)* &
                  (spx(i)*spx(j) +spy(i)*spy(j) +spz(i)*spz(j))
!
      fx(i)= fx(i) - djdr*xx/r
      fy(i)= fy(i) - djdr*yy/r
      fz(i)= fz(i) - djdr*zz/r
!
      e_sp= e_sp +djdr
      end do
      end do
!
!* Coulomb(LJ) forces
! #######################################
      do i= i3(rank),i4(rank)
      do k= 1,nintC(i)
      j= lintC(k,i)
!
      xx= x(i) - x(j)
      yy= y(i) - y(j)
      zz= z(i) - z(j)
      xx= xx -nint(xx/xleng)*xleng
      yy= yy -nint(yy/yleng)*yleng
      zz= zz -nint(zz/zleng)*zleng
!
      r2= xx**2 +yy**2 +zz**2
      r = sqrt(r2)
!
      rsi = 1/r2
      snt = rsi*rsi*rsi
!
!*  Lennard-Jones force only for the Fe-O pair
!     ref. J.Rustad, B.Hay, and J.Halley, J.Chem.Phys. 102, 427 (1995).
!             ***
!
      if(spec(i).eq.0) then           ! O-
        if(site(j).eq.1) go to 300    !  -Fe(A)
        if(site(j).eq.2) go to 400    !  -Fe(B) +3, +2
      else                            ! Fe-
        if(spec(j).eq.0) then         !   -O
          if(site(i).eq.1) go to 300  ! Fe on A-site
          if(site(i).eq.2) go to 400  ! Fe on B-site
        end if
      end if
      go to 500
!
  300 continue
      e_LJ = e_LJ +fc2*( a*exp(-b*r) +(c +d*snt)*snt )
      ccel =  fc2*( a*b*exp(-b*r)/r +(6*c +12*d*snt)*snt/r2 )
      go to 600
!
  400 continue
      if(.true.) go to 300  !<-- 
!     +++++++++ 
!
      if(spec(i).eq.1 .or. spec(j).eq.1) then
        bb = 0.85d0*b  ! +3
      else
        bb = 0.9d0*b   ! +2
      end if
!
      e_LJ = e_LJ +fc2*( a*exp(-bb*r) +(c +d*snt)*snt )
      ccel =  fc2*( a*bb*exp(-bb*r)/r +(6*c +12*d*snt)*snt/r2 )
      go to 600
!
!* Fe-Fe or O-O pair: LJ forces
  500 continue
!
      ccel= 0
      if(if_LJ(i).eq.1) then   ! Limit the pair to r(t=0) < rcutC
!        ++++++++
!
        snt2= 2**(1.d0/6.d0)
        rLJ = r/(rintC(k,i)/snt2)   ! Minimum at rintC(k,i) for this pair
!                *****
        rsi = 1/dmax1(rLJ**2, 0.81d0)
        snt = rsi*rsi*rsi
!
        e_LJ = e_LJ +fc2*(d/12.d0)*snt*(snt -1.d0)
        ccel =       fc2*d*snt*(snt -0.5d0)/r2
      end if
!            <-- also go to 600 !
!
  600 continue
!*  Coulomb force in Ewald sum
!
      ar   = alpha*r
      tt   = 1.d0 / ( 1.d0 + pp * ar )
      erfc = tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))
!
      e_c_r  = e_c_r + fc2*ch(i)*ch(j)* erfc*dexp(-ar**2)/r
      forceV = fc2*ch(i)*ch(j)*(erfc/r +2*alpha/sqrtpi)*dexp(-ar**2)/r2
!        1         2         3         4         5         6         71
!
      fxC(i)= fxC(i) +(ccel +forceV)*xx 
      fyC(i)= fyC(i) +(ccel +forceV)*yy
      fzC(i)= fzC(i) +(ccel +forceV)*zz
      end do
      end do
! #######################################
!
      return
      end subroutine forces 
!
!
!-----------------------------------------------------------------
      subroutine init (x,y,z,ch,spx,spy,spz,sp2,spec,site,     &
                       xleng0,yleng0,zleng0,np1,np2,np10,np20, &
                       np000,i1x,i2x,i1y,i2y,i1z,i2z,rank,size)
!-----------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'mpif.h'
      include    'param-spinRL5.h'
!
      real(C_DOUBLE) xo(100),yo(100),zo(100),cho(100)
      real(C_DOUBLE) x(np0),y(np0),z(np0),ch(np0),spx(np0),spy(np0),&
                 spz(np0),sp2(np0),spin2,spin3,Jaa,Jbb,Jab,     &
                 Bapx,Bapy,Bapz,tau_b,tau_diss,Temp,TCurie,     &
                 xleng,yleng,zleng,pi,pi2,th,ph,xx,yy,zz,       &
                 t8,dt,tmax,cptot,xleng0,yleng0,zleng0,         &
                 a(3),b(3),c(3),qFe2,qFe3,qo
      real(C_DOUBLE) eps,ranff
      integer(C_INT) spec(np0),site(np0),np1,np2,rank,size,  &
                 np000,np10,np20,mx,my,mz,i,l,j,lo,          & 
                 ia,ja,ka,ic,jc,kc,i1,i2,j1,j2,k1,k2,        &
                 nFe2,nFe3,ierror,i1x,i2x,i1y,i2y,i1z,i2z
      character  char*2,text1*173  !text1*125
!
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b, &
                    tau_diss,Temp,TCurie
      common/parm3/ t8,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
      common/parm4/ mx,my,mz
!
!* Three axes of the crystal
      a(1)= xleng0
      a(2)= 0.000000000 
      a(3)= 0.000000000 
      b(1)= 0.000000000 
      b(2)= yleng0
      b(3)= 0.000000000 
      c(1)= 0.000000000 
      c(2)= 0.000000000 
      c(3)= zleng0
!
      if(rank.eq.0) then
!       open (unit=17,file='magnetite.xyz', & ! *125
        open (unit=17,file='magnetite8.xyz', & ! *173
                              status='old',form='formatted')
!
        read(17,'(i6)') np000   !<<-- the number in a cell
        read(17,'(a125)') text1
      end if
!
      call mpi_bcast (np000,  1,mpi_integer,  0,mpi_comm_world,ierror)
!     call mpi_bcast (text1,125,mpi_character,0,mpi_comm_world,ierror)
!
      pi= 4.d0*atan(1.d0)
      pi2= 2*pi
!
      spin2= 2.0d0  ! Fe 2+
      spin3= 2.5d0  ! Fe 3+
!
      qFe2= 2.0d0
      qFe3= 3.0d0
      qo = -2.0d0
!
      i= 0
      j= 0
      l= 0
!
!  magnet3.xyz
!  FT17: O(1)...O(32) B(1)...B(16) A(1)...A(8) 
!
  100 l= l +1
      if(l.gt.np000) go to 700
!
!* Balance charges on cubic sub-lattice
!
      if(l.gt.32 .and. mod(l,4).eq.1) then
        nFe2= 0
        nFe3= 0
      end if
!
      if(rank.eq.0) then
        read(17,101) char,xx,yy,zz
  101   format(a2,f12.0,f12.0,f12.0)
      end if
!
      call mpi_bcast (char,2,mpi_character,     0,mpi_comm_world,ierror)
      call mpi_bcast (xx,1,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast (yy,1,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast (zz,1,mpi_double_precision,0,mpi_comm_world,ierror)
!
! Oxygen
      if(l.le.32) then  
        j= j +1
!       +++++++
!
        xo(j)= xx
        yo(j)= yy
        zo(j)= zz
        cho(j)= qo
        spec(j)= 0     ! oxygen
        go to 100
      end if
!
! Fe
      i= i+1
!     +++++++
!
      x(i)= xx
      y(i)= yy
      z(i)= zz
!
      if(l.ge.49 .and. l.le.56) then
        ch(i)= qFe3
        sp2(i)= spin3
        spec(i)= 1
        site(i)= 1      ! A site (8)
      else
        eps= ranff(0.d0)
        if(eps.gt.0.5d0) then
          if(nFe3.lt.2) go to 110  ! ...2 Fe(3+) are allowed on sub lattice
          go to 111
        else
          if(nFe2.lt.2) go to 111  ! ...2 Fe(2+) are allowed
          go to 110
        end if
!
  110   nFe3= nFe3 +1
        ch(i)= qFe3
        sp2(i)= spin3
        spec(i)= 1       ! 3+
        site(i)= 2       ! B site (8+8)
        go to 130
!
  111   nFe2= nFe2 +1
        ch(i)= qFe2
        sp2(i)= spin2
        spec(i)= 2       ! 2+
        site(i)= 2       ! b site (8+8)
!
  130   continue
      end if
!  
      th=  pi*ranff(0.d0)
      ph= pi2*ranff(0.d0)
      spx(i)= sp2(i)*sin(th)*cos(ph)
      spy(i)= sp2(i)*sin(th)*sin(ph)
      spz(i)= sp2(i)*cos(th)
      go to 100
!
  700 np10= i
      np20= j

      if(rank.eq.0) then
        close (17)  !!<<--
!
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,703) np10,np20,nFe2,nFe3+8
  703   format(/,' Fundamental domain of irons and oxygens ',/, &
                 '   #Fe,:#O, Fe(2+), Fe(3+)=',2i4,2i4)
!
        close (11)
      end if
!
!------------------------
!* Expand the system
!------------------------
!  (mx,my,mz) the number of domains - must be odd numbers
!
      l= np10  !<-- the number of Fe in the cell
!
      xleng= mx*xleng0
      yleng= my*yleng0
      zleng= mz*zleng0
!
      if(mx.eq.1) then
        i1= 0
        i2= 0
      else
        if(mod(mx,2).eq.1) then
          i1= -mx/2
          i2=  mx/2
        else
          i1= -mx/2 +1
          i2=  mx/2
        end if
      end if
!
      if(my.eq.1) then
        j1= 0
        j2= 0
      else
        if(mod(my,2).eq.1) then
          j1= -my/2
          j2=  my/2
        else
          j1= -my/2 +1
          j2=  my/2
        end if
      end if
!
      if(mz.eq.1) then
        k1= 0
        k2= 0
      else
        if(mod(mz,2).eq.1) then
          k1= -mz/2
          k2=  mz/2
        else
          k1= -mz/2 +1
          k2=  mz/2
        end if
      end if
!
      i1x= i1
      i2x= i2
      i1y= j1
      i2y= j2
      i1z= k1
      i2z= k2
!
      if(mx.eq.1) go to 301
!     +++++++++++++++++++++  i1=i2=0 n/g on infiniband 
!
      do 300 ia= i1,i2
      do 300 ja= j1,j2
      do 300 ka= k1,k2
      if(ia.eq.0 .and. ja.eq.0 .and. ka.eq.0) go to 300
!        ++++++++ Fundamental domain - skipped !
!
! Fe
      do 400 i= 1,np10 
      l = l +1
!
      x(l)= x(i) +ia*a(1) +ja*b(1) +ka*c(1)
      y(l)= y(i) +ia*a(2) +ja*b(2) +ka*c(2)
      z(l)= z(i) +ia*a(3) +ja*b(3) +ka*c(3)
!
!  Fe: b(1)...b(16) a(17)...a(24) 
      if(i.le.16 .and. mod(i,4).eq.1) then
        nFe2= 0
        nFe3= 0
      end if
!
! A site
      if(i.gt.16) then  
        ch(l) = ch(i)
        sp2(l) = sp2(i)
        spec(l)= spec(i)
        site(l)= site(i)
        go to 430
      else
!
! B site
!* shuffle Fe(2+)/Fe(3+)= 8/8
!
        if(ranff(0.d0).gt.0.5d0) then
          if(nFe3.lt.2) go to 410  ! Accomodate 2 on sub-lattice
          go to 411
        else
          if(nFe2.lt.2) go to 411
          go to 410
        end if
      end if
!
  410 nFe3= nFe3 +1
      ch(l) = qFe3
      sp2(l)= spin3
      spec(l)= 1
      site(l)= 2       ! B site (16)
      go to 430
!
  411 nFe2= nFe2 +1
      ch(l) = qFe2
      sp2(l)= spin2
      spec(l)= 2
      site(l)= 2       ! B site (16)
!
  430 continue
      th=  pi*ranff(0.d0)
      ph= pi2*ranff(0.d0)
      spx(l)= sp2(i)*sin(th)*cos(ph)
      spy(l)= sp2(i)*sin(th)*sin(ph)
      spz(l)= sp2(i)*cos(th)
  400 continue
  300 continue
!
  301 continue
! --------------
      np1= l     !<<-- the total number of irons
! --------------
!
!% mpif90 -mcmodel=medium -fpic -o a.out @spin_SMD5a.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3 &> log
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
        write(11,*) 'np000,np1=',np000,np1
        close (11)
      end if
!
! Oxygen
!
      j= 0
!
      do 500 i= 1,np20
      l= l +1
      j= j +1
      x(l)= xo(i)
      y(l)= yo(i)
      z(l)= zo(i)
      ch(l)= qo
      spec(l)= 0
  500 continue
!
      if(mx.eq.1) go to 601
!     +++++++++++++++++++++  i1=i2=0 n/g on infiniband 
      do 600 ia= i1,i2
      do 600 ja= j1,j2
      do 600 ka= k1,k2
      if(ia.eq.0 .and. ja.eq.0 .and. ka.eq.0) go to 600
!
      do 630 i= 1,np20
      l= l +1
      j= j +1
      x(l)= xo(i) +ia*a(1) +ja*b(1) +ka*c(1)
      y(l)= yo(i) +ia*a(2) +ja*b(2) +ka*c(2)
      z(l)= zo(i) +ia*a(3) +ja*b(3) +ka*c(3)
      ch(l)= qo
      spec(l)= 0
  630 continue
!
  600 continue
  601 continue
!
! --------------
      np2= j     !<<-- the total number of oxygens
! --------------
!
!******************************************************
!*  File I/O is allowed only by the master node       *
!******************************************************
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,707) size,np1,np2
  707   format(' /init/: Number of processors=',i3,/, &
               '         the number of np(Fe),np(O)=',2i6,/)
        close (11)
! +++++
!
        open (unit=18,file='mag3.xyz',form='formatted')
!
        write(18,'(i6)') np1,np2
        write(18,'(a80)') 'magnetite in 3d'
!
        do 105 i= 1,np1
        j= np1 +i
        write(18,107) 'Fe ',i,x(i),y(i),z(i),ch(i), &
                      ' O ',j,x(j),y(j),z(j),ch(j)
  107   format(a3,'i=',i4,3f10.5,' q=',f5.1,2x, &
               a3,'i=',i4,3f10.5,' q=',f5.1)
  105   continue
!
        do 115 i= np1+1,np2
        j= np1 +i
        write(18,108) ' o ',j,x(j),y(j),z(j),ch(j)
  108   format(49x,a3,'i=',i4,3f10.5,' q=',f5.1)
  115   continue
!
        close (18)
!
!* make xyz file for ds viewer
!
        open (unit=23,file=praefixc//suffix2//'mg.xyz', &
                                 status='unknown',form='formatted')
!
        write(23,'(i6)') np1+np2
        write(23,'(a30)') 'all atoms in the entire system'
!
        do i= 1,np1
        write(23,123) 'Fe',x(i),y(i),z(i)
  123   format(a2,3f12.6)
        end do
!
        do i= np1+1,np1+np2
        write(23,123) 'o ',x(i),y(i),z(i)
        end do
!
        close (23)
! +++++
      end if
!
      return
      end subroutine init
!
!
!cccccccccccccc   p3m (fortran 77)  ccccccccccccccccccccccccccccccccccc
!                                 26.10.99 
!                                     Motohiko Tanaka, Christian Holm
!  Caling sequences:
!
!     call  p3m_init (length,alpha0,ip)
!     call  p3m_perform (coox,cooy,cooz,q,fx,fy,fz,e_Coulomb_p3m)
!
!/*---------------------------------------------------------------------
! subunit:  p3m_v2.f  (fortran 90)
! 
!       version   20.10.1999 (ch) 
!       corrected 26.10.1999 (mt)
!                 23.06.2000 (mt)
!
! version:  20 january 1999
! Author:   Markus Deserno
!
!    brillouin is now a parameter
!    maxinterpol --> mintpol
!    floor --> defined: must be stepwise at x= 0.
!    dround --> = nint (round off to a nearest integer). 
!/*---------------------------------------------------------------------
      subroutine p3m_perform (x,y,z,ch,fkx,fky,fkz,np,iwrt1, &
                              rank,e_Coulomb_p3m)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit    none
!
      include     'mpif.h'
      include     'fftw3.f'
!
      include     'param-spinRL5.h'
      integer(C_INT) np,iwrt1,rank,ierror
!
      real(C_DOUBLE) x(0:np0-1),y(0:np0-1),z(0:np0-1),ch(0:np0-1), &
                  fkx(0:np0-1),fky(0:np0-1),fkz(0:np0-1),          &
                  ecl,e_Coulomb_p3m,wupi,pi
!
      integer(C_INT) xarg(0:np0-1),yarg(0:np0-1),zarg(0:np0-1), &
                  i,j,k,m,meshmaskx,meshmasky,meshmaskz, &   
                  gi0,gi1,gi2,xpos,ypos,zpos,            &
                  assignshiftx,assignshifty,assignshiftz,&
                  qzahl,m0
      real(C_DOUBLE) d1,hix,hiy,hiz,mi2,modadd1,modadd2,t1,t2,t3, &
                  sum_q_2, sum_q2
!     
!-----------
      real(C_DOUBLE) alpha,fc2,xleng,yleng,zleng
      real(C_DOUBLE) meshiftx,meshifty,meshiftz,dnx,dny,dnz
      integer(C_INT) PP
      common/ewald1/ alpha,fc2
      common/ewald2/ PP 
      common/ewald3/ xleng,yleng,zleng
!
      real(C_DOUBLE) coop,qp,ql,ghat,intcaf
      integer(C_INT) global,g
!-----------
      common/coordi/ coop(0:np0-1,0:2),qp(0:np0-1)
      common/pindex/ global(0:np0-1),g(0:np0-1,0:2)
      common/prmesh/ ql(0:p_max**3-1,0:np0-1)
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy), &
                     meshiftz(0:meshz),dnx(0:meshx),dny(0:meshy), &
                     dnz(0:meshz)
!
      common/influf/ ghat(0:meshx,0:meshy,0:meshz)
      common/intcaf/ intcaf(0:p_max,0:2*mintpol+1)
!
      complex*16         qq_c,phi_x_c,phi_y_c,phi_z_c
      common/fft001/     qq_c(0:meshx/2, 0:meshy-1,0:meshz-1), & ! fftw-style
                      phi_x_c(0:meshx/2, 0:meshy-1,0:meshz-1), & !  arrays
                      phi_y_c(0:meshx/2, 0:meshy-1,0:meshz-1), &
                      phi_z_c(0:meshx/2, 0:meshy-1,0:meshz-1)
!
      real(C_DOUBLE)  qq,phi_x,phi_y,phi_z
      common/fft002/     qq(0:meshx-1,0:meshy-1,0:meshz-1), & 
                      phi_x(0:meshx-1,0:meshy-1,0:meshz-1), & 
                      phi_y(0:meshx-1,0:meshy-1,0:meshz-1), & 
                      phi_z(0:meshx-1,0:meshy-1,0:meshz-1)  
!-----------
!
! -------------------------
!*  FFTW3 - Fortran header
! -------------------------
!       em64t: 64-bit pointer
!                 **** ****
      integer*8   plan,pinv
!
      logical     first
      common/fftsav/ plan,pinv,first
      data        first /.true./
!
      real(C_DOUBLE) fft_scale2
      complex*16  ei,eigq
!
! ---------------------------
!*  Prepare for fftw calls
! ---------------------------
!
      if(first) then
        first= .false.
!
!       call dfftw_init_threads (ierror)
!       call dfftw_plan_with_nthreads (nthreads)
!
        call dfftw_plan_dft_r2c_3d (plan,meshx,meshy,meshz,qq,qq_c, &
                                    fftw_estimate)
!
        call dfftw_plan_dft_c2r_3d (pinv,meshx,meshy,meshz,qq_c,qq, &
                                    fftw_estimate)
      end if
!
!
      ei = dcmplx(0.d0, 1.d0)
! --------------------------------
      pi = 4.d0*atan(1.d0)
      meshmaskx = meshx-1     
      meshmasky = meshy-1     
      meshmaskz = meshz-1     
!
      hix = meshx / xleng
      hiy = meshy / yleng
      hiz = meshz / zleng
!
      mi2 = 2.d0*dfloat(mintpol)
      assignshiftx = meshx -(pp -1)/2
      assignshifty = meshy -(pp -1)/2
      assignshiftz = meshz -(pp -1)/2
!
      qzahl = 0
      sum_q_2 = 0.d0
      sum_q2  = 0.d0
!
      do 100 i= 0,np-1
      if (dabs(ch(i)) .gt. 1.d-5) then 
        coop(qzahl, 0) = x(i) - nint(x(i)/xleng -0.5d0)*xleng 
        coop(qzahl, 1) = y(i) - nint(y(i)/yleng -0.5d0)*yleng
        coop(qzahl, 2) = z(i) - nint(z(i)/zleng -0.5d0)*zleng
!
        qp(qzahl) = ch(i)
        sum_q_2 = sum_q_2 + qp(qzahl)
        sum_q2  = sum_q2  + qp(qzahl)**2
!
        global(qzahl) = i
        qzahl= qzahl + 1
      end if
!
      fkx(i)= 0
      fky(i)= 0
      fkz(i)= 0
  100 continue
!
      sum_q_2 = sum_q_2 **2
!
      do 200 i= 0,meshx-1 
      do 200 j= 0,meshy-1 
      do 200 k= 0,meshz-1 
      qq(i,j,k) = 0.d0
  200 continue
!
      if(pp.eq.3) then
         modadd1 = 0.0d0
         modadd2 = 0.5d0
      else
         write(6,*) "Error in function 'p3m_perform':"
         write(6,*) "charge assignment order p=",pp," unknown."
         write(6,*) "program terminated."
         call exit(1)
      end if
!
!
      do 300 i= 0, qzahl-1
!*
      d1  = coop(i,0)*hix + modadd1
      gi0 = int(d1 + modadd2) + assignshiftx
      g(i,0) = gi0 
      xarg(i) = int( (d1 - nint(d1) + 0.5d0)*mi2 )
!      
      d1  = coop(i,1)*hiy + modadd1 
      gi1 = int(d1 + modadd2) + assignshifty
      g(i,1) = gi1
      yarg(i) = int( (d1 - nint(d1) + 0.5d0)*mi2 )
!      
      d1  = coop(i,2)*hiz + modadd1 
      gi2 = int(d1 + modadd2) + assignshiftz
      g(i,2) = gi2 
      zarg(i) = int( (d1 - nint(d1) + 0.5d0)*mi2 )
!
      m0= -1
!*
      do 330 j = 0, PP-1
      do 330 k = 0, PP-1
      do 330 m = 0, PP-1
      t1 = qp(i) * intcaf(j,xarg(i))
      t2 = t1 *    intcaf(k,yarg(i))
      t3 = t2 *    intcaf(m,zarg(i))
!    
      m0= m0 + 1
      ql(m0,i) = t3          ! Assignment factor.
!      
      xpos = iand( (g(i,0) + j), meshmaskx)
      ypos = iand( (g(i,1) + k), meshmasky)
      zpos = iand( (g(i,2) + m), meshmaskz)
!    
      qq(xpos,ypos,zpos) = qq(xpos,ypos,zpos) + ql(m0,i)
  330 continue
  300 continue
!
!
      fft_scale2= dfloat(meshx*meshy*meshz)
      call dfftw_execute_dft_r2c (plan,qq,qq_c)
!
! 
      ecl = 0.d0
!
      do 700  i= 0, meshx/2
      do 700  j= 0, meshy-1
      do 700  k= 0, meshz-1
      ecl = ecl + ghat(i,j,k)*cdabs(qq_c(i,j,k))**2/fft_scale2 
  700 continue
!
!-----------------------------------------------------------------
!*  Prefactor is for fft3r.    
!      checked by Madelung constant                 10/06/2000
!-----------------------------------------------------------------
!
      e_Coulomb_p3m = fc2 * ecl* (1.d0/meshx) *xleng/(4.d0*pi)
!
!     wupi = dsqrt(pi)
!     e_Coulomb_p3m = e_Coulomb_p3m - 
!    *         l *prefactor *( sum_q2 *alpha /wupi
!    *                       + sum_q_2 *pi /(2.d0*l**3*alpha**2) )
!
!
      do 800 j= 0,meshy-1
      do 800 k= 0,meshz-1
      phi_x_c(0,j,k) = 0.d0
      phi_y_c(0,j,k) = 0.d0
      phi_z_c(0,j,k) = 0.d0
  800 continue
!
      do 810 i= 0,meshx/2
      do 810 j= 0,meshy-1
      do 810 k= 0,meshz-1
      eigq = ei*ghat(i,j,k)*qq_c(i,j,k)/fft_scale2
!
      phi_x_c(i,j,k) = eigq*dnx(i)
      phi_y_c(i,j,k) = eigq*dny(j)
      phi_z_c(i,j,k) = eigq*dnz(k)
  810 continue
!    
!*  For diagnosis.
!    
      do 860 i= 0,meshx/2
      do 860 j= 0,meshy-1
      do 860 k= 0,meshz-1
      qq_c(i,j,k) = ghat(i,j,k)*qq_c(i,j,k)/fft_scale2
  860 continue
!
      if(iwrt1.eq.0) then
        call dfftw_execute_dft_c2r (plan,qq_c,qq)
      end if
!
      call dfftw_execute_dft_c2r (plan,phi_x_c,phi_x)
      call dfftw_execute_dft_c2r (plan,phi_y_c,phi_y)
      call dfftw_execute_dft_c2r (plan,phi_z_c,phi_z)
!
!
      do 900 i = 0, qzahl-1
!
      m0= -1
      do 930 j = 0, PP-1
      do 930 k = 0, PP-1
      do 930 m = 0, PP-1
      xpos = iand( (g(i,0) + j), meshmaskx)
      ypos = iand( (g(i,1) + k), meshmasky)
      zpos = iand( (g(i,2) + m), meshmaskz)
!    
      m0 = m0 +1
      d1 = fc2 * ql(m0,i)
!    
      fkx(global(i)) = fkx(global(i)) - d1*phi_x(xpos,ypos,zpos)
      fky(global(i)) = fky(global(i)) - d1*phi_y(xpos,ypos,zpos)
      fkz(global(i)) = fkz(global(i)) - d1*phi_z(xpos,ypos,zpos)
  930 continue
  900 continue
!
      return
      end subroutine p3m_perform 
!
!
!/*---------------------------------------------------------------------
      subroutine perform_aliasing_sums (nx,ny,nz,nominatorx,nominatory,& 
                                        nominatorz,denominator)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit    none
!
      include    'param-spinRL5.h'
!
      integer(C_INT) brillouin
      parameter  (brillouin=1)
!
      integer(C_INT) PP
      real(C_DOUBLE) alpha,fc2,xleng,yleng,zleng
      real(C_DOUBLE) meshiftx,meshifty,meshiftz,dnx,dny,dnz, &
                     dmeshx,dmeshy,dmeshz,pi
      common/ewald1/ alpha,fc2
      common/ewald2/ PP
      common/ewald3/ xleng,yleng,zleng
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy), &
                     meshiftz(0:meshz),dnx(0:meshx),dny(0:meshy), &
                     dnz(0:meshz)
!
      integer(C_INT) nx,ny,nz,mx,my,mz
      real(C_DOUBLE) nominatorx,nominatory,nominatorz,denominator, &
                  s1,s2,s3,fak1x,fak1y,fak1z,fak2,fak3, & 
                  nmx,nmy,nmz,nm2,expo,exponent_limit,sinc
!
      exponent_limit = 30.d0
!
      dmeshx= dfloat(meshx)
      dmeshy= dfloat(meshy)
      dmeshz= dfloat(meshz)
!
      pi = 4.d0*atan(1.d0)
!
!     fak2 = (pi/alpha)**2/leng**2
      fak2 = (pi/alpha)**2
!
      nominatorx = 0.d0
      nominatory = 0.d0
      nominatorz = 0.d0
      denominator = 0.d0
!
      do 100 mx = -brillouin, brillouin
      nmx = meshiftx(nx) + dmeshx*mx
      s1  = sinc(nmx/dmeshx)**(2*pp)
!
        do 200 my = -brillouin, brillouin
        nmy = meshifty(ny) + dmeshy*my
        s2  = s1* sinc(nmy/dmeshy)**(2*pp)
!
          do 300 mz = -brillouin, brillouin
          nmz = meshiftz(nz) + dmeshz*mz
          s3  = s2* sinc(nmz/dmeshz)**(2*pp)
!
          denominator = denominator + s3
          nm2 = nmx**2 + nmy**2 + nmz**2
!
          expo= fak2*((nmx/xleng)**2 +(nmy/yleng)**2 +(nmz/zleng)**2)
          if (expo .lt. exponent_limit) then
              fak3 =  s3* dexp(-expo)/nm2 
          else
              fak3 = 0.d0
          end if
!
          nominatorx = nominatorx + fak3 * nmx
          nominatory = nominatory + fak3 * nmy
          nominatorz = nominatorz + fak3 * nmz
  300     continue
  200   continue
  100 continue
!
      return
      end subroutine perform_aliasing_sums 
!
!
!/*---------------------------------------------------------------------
      subroutine calculate_differential_operator (rank)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit    none
!
      include     'param-spinRL5.h'
      integer(C_INT) rank,i,j,k
!-----------
      real(C_DOUBLE) meshiftx,meshifty,meshiftz,dnx,dny,dnz, &
                     dmeshx,dmeshy,dmeshz
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy), &
                     meshiftz(0:meshz),dnx(0:meshx),dny(0:meshy), &
                     dnz(0:meshz)
!                                       ******  defined here.
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,*) ' - Calculating calculate_differential_operator'
        close (11)
      end if
!
      dmeshx= dfloat(meshx)
      dmeshy= dfloat(meshy)
      dmeshz= dfloat(meshz)
!
      do i= 0,meshx-1
      dnx(i) = dfloat(i) - nint(dfloat(i)/dmeshx)*dmeshx
      end do
!
      do j= 0,meshy-1
      dny(j) = dfloat(j) - nint(dfloat(j)/dmeshy)*dmeshy
      end do
!
      do k= 0,meshz-1
      dnz(k) = dfloat(k) - nint(dfloat(k)/dmeshz)*dmeshz
      end do
!
      dnx(meshx/2) = 0.d0
      dny(meshy/2) = 0.d0
      dnz(meshz/2) = 0.d0
!
      return
      end subroutine calculate_differential_operator
!
!
!/*---------------------------------------------------------------------
      subroutine calculate_influence_function (rank)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit    none
!
      include     'param-spinRL5.h'
      integer(C_Int) rank
!-----------
      real(C_DOUBLE) alpha,fc2,ghat,xleng,yleng,zleng,leng2
      real(C_DOUBLE) meshiftx,meshifty,meshiftz,dnx,dny,dnz
      integer(C_Int) PP
      common/ewald1/ alpha,fc2
      common/ewald2/ PP
      common/ewald3/ xleng,yleng,zleng
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy), &
                     meshiftz(0:meshz),dnx(0:meshx),dny(0:meshy), &
                     dnz(0:meshz)
      common/influf/ ghat(0:meshx,0:meshy,0:meshz)
!                    **************************  defined here.
!
      integer(C_Int) nx,ny,nz
      real(C_DOUBLE) ddnx,ddny,ddnz,dn2,fak1,fak2,fak3, &
                  nominatorx,nominatory,nominatorz,denominator
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,*) ' - Calculating influence function with parameters'
        write(11,601) meshx,meshy,meshz,pp,alpha,xleng,yleng,zleng
  601   format('   meshx, meshy, meshz=',3i5,/,  &
               '   p=',i5,/,'   alpha=',d20.12,/,&
               ',  xleng, yleng, zleng=',3d20.12,/)
        close (11)
      end if
! 
!     fak1= 2.d0*dfloat(mesh*mesh*mesh)/leng**2
!
!         write(11,*) 'p3m - aliasing sums for different nx,ny,nz:'
!  
      do 100 nx = 0, meshx-1
       do 200 ny = 0, meshy-1
        do 300 nz = 0, meshz-1
        if ( (nx.eq.0).and.(ny.eq.0).and.(nz.eq.0)) then
           ghat(nx,ny,nz)= 0.d0 
        else
           call perform_aliasing_sums (nx,ny,nz,nominatorx,nominatory, &
                                       nominatorz,denominator)
           ddnx = dnx(nx)
           ddny = dny(ny)
           ddnz = dnz(nz)
!
           dn2 = ddnx**2 + ddny**2 + ddnz**2
!  
           if (dn2 .gt. 1.d-7) then
             ghat(nx,ny,nz) = 2.d0*( &
                           meshx*ddnx*nominatorx/xleng**2 &
                         + meshy*ddny*nominatory/yleng**2 &
                         + meshz*ddnz*nominatorz/zleng**2 ) &
                         /( dn2 * denominator**2 )
           else 
             ghat(nx,ny,nz) = 0.d0
           end if
        end if
!
  300   continue
  200  continue
  100 continue
!
      return
      end subroutine calculate_influence_function 
!
!
!/*---------------------------------------------------------------------
      subroutine interpol_charge_assign_function (rank)
!/*---------------------------------------------------------------------
!*  p = 3 only
! ---------------
      use, intrinsic :: iso_c_binding
      implicit    none
!
      include       'param-spinRL5.h'
      integer(C_INT) PP,rank
      common/ewald2/ PP
!
      integer(C_INT) i
      real(C_DOUBLE) intcaf, dinterpol, x
      common/intcaf/ intcaf(0:p_max,0:2*mintpol+1)
!                    ***************************** defined here.
!
      if(PP.ne.3) then
        write(11,*) ' ## p = 3 only ###'
        stop
      end if
!
      dinterpol= dfloat(mintpol)
!
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,601) PP
  601   format(/,' - interpolating the order-',i1,' charge assignment', &
               ' function')
        close (11)
      end if
!
      do i= -mintpol, mintpol
      x= i/(2.d0*dinterpol)
      intcaf(0, i+mintpol) = 0.50d0*(0.5d0 - x)**2
      intcaf(1, i+mintpol) = 0.75d0 - x*x
      intcaf(2, i+mintpol) = 0.50d0*(0.5d0 + x)**2
      end do
!
      return
      end subroutine interpol_charge_assign_function
!
!
!/*---------------------------------------------------------------------
      subroutine calculate_meshift (rank)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit   none
!
      include    'param-spinRL5.h'
      integer(C_INT) rank,i,j,k
!-----------
      real(C_DOUBLE) meshiftx,meshifty,meshiftz,dnx,dny,dnz, &
                     dmeshx,dmeshy,dmeshz
      common/mesh01/ meshiftx(0:meshx),meshifty(0:meshy), &
                     meshiftz(0:meshz),dnx(0:meshx),dny(0:meshy), &
                     dnz(0:meshz)
!-----------
! 
      if(rank.eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2, &
                    status='unknown',position='append',form='formatted')
!
        write(11,*) ' - Calculating mesh-shift'
        close (11)
      end if
!  
      dmeshx= dfloat(meshx)
      dmeshy= dfloat(meshy)
      dmeshz= dfloat(meshz)
!
      do i= 0, meshx-1
      meshiftx(i) = i - nint(i/dmeshx)*dmeshx
      end do
!
      do j= 0, meshy-1
      meshifty(j) = j - nint(j/dmeshy)*dmeshy
      end do
!
      do k= 0, meshz-1
      meshiftz(k) = k - nint(k/dmeshz)*dmeshz
      end do
!
      return
      end subroutine calculate_meshift
!
!
!/*---------------------------------------------------------------------
      function sinc (d)
!/*---------------------------------------------------------------------
      implicit none
      real*8   sinc,d,epsi,c2,c4,c6,c8,pi,pid,pid2
!
      epsi =  0.1d0
      c2 = -0.1666666666667d-0
      c4 =  0.8333333333333d-2
      c6 = -0.1984126984127d-3
      c8 =  0.2755731922399d-5
!
      pi = 4.d0*atan(1.d0)
      pid = pi*d
!
      if (dabs(d).gt.epsi) then
         sinc = sin(pid) / pid
      else 
         pid2 = pid*pid
         sinc = 1.d0 + pid2*( c2 + pid2*( c4 + pid2*(c6 + pid2*c8) ) )
      end if
!
      if(dabs(sinc).lt.1.d-100) sinc= 0.d0
!
      return
      end
!
!
!---------------------------------------------------------------------
       function floor (x)
!---------------------------------------------------------------------
       real*8   floor, x, xlim
!
       xlim= 100000.d0
       if(abs(x).lt.xlim) then
         floor = int(x + xlim) - xlim
       else
         write(11,*) " floor: argument too large -- run terminated."
         call exit(1)
       end if
!
       return
       end
!
!
!--------------------------------------------------------------------
      subroutine magnetiz (spx,spy,spz,g,wx,wy,wz,wn,u,uav,wt,np,ic)
!----------------------------------------------------*-***-------**--
!  Magnetization
!
      use, intrinsic :: iso_c_binding
      implicit  none
      include  'param-spinRL5.h'
!
      real(C_DOUBLE) spx(np0),spy(np0),spz(np0),g,wx,wy,wz,wn,u,uav,wt
      integer(C_INT) np,ic,i
!
! ic< 0: Reset and get instantaneous values
! ic= 0: Reset
! ic= 7: Get averaged value
!   one cycle must be ic= 1 -> 7 -> 0
!
      if(ic.le.0) then
        wx= 0
        wy= 0
        wz= 0
        wn= 0
!
        uav= 0
        wt = 0
        if(ic.eq.0) return
      end if
!
      do 100 i= 1,np
      wx= wx +spx(i)
      wy= wy +spy(i)
      wz= wz +spz(i)
      wn= wn +1
  100 continue
!
      uav= uav +u      ! System total energy (not per atom)
      wt = wt +1
!
      if(ic.lt.0 .or. ic.eq.7) then
        wx= -g*wx/wn   ! in mue_B unit
        wy= -g*wy/wn   !  -g*<spx>
        wz= -g*wz/wn
        uav= uav/wt    !  <uav>
      end if
!
      return
      end subroutine magnetiz 
!
!
!------------------------------------------------------------
      subroutine READ_CONF (xleng0,yleng0,zleng0,rank)
!------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include     'param-spinRL5.h'
!
      integer(C_INT) rank,mx,my,mz,itermax,n_MCsteps,nt_P3M
      real(C_DOUBLE) xleng0,yleng0,zleng0,t8,pi,dt,tmax,cptot, &
                     spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,   &
                     tau_b,tau_diss,Temp,TCurie,toler,         &
                     dtwr,dtwr2
      character*4    text1*40 ! param.h -> praefix*6,suffix2*2
!
      real(C_DOUBLE) rad_Fe,rad_O,elj_Fe,elj_O,rcut,rcutC
      common/atoms/  rad_Fe,rad_O,elj_Fe,elj_O,rcut,rcutC
!
      common/parm2/ dtwr,dtwr2
      common/parm3/ t8,pi,dt,tmax,cptot
      common/parm4/ mx,my,mz
      common/parm5/ n_MCsteps,nt_P3M
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b, &
                    tau_diss,Temp,TCurie
      common/itera/ toler,itermax
!
      if(rank.eq.0) then
        OPEN (unit=11,file=praefixc//'.11'//suffix2, &
                status='unknown',position='append',form='formatted')
!                                                            suffix0
        write(11,*) ' read: .',praefixs//'_config.START'//suffix0 
        write(11,*) ' write: ',praefixc//'.13'//suffix2
        write(11,*) 'READ_CONF: Parameter read... start'
!
        close (11)
      end if
!
!
      OPEN (unit=08,file=praefixs//'_config.START'//suffix0,       &
                                     status='old',form='formatted')
!* ------------------------------------------------------
      read (08,'(a40,a6)') text1, praefix   ! String of identification
      read (08,'(a40,f20.0)') text1,cptot   ! Maximum wall time for each run 
      read (08,'(a40,f20.0)') text1,tmax    ! Maximum physical time (ps)
      read (08,'(a40,f20.0)') text1,dt      ! Time step (ps)
      read (08,'(a40,f20.0)') text1,dtwr    ! Energy output interval
      read (08,'(a40,f20.0)') text1,dtwr2   ! Data output interval 
      read (08,'(a40,f20.0)') text1,xleng0  ! Length x of unit cell (Ang)
      read (08,'(a40,f20.0)') text1,yleng0  ! Length y of unit cell 
      read (08,'(a40,f20.0)') text1,zleng0  ! Length z of unit cell 
!
      if(rank.eq.0) then 
        OPEN (unit=11,file=praefixc//'.11'//suffix2, &
                status='unknown',position='append',form='formatted')
!
        write(11,*) '### READ_CONF ###'
        write(11,*) 
!
        write(11,600) cptot,tmax,dt,dtwr,dtwr2,xleng0,yleng0,zleng0
  600   format('cptot,tmax=',2f8.1,/, &
               'dt,dtwr,dtwr2=',1p3d11.1,/, &
               'xleng0,yleng0,zleng0=',0p3f8.1)
!
        close (11)
      end if 
!        
!     read (08,'(a40,i10)')   text1,mx      ! Number of domains
!     read (08,'(a40,i10)')   text1,my      !    <--param.h, meshx= 5
!     read (08,'(a40,i10)')   text1,mz      ! 
      read (08,'(a40,i10)')   text1,n_MCsteps ! Maximum iteration count of MC
      read (08,'(a40,f20.0)') text1,rcut    ! Cutoff radius of spin interaction
      read (08,'(a40,f20.0)') text1,rcutC   ! Cutoff radius of ES interaction
      read (08,'(a40,i10)')   text1,itermax ! Maximum iteration of nonlinear solver 
      read (08,'(a40,f20.0)') text1,toler   ! Convergence or iteration tolerance
      read (08,'(a40,f20.0)') text1,Jaa     ! Exchange integral A-A (10meV)
      read (08,'(a40,f20.0)') text1,Jbb     ! Exchange integral B-B 
      read (08,'(a40,f20.0)') text1,Jab     ! Exchange integral A-B 
!
      if(rank.eq.0) then
        OPEN (unit=11,file=praefixc//'.11'//suffix2, &
                status='unknown',position='append',form='formatted')
!
!       write(11,602) mx,my,mz,n_MCsteps,rcut,rcutC, &
        write(11,602) n_MCsteps,rcut,rcutC,   &
                      itermax,toler,Jaa,Jbb,Jab
  602   format('n_MCsteps=',i8,/, &
               'rcut,rcutCi=',2f8.1,/, &
               'itermax,toler=',i8,1pd11.3,/, & 
               'Jaa,Jbb,Jab=',3d11.3)
!
        close (11)
      end if  
!
      read (08,'(a40,f20.0)') text1,Bapx    ! Magnetic field Bx (100gauss)
      read (08,'(a40,f20.0)') text1,Bapy    ! Magnetic field By
      read (08,'(a40,f20.0)') text1,Bapz    ! Magnetic field Bz
      read (08,'(a40,f20.0)') text1,tau_b   ! Period of microwave (ps)
      read (08,'(a40,f20.0)') text1,tau_diss! Dissipation relaxation time (ps)
      read (08,'(a40,f20.0)') text1,Temp    ! Temperature in Kelvin
      read (08,'(a40,f20.0)') text1,TCurie  ! Curie Temperature in Kelvin
!
      if(rank.eq.0) then
        OPEN (unit=11,file=praefixc//'.11'//suffix2, &
                status='unknown',position='append',form='formatted')
!
        write(11,604) Bapx,Bapy,Bapz,tau_b,tau_diss,Temp,TCurie
  604   format('Bapx,Bapy,Bapz=',1p3d11.3,/, &
               'tau_b,tau_diss=',2d11.3,/, &
               'Temp,TCurie=',0p2f8.1)
!
        close (11)
      end if
!  
      return
      end subroutine READ_CONF
!
!
!-----------------------------------------------------------------------
      subroutine Temp_fix (spx,spy,spz,sp2,atsz,tsz0,spec,site,np, &
                           ifcmp,nstep_MCt)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'mpif.h'
      include    'param-spinRL5.h'
!
      real(C_DOUBLE) spx(np0),spy(np0),spz(np0),sp2(np0), &
                 atsz,tsz0,ss,sp,sc,alpha,pi,th,ph,ranff,th0, &
                 sa3,sb3,sb2
      integer(C_INT) spec(np0),site(np0),np,ifcmp,nstep_MCt, &
                 i_MC,i,j,na3,nb3,nb2
!
! ---------------------------------
      if(ifcmp.eq.1) go to 2000
! ---------------------------------
!
!* Radomization...
!
      i_MC= 0
      pi= 4.d0*atan(1.d0)
!
 1000 i_MC= i_MC +1
      if(i_MC.gt.nstep_MCt) return   ! No operation if nstep_MCt= 0
!
      j= mod(int(np*ranff(0.d0)+0.00001),np) +1
      th=   pi*ranff(0.d0)
      ph= 2*pi*ranff(0.d0)
!
      th0= dacos(spz(j)/sp2(j))
      th = 0.7d0*th0 +0.3d0*th    ! weighted average 
! 
!     spx(j)= sp2(j)*sin(th)*cos(ph)
!     spy(j)= sp2(j)*sin(th)*sin(ph)
      spz(j)= sp2(j)*cos(th)
!
      ss= sqrt(sp2(j)**2 -spz(j)**2)
      sp= sqrt(spx(j)**2 +spy(j)**2)
!
      spx(j)= ss *spx(j)/sp       ! keep the rotation phase
      spy(j)= ss *spy(j)/sp
!
      go to 1000
! ----------------------------
!
 2000 continue
!
!* Compression toward its axis ...
!
      sa3= 0     ! A site,3+
      sb3= 0     ! B site,3+
      sb2= 0     ! B site,2+
      na3= 0
      nb3= 0
      nb2= 0
!
      do 300 i= 1,np
      if(site(i).eq.1) then
        sa3= sa3 +spz(i)
        na3= na3 +1
      else
        if(spec(i).eq.1) then
          sb3= sb3 +spz(i)
          nb3= nb3 +1
        else
          sb2= sb2 +spz(i)
          nb2= nb2 +1
        end if
      end if
  300 continue
!
      sa3= sa3/na3   ! Axis of distribution
      sb3= sb3/nb3
      sb2= sb2/nb2
!
!* Compress the length  -s.......0.......s
!                        i     s         i
!               -s+alp*(s+s)         s-alp*(s-s)
!               =(-s)*(1-alp)+alp*s  =s*(1-alp)+alp*s
!
      alpha= tsz0/atsz
!
      do 400 i= 1,np
      if(site(i).eq.1) then
        if(sa3.lt.0.) sc= -sp2(i)
        if(sa3.gt.0.) sc=  sp2(i)
      else
        if(spec(i).eq.1) then
          if(sb3.lt.0.) sc= -sp2(i)
          if(sb3.gt.0.) sc=  sp2(i)
        else
          if(sb2.lt.0.) sc= -sp2(i)
          if(sb2.gt.0.) sc=  sp2(i)
        end if
      end if
!
      spz(i)= sc*(1.d0 -alpha) + alpha*spz(i)
!
      ss= sqrt(sp2(i)**2 -spz(i)**2)
      sp= sqrt(spx(i)**2 +spy(i)**2)
!
      spx(i)= ss *spx(i)/sp
      spy(i)= ss *spy(i)/sp
  400 continue
!
      return
      end subroutine Temp_fix 
!
!
!-------------------------------------------------------------
      subroutine plot_spin (x,y,z,spx,spy,spz,spec,site,np)
!-------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      include  'param-spinRL5.h'
!
      real(C_DOUBLE) x(np0),y(np0),z(np0),spx(np0),spy(np0),spz(np0)
      integer(C_INT) spec(np0),site(np0)
      real(C_DOUBLE) t8,xleng,yleng,zleng,pi,dt,tmax,cptot, &
                 spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b,tau_diss,&
                 Temp,TCurie
      common/parm3/ t8,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b, &
                    tau_diss,Temp,TCurie
!
      character*8    label,cdate*10,ctime*8,cax*1
      common/headr1/ label,cdate
      common/headr2/ time
!
      hh= 0.7
      call symbol ( 0.5,17.0,hh,label,0.,8)
      call symbol ( 5.0,17.0,hh,praefixc,0.,29)
      call symbol (15.0,17.0,hh,'Spin dynamics (MD)',0.,18)
!
      call symbol ( 0.0,1.0,hh,cdate, 0.,10)
      call symbol (19.0,1.0,hh,'time=', 0.,5)
      call values (20.5,1.0,hh,time,0.,101)
!
      spin24= spin2
      spin34= spin3
      TCurie4= TCurie
      call symbol (0.5,16.0,hh,'spin(A)=', 0.,8)
      call values (3.5,16.0,hh,spin24,0.,101)
      call symbol (0.5,15.2,hh,'spin(B)=', 0.,8)
      call values (3.5,15.2,hh,spin34,0.,101)
      call symbol (0.5,14.4,hh,'TCurie= ', 0.,8)
      call values (3.5,14.4,hh,TCurie4,0.,101)
!
      qaa4= Jaa  ! in 1 MeV - 10 MeV 
      qbb4= Jbb
      qab4= Jab
      call symbol (7.5,16.0,hh,'Jaa=', 0.,4)
      call values (9.0,16.0,hh,qaa4,0.,101)
      call symbol (7.5,15.2,hh,'Jbb=', 0.,4)
      call values (9.0,15.2,hh,qbb4,0.,101)
      call symbol (7.5,14.4,hh,'Jab=', 0.,4)
      call values (9.0,14.4,hh,qab4,0.,101)
!
      Temp4= Temp
      Bext4= sqrt(Bapx**2 +Bapy**2 +Bapz**2)  ! 100 gauss
      call symbol (13.0,16.0,hh,'Temp=', 0.,5)
      call values (15.0,16.0,hh,Temp4,0.,101)
      call symbol (13.0,15.2,hh,'Bapp=', 0.,5)
      call values (15.0,15.2,hh,Bext4,0.,101)
!
!     xleng4= xleng4
      dt4 = dt
      tau4= tau_b    !<-- 400 ps
      tau_diss4= tau_diss
!     call symbol (19.0,16.0,hh,'Leng=', 0.,5)
!     call values (21.0,16.0,hh,xleng4,0.,101)
      call symbol (19.0,16.0,hh,'tau_dis=', 0.,8)
      call values (21.0,16.0,hh,tau_diss4,0.,101)
      call symbol (19.0,15.2,hh,'np=', 0.,3)
      call values (21.0,15.2,hh,float(np),0.,101)
      call symbol (19.0,14.4,hh,'dt=', 0.,3)
      call values (21.0,14.4,hh,dt4,0.,101)
      call symbol (19.0,13.6,hh,'tau_b=', 0.,6)
      call values (21.0,13.6,hh,tau4,0.,101)
!
      fsize= 4.
      hl=  11.
      vd=   6.
      phi= -60.
      tht=  15.
!
      pi=  4.e0*atan(1.0e0)
      pha= pi*phi/180.
      tha= pi*tht/180.
!
      cph= cos(pha)
      sph= sin(pha)
      cth= cos(tha)
      sth= sin(tha)
!
      xmax= 0.5*xleng
      ymax= 0.5*yleng
      zmax= 0.5*zleng
!
      xp= xmax*cph -ymax*sph
      yp= xmax*sph +ymax*cph
      zp= zmax
!
      ypp=  yp
      zpp= -xp*sth +zp*cth
!
      rmax1= sqrt(ypp**2 +zpp**2)
      ps= fsize/rmax1
!
!**********************
!*  Draw axes.        *
!**********************
!
      do 100 i= 1,3
      if(i.eq.1) then
         x1= xmax
         y1= 0.
         z1= 0.
         cax='x'
      else if(i.eq.2) then
         x1= 0.
         y1= ymax
         z1= 0.
         cax='y'
      else if(i.eq.3) then
         x1= 0.
         y1= 0.
         z1= zmax
         cax='z'
      end if
!
      xp= x1*cph -y1*sph
      yp= x1*sph +y1*cph
      zp= z1
!
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp  +hl
      yy= ps*zpp  +vd
      call plot (hl,vd,3)
      call plot (xx,yy,2)
!
      call symbol (xx-0.7,yy-0.5,hh,cax,0.,1)
  100 continue
!
!-------------------
      dd= 0.2
      ff= 0.25
!-------------------
      do 400 i= 1,np
      xp= (x(i) -ff*spx(i))*cph -(y(i) -ff*spy(i))*sph
      yp= (x(i) -ff*spx(i))*sph +(y(i) -ff*spy(i))*cph
      zp=  z(i) -ff*spz(i) 
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
      xx1= ps*ypp +hl
      yy1= ps*zpp +vd
      if(xx1.lt.0. .or. xx1.gt.20.) go to 400
      if(yy1.lt.0. .or. yy1.gt.18.) go to 400
!
      xp= (x(i) +ff*spx(i))*cph -(y(i) +ff*spy(i))*sph
      yp= (x(i) +ff*spx(i))*sph +(y(i) +ff*spy(i))*cph
      zp=  z(i) +ff*spz(i) 
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
      xx2= ps*ypp +hl
      yy2= ps*zpp +vd
!
      if(spec(i).eq.1) then      !* Fe(3+)
        if(site(i).eq.1) then    ! A-site
          call newcolor (3,0.,1.,0.)  ! 3-color (r,g,b)
          call circle (xx1-0.01,yy1-0.05,dd,2)
        else                     ! B-site
          call newcolor (3,0.,0.,1.)  ! 3-color (r,g,b)
          call circle (xx1-0.01,yy1-0.05,dd,2)
        end if
      end if
      if(spec(i).eq.2) then     !* Fe(2+)
        call newcolor (3,1.,0.,0.)  ! 3-color (r,g,b)
        call triang (xx1-0.01,yy1-0.05,dd,2)
      end if
!
      call plot (xx1,yy1,3)
      call plot (xx2,yy2,2)
      call plot (xx2,yy2,3)
  400 continue
!
      call newcolor (0,1.,0.,0.)  ! gray scale
!---------------------
      call chart
!---------------------
!
      return
      end subroutine plot_spin 
!
!
!-----------------------------------------------------------
      subroutine distr_spin (spx,spy,spz,spec,site,np1)
!-----------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include  'param-spinRL5.h'
!
      integer(C_INT) spec(np0),site(np0),np1,ir,ILN,ILG,i,k,      &
                     nxtick,nytick
      real(C_DOUBLE) spx(np0),spy(np0),spz(np0),ss
      real(C_float)  ang1x(101),ang2x(101),ang1y(101),ang2y(101), &
                     ang1z(101),ang2z(101),angax(101),angay(101), &
                     angaz(101),hha(101),costhx,costhy,costhz,    &
                     famax1,famin1,famax2,famin2,famax3,famin3,   &
                     famax,famin
!
      do 100 k= 1,101
      ang1x(k)= 0
      ang2x(k)= 0
      ang1y(k)= 0
      ang2y(k)= 0
      ang1z(k)= 0
      ang2z(k)= 0
      angax(k)= 0
      angay(k)= 0
      angaz(k)= 0
      hha(k)= (k-51)/50.e0
  100 continue
!
      do 200 i= 1,np1
      ss= sqrt(spx(i)**2+spy(i)**2+spz(i)**2)
!
      costhx= spx(i)/ss
      costhy= spy(i)/ss
      costhz= spz(i)/ss
!
      ir= 50.*(costhx +1.) +1.01
      if(ir.gt.0 .and. ir.lt.101) then
        if(site(i).eq.1) then
          angax(ir)= angax(ir) +1.                    ! A:3+
        else
          if(spec(i).eq.1) ang1x(ir)= ang1x(ir) +1.   ! B:3+
          if(spec(i).eq.2) ang2x(ir)= ang2x(ir) +1.   ! B:2+
        end if
      end if
!
      ir= 50.*(costhy +1.) +1.01
      if(ir.gt.0 .and. ir.lt.101) then
        if(site(i).eq.1) then
          angay(ir)= angay(ir) +1.
        else
          if(spec(i).eq.1) ang1y(ir)= ang1y(ir) +1.
          if(spec(i).eq.2) ang2y(ir)= ang2y(ir) +1.
        end if
      end if
!
      ir= 50.*(costhz +1.) +1.01
      if(ir.gt.0 .and. ir.lt.101) then
        if(site(i).eq.1) then
          angaz(ir)= angaz(ir) +1.
        else
          if(spec(i).eq.1) ang1z(ir)= ang1z(ir) +1.
          if(spec(i).eq.2) ang2z(ir)= ang2z(ir) +1.
        end if
      end if
  200 continue
!
!                                                                *
      ILN=  1
      ILG=  2
      call lplmax (ang1x,famax1,famin1,101)
      call lplmax (ang2x,famax2,famin2,101)
      famax= amax1(famax1,famax2)
      famin= amin1(famin1,famin2)
!
        OPEN (unit=11,file=praefixc//'.11'//suffix2, &
                status='unknown',position='append',form='formatted')
!
        write(11,932) famax1,famax2,famin1,famin2
  932   format('max,min=',1p4e11.2)
!
        close (11)
!
      nxtick= 3
      nytick= 3
      call hplot1 (2,2,101,hha,ang1x,famax,famin,ILN,nxtick,nytick,&
                     '        ',8,'cos(thx)',8,' (3 B   ',8,0)
      call hplot1 (2,2,101,hha,angax,famax,famin,ILN,nxtick,nytick,&
                     '        ',8,'cos(thx)',8,' (3A    ',8,1)
      call hplot1 (2,3,101,hha,ang2x,famax,famin,ILN,nxtick,nytick,&
                     '        ',8,'cos(thx)',8,' (2B)   ',8,0)
!
      call lplmax (ang1y,famax1,famin,101)
      call lplmax (ang2y,famax2,famin,101)
      famax= amax1(famax1,famax2)
      call hplot1 (3,2,101,hha,ang1y,famax,famin,ILN,nxtick,nytick,&
                     '        ',8,'cos(thy)',8,' (3 B   ',8,0)
      call hplot1 (3,2,101,hha,angay,famax,famin,ILN,nxtick,nytick,&
                     '        ',8,'cos(thy)',8,' (3A    ',8,1)
      call hplot1 (3,3,101,hha,ang2y,famax,famin,ILN,nxtick,nytick,&
                     '        ',8,'cos(thy)',8,' (2B)   ',8,0)
!---------------------
      call chart
!---------------------
!
      call lplmax (ang1z,famax1,famin,101)
      call lplmax (ang2z,famax2,famin,101)
      famax= amax1(famax1,famax2)
      call hplot1 (2,2,101,hha,ang1z,famax,famin,ILN,nxtick,nytick,&
                     '        ',8,'cos(thz)',8,' (3 B   ',8,0)
      call hplot1 (2,2,101,hha,angaz,famax,famin,ILN,nxtick,nytick,&
                     '        ',8,'cos(thz)',8,' (3A    ',8,1)
      call hplot1 (2,3,101,hha,ang2z,famax,famin,ILN,nxtick,nytick,&
                     '        ',8,'cos(thz)',8,' (2B)   ',8,0)
!---------------------
      call chart
!---------------------
      return
      end subroutine distr_spin 
!
!
!------------------------------------------------------
      subroutine plot_disp (axis,freq,modes,plot_ch)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      integer(C_INT) modes,ILN,nxtick,nytick
      real(C_float) axis(100),freq(100),emax,emin
      character*8   plot_ch
!
      ILN= 1
      nxtick= 4
      nytick= 4
      call lplmax (freq,emax,emin,modes)
        write(11,*) 'plot_disp: max,min=',emax,emin
!
      call lplot1 (1,1,modes,axis,freq,emax,emin,ILN,nxtick, &
                   nytick,plot_ch,8,'wave num',8,'        ',8,0)
!---------------------
      call chart
!---------------------
      return
      end subroutine plot_disp 
!
!
!------------------------------------------------------
      subroutine lplots 
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit      none
      include      'param-spinRL5.h'
!
      real(C_float) spinx,spinz,spin7,Bextx,Bextz,magx,magy,magz, &
                    Usys,conv,aitr,psdt,tfix,uss,usb,tsx,tsy,tsz, &
                    sum_mb,U_Fe,U_O,ds_Fe,ds_O,fdt4,vdt4,idt4,timeh,&
                    cosd,sind,B00,Bmw,mue_B,hh,av_mz(nhs),ss,s2,  &
                    emax,emin,emax1,emax2,emax3,emin1,emin2,emin3
      common/ehist/ spinx(nhs),spinz(nhs),spin7(nhs), &
                    Bextx(nhs),Bextz(nhs),magx(nhs),magy(nhs),magz(nhs),&
                    Usys(nhs),conv(nhs),aitr(nhs),psdt(nhs),tfix(nhs),  &
                    uss(nhs),usb(nhs),tsx(nhs,3),tsy(nhs,3),tsz(nhs,3), &
                    sum_mb(nhs),U_Fe(nhs),U_O(nhs),ds_Fe(nhs),ds_O(nhs),&
                    fdt4(nhs),vdt4(nhs),idt4(nhs),timeh(nhs)
!
      real(C_DOUBLE) spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b, &
                    tau_diss,Temp,TCurie
      common/spins/ spin2,spin3,Jaa,Jbb,Jab,Bapx,Bapy,Bapz,tau_b, &
                    tau_diss,Temp,TCurie
      real(C_DOUBLE) t8,xleng,yleng,zleng,pi,dt,tmax,cptot
      common/parm3/ t8,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
!
      character*8    label,cdate*10,ctime*8
      integer(C_INT) i,k,it,is,is0,ns,ILN,ILG,nxtick,nytick
      common/headr1/ label,cdate
      common/parm1/  it,is
!
      ILN= 1
      ILG= 2
      hh= 0.7
      call symbol (1.0,18.5,hh,'Spin dynamics (MD)',0.,18)
!
!-------------------------------------------------
!* One-period averaged mz(t)= mz - <mz>
!   40 data points are generated per period
!-------------------------------------------------
!
      is0= av_start*is  ! Start averaging at this step
!
      do 103 k= 1,is
      av_mz(k)= magz(k)
  103 continue
!
      do 100 k= 40,is-40
      ss= 0
      ns= 0
!
      do 130 i= k-39,k+40
      ss= ss +magz(i)
      ns= ns +1
  130 continue
!
      av_mz(k)= ss/ns
  100 continue
!
!
      ss= 0
      s2= 0
      ns= 0
!
      do 200 k= is0,is-40
      ss= ss +(magz(k) -av_mz(k))*Bextz(k)   ! magz= -g*<sz>
      s2= s2 +(magz(k) -av_mz(k))**2         ! <dm**2>
      ns= ns +1
  200 continue
!
      mue_B= 9.27410e-21   ! e*hbar/2mc
      B00= 100             ! gauss
      Bmw= B00*sqrt(Bapx**2 +Bapy**2 +Bapz**2)
!
      call lplmax (magx,emax1,emin1,is)
      call lplmax (magy,emax2,emin2,is)
      call lplmax (magz,emax3,emin3,is)
      emax = amax1(emax1,emax2,emax3,-emin1,-emin2,-emin3)
!
!
      call lplmax (magx,emax1,emin1,is)
      call lplmax (magy,emax2,emin2,is)
      call lplmax (magz,emax3,emin3,is)
      emax = amax1(emax1,emax2,emax3,-emin1,-emin2,-emin3)
      emin = -emax
      nxtick= 3
      nytick= 3
      call lplot1 (2,4,is,timeh,magx,emax,emin,ILN,nxtick,nytick,&
                   'magnet-x',8,'        ',8,'        ',8,0)
      call lplot1 (2,5,is,timeh,magy,emax,emin,ILN,nxtick,nytick,&
                   'magnet-y',8,'        ',8,'        ',8,0)
      call lplot1 (2,6,is,timeh,magz,emax3,emin3,ILN,nxtick,nytick,&
                   'magnet-z',8,' time   ',8,'        ',8,0)
!
!     call lplmax (spinx,emax1,emin1,is)
      call lplmax (spinz,emax2,emin2,is)
      emax = amax1(emax1,emax2,-emin1,-emin2)
      emin = -emax
!     call lplot1 (3,4,is,timeh,spinx,emax,emin,ILN,nxtick,nytick,&
!                  'spin-x.7',8,'        ',8,'        ',8,0)
      call lplot1 (3,5,is,timeh,spinz,emax,emin,ILN,nxtick,nytick,&
                   'spin-z.7',8,'        ',8,'        ',8,0)
!
      call lplmax (spin7,emax,emin,is)
      emax= 1.2*emax
      call lplot1 (3,6,is,timeh,spin7,emax,0.,ILN,nxtick,nytick,&
                   'spin-7  ',8,' time   ',8,'        ',8,0)
!------------------------
      call chart
!------------------------
!
      call lplmax (Usys,emax,emin,is)
      call lplot1 (2,4,is,timeh,Usys,emax,emin,ILN,nxtick,nytick,&
                   'Usys    ',8,'        ',8,'        ',8,0)
!
      call lplmax (Bextx,emax,emin,is)
      emax = amax1(emax,-emin)
      emin = -emax
      call lplot1 (2,5,is,timeh,Bextx,emax,emin,ILN,nxtick,nytick,&
                   'Bextx.7 ',8,'        ',8,'        ',8,0)
!
      call lplmax (Bextz,emax,emin,is)
      emax = amax1(emax,-emin)
      emin = -emax
      call lplot1 (2,6,is,timeh,Bextz,emax,emin,ILN,nxtick,nytick,&
                   'Bextz.7 ',8,'  time  ',8,'        ',8,0)
!
      call lplmax (uss,emax,emin,is)
      call lplot1 (3,4,is,timeh,uss,emax,emin,ILN,nxtick,nytick,&
                   ' us*s   ',8,'        ',8,'        ',8,0)
!
      call lplmax (usb,emax,emin,is)
      call lplot1 (3,5,is,timeh,usb,emax,emin,ILN,nxtick,nytick,&
                   ' us*b   ',8,'        ',8,'        ',8,0)
!------------------------
      call chart
!------------------------
!
      call lplmax (ds_Fe,emax,emin,is)
      emax = amax1(emax,-emin)
      emin = -emax
      call lplot1 (2,4,is,timeh,ds_Fe,emax,emin,ILN,nxtick,nytick,&
                   ' ds_Fe  ',8,'        ',8,'        ',8,0)
!
      call lplmax (ds_O,emax,emin,is)
      emax = amax1(emax,-emin)
      emin = -emax
      call lplot1 (2,5,is,timeh,ds_O,emax,emin,ILN,nxtick,nytick,&
                   ' ds_O   ',8,'        ',8,'        ',8,0)
!
      call lplmax (U_O,emax,emin,is)
      call lplot1 (3,4,is,timeh,U_O,emax,emin,ILN,nxtick,nytick,&
                   ' Kin_O  ',8,'        ',8,'        ',8,0)
!
      call lplmax (U_Fe,emax,emin,is)
      call lplot1 (3,5,is,timeh,U_Fe,emax,emin,ILN,nxtick,nytick,&
                   ' Kin_Fe ',8,'        ',8,'        ',8,0)
!
      call lplmax (fdt4,emax,emin,is)
      call lplot1 (2,6,is,timeh,fdt4,emax,emin,ILN,nxtick,nytick,&
                   'f*dt/mvo',8,'  time  ',8,'        ',8,0)
!
      call lplmax (vdt4,emax,emin,is)
      call lplot1 (3,6,is,timeh,vdt4,emax,emin,ILN,nxtick,nytick,&
                   ' v*dt   ',8,'  time  ',8,'        ',8,0)
!------------------------
      call chart
!------------------------
!
      return
      end subroutine lplots 
!
!
!------------------------------------------------------
      subroutine rehist (rank)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit      none
      include       'param-spinRL5.h'
!
      integer(C_INT) rank,i,j,it,is,iss,iwa,iwb
      common/parm1/ it,is
      common/imemo/ iwa,iwb
!
      real(C_float) array
      common/ehist/ array(nhs,33)
!
      real(C_DOUBLE) dtwr,dtwr2
      common/parm2/  dtwr,dtwr2
!
      real(C_DOUBLE) t8,xleng,yleng,zleng,pi,dt,tmax,cptot
      common/parm3/ t8,pi,dt,tmax,cptot
      common/ewald3/ xleng,yleng,zleng
!     character     suffix2*2
!
      do 100 j= 1,33
      iss= 0
!
      do 100 i= 1,is,2
      iss= iss +1
      array(iss,j)= 0.5*(array(i,j) +array(i+1,j))
  100 continue
!
      is= iss
      dtwr= 2*dtwr
      iwa= 0 ! t/dtwr
!
      if(rank.eq.0) then
        write(11,*) ' rehist is called: it, is=',it,is
      end if
!
      return
      end subroutine rehist 
!
!
!---------------------------------------------------------------
      subroutine ggauss
!---------------------------------------------------------------
!  (-3.,3.) case
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) fv,vv0,dv,vv,s,sdv,fun
      integer(C_INT) j,ns,k2,k,i
      common/gaus1/ fv(51),vv0,dv
!
      fv(1)= 0
!
      vv0= -3.d0
      dv= 2.d0*abs(vv0)/50.d0
!
      vv= vv0
      do 100 j=1,50
      s= 0
      ns=1000
      k2=ns/2
      sdv=dv/ns
!
      do 130 k=1,k2
      vv=vv +2.d0*sdv
      s=s +4.d0*fun(vv-sdv) +2.d0*fun(vv)
  130 continue
      s= (s +4.d0*fun(vv+sdv) +fun(vv+2.d0*sdv))*sdv/3.d0
      fv(j+1)= fv(j)+s
  100 continue
!
      do 200 i=1,51
      fv(i)= fv(i)/fv(51)
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
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) dgaus2,vmax,fv,vv0,dv,s,sdv,fun,eps,x2,y1,y2
      real(C_DOUBLE) ranff
      integer(C_INT) j,ns,k2,k
!
      common/gaus1/ fv(51),vv0,dv
!
      eps= ranff(0.d0)
!
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
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) fun,v
!
      fun= exp(-v**2/2.)
!
      return
      end function fun
!
!
!------------------------------------------------------
      subroutine averg1 (q,qav,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      integer(C_INT) is
      real(C_float) q(is),qav
!
      qav= 0.
!
      do 100 i= is-9,is
      qav= qav +q(i)
  100 continue
!
      qav= qav/10.
!
      return
      end subroutine averg1 
!
!
!------------------------------------------------------
      subroutine clocks (cputime,walltime)
!------------------------------------------------------
!* get the cpu and elapsed times on each node.
!
      use, intrinsic :: iso_c_binding
      include 'mpif.h'
!
      real(C_float) ct,tm(2)
      real(C_DOUBLE) cputime,walltime,walltime0
      logical first_clk
      common/saveif2/ walltime0,first_clk
!               - defined as .true. in /block data/
!
      if(first_clk) then
         walltime0= mpi_wtime()
         first_clk= .false.
      end if
!
!     ct = etime(tm)   ! pentium and unix general.
      cputime = 0      ! tm(1)
      walltime= mpi_wtime() -walltime0  ! sec
!
      return
      end subroutine clocks 
!
!
!------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      include    'param-spinRL5.h'
!
      integer(C_INT) is
      real(C_float) f(is),fmax,fmin
!
      fmax= -1.e10
      fmin=  1.e10
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
!-------------------------------------------------
      subroutine circle (x,y,d,ic)
!-------------------------------------------------
!*  open circle centered at (x,y) /or outer edge.
!
      use, intrinsic :: iso_c_binding
      integer(C_INT) ic
      real(C_float) x,y,d
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
!*  filled circle centered at (x,y).
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
      end subroutine circle 
!
!
!-------------------------------------------------
      subroutine triang (x,y,d,ic)
!-------------------------------------------------
!*  open triangle centered at (x,y) /or outer edge.
!
      use, intrinsic :: iso_c_binding
      integer(C_INT) ic
      real(C_float) x,y,d
!
      write(77,*) " 3.0 setlinewidth"
!
      a= d/2.
      b= a/1.732
      c= 2.*b
!
      call plot (x-a,y-b,3)
      call plot (  x,y+c,2)
      call plot (x+a,y-b,2)
      call plot (x-a,y-b,2)
      call plot (x+a,y-b,3)
!
      write(77,*) " 1.0 setlinewidth"
      if(ic.eq.1) return
!
!------------------------
!*  fill the triangle.
!------------------------
      nc=7
      nch= (nc+1)/2
      dx= d/nc
      y2= y -b
!
      do 100 j= 1,nc
      x1= x-a +dx*(j-0.5)
      if(j.le.nch) then
         y1= y +1.732*a*(j-0.5)/nch
      else
         y1= y +1.732*a*(nc-j)/nch 
      end if
!
      call plot (x1,y1,3)
      call plot (x1,y2,2)
  100 continue
!
      call plot (x1,y2,3)
      write(77,*) " 1.0 setlinewidth"
!
      return
      end subroutine triang 
!
!
!------------------------------------------------
      function iwrta (t8,dtwr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding
!
      integer(C_INT) iwrta,iw
      real*8  t8,dtwr
      common/imemo/ iwa,iwb
!
      iw= t8/dtwr
      if(iw.gt.iwa) then
        iwa= iw
        iwrta= 0
      else
        iwrta= 1
      end if
!
      return
      end function iwrta 
!
!
!------------------------------------------------
      function iwrtb (t8,dtwr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding
!
      integer(C_INT) iwrtb,iw
      real*8  t8,dtwr
      common/imemo/ iwa,iwb
!
      iw= t8/dtwr 
      if(iw.gt.iwb) then
        iwb= iw
        iwrtb= 0
      else
        iwrtb= 1
      end if
!
      return
      end function iwrtb 
!
!
!------------------------------------------------
      block data
!------------------------------------------------
      real*8  walltime0
      integer*4 ir0
      logical first_shf,first_clk
!
      common/saveif2/ walltime0,first_shf,first_clk
      common/ranfff/ ir0
!
      data  first_shf/.true./, first_clk/.true./
      data  ir0/17331/
      end
!
!
!------------------------------------------------
      function ranff (x)
!------------------------------------------------
!*  ranf= (0,1)
!
      use, intrinsic :: iso_c_binding
      real(C_DOUBLE) ranff,x
      integer(C_INT) ir0
      common/ranfff/ ir0
!
      real*8     invm
      parameter  (mask=2**30+(2**30-1),invm= 0.5d0**31)
      parameter  (lambda=48828125)
!
      ir0= iand( lambda*ir0, mask)
      ranff= ir0*invm
!
!     ask= 371597.
!     ambda= sqrt(ask)
!     qq= 0.3713*ask
!
!     ir0= amod( ambda*ir0 +qq, ask)
!     ranff= ir0/ask
!
      return
      end function ranff 
!
!
!-----------------------------------------------------------------------
      subroutine lplot1 (ix,iy,npt1,x,y,ymax,ymin,IL,nxtick, &
                         nytick,lab1,n1,lab2,n2,lab3,n3,iskip)
!-----------------------------------------------------------------------
!  <<warning>>  order and number of arguments /lplot/ have been changed.
!               also, x (time) is defined for all range.
!               date: 5/18/96 at mit.
!***********************************************************************
!   il=1................ linear plot of (x,y)
!   il=2................ log10 plot of (x,log y)
!***********************************************************************
!
      use, intrinsic :: iso_c_binding
!
      integer(C_INT) ix,iy,npt1,IL,nxtick,nytick,n1,n2,n3,iplot
      real(C_float)  x(npt1),y(npt1),u(npt1),v(npt1),   &
                     xmax,xmin,ymax,ymin,time
      dimension  xcm(6),ycm(7),pl(6),pr(6),ql(7),qr(7)
!
      character*8    lab1,lab2,lab3,label,cdate*10
      common/headr1/ label,cdate
      common/headr2/ time
      common/pplcom/ nfine,pl1(10),pr1(10),ql1(10),qr1(10), &
                     xmin1(10),xmax1(10),ymin1(10),ymax1(10)
!
      data  xcm/21.0, 2*10.00, 3*6.00/,       &
            ycm/15.0, 2*6.80, 4*4.00/,        &
            pl/2.0,  2.0,14.0, 2.0,9.0,16.0/, &
            ql/2.3, 10.5,2.3, 14.0,9.5,5.0,0.5/
!
      call set_width (1.0)
!
      iplot=1
      go to 1
!
!-----------------------------------------------------------------------
      entry hplot1 (ix,iy,npt1,x,y,ymax,ymin,IL,nxtick,nytick,  &
                    lab1,n1,lab2,n2,lab3,n3,iskip)
!-----------------------------------------------------------------------
      iplot=2
!
    1 npt= npt1
      isc= 1
!
      do 5 i=1,6
    5 pr(i)= pl(i) +xcm(i)
!
      do 6 j=1,7
    6 qr(j)= ql(j) +ycm(j)
!
!                 ******************************************************
!*                **  Make a copy before the top-left frame is drawn. **
!                 ******************************************************
      hh= 0.70
      i1= iabs(ix)
      j1= iabs(iy)
      if(i1.ge.3) go to 10
      if(j1.eq.3.or.j1.ge.5) go to 10
!                                              ************************
!                                              ** label of the page. **
!                                              ************************
      call symbol (20.3,0.1,hh,'t=',0.,2)
      call values (21.5,0.1,hh,time,0.,101)
!
   10 continue
!
      do 23 i=1,npt
   23 u(i)= x(i)
      xmax= u(npt)
      xmin= u(1)
!                             ************************************
!                             ** three-point average if il > 0  **
!                             ************************************
      if(il.gt.0) then
        v(1)=   y(1)
        v(npt)= y(npt)
        v(npt-1)= y(npt-1)
!
        do 37 i=2,npt-2
   37   v(i)= 0.33333*(y(i-1)+y(i)+y(i+1))
      else
        do 38 i=1,npt
   38   v(i)= y(i)
      end if
!                                                *****************
!                                                **  log. scale **
!                                                *****************
      if(iabs(il).eq.2) then
         do 40 i=1,npt
         if(v(i).gt.0.) then
            v(i)= alog10(v(i))
         else
            v(i)= -10.
         end if
   40    continue
      end if
!                                **************************************
!                                ** Set a new scale and draw a frame.**
!                                **************************************
      if(iplot.eq.2) then
         ymax= -1.e10
         ymin=  1.e10
!
         do 50 i= 1,npt
         ymax= amax1(ymax,v(i))
         ymin= amin1(ymin,v(i))
   50    continue
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
!                                                      **  Frame. **
!                                                      *************
      call plot (pl(i1),ql(j1),3)
      call plot (pl(i1),qr(j1),2)
      call plot (pr(i1),qr(j1),2)
      call plot (pr(i1),ql(j1),2)
      call plot (pl(i1),ql(j1),2)
!                                                    ******************
!                                                    **  Tick marks. **
!                                                    ******************
      scx= xcm(i1)/(nxtick+1)
      scy= ycm(j1)/(nytick+1)
!
      x0= pl(i1)
      y1= ql(j1)
      y4= qr(j1)
      y2= y1 +0.25
      y3= y4 -0.25
!
      do 62 k=1,nxtick
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
      do 63 k=1,nytick
      y0= y0 +scy
      call plot (x1,y0,3)
      call plot (x2,y0,2)
      call plot (x3,y0,3)
      call plot (x4,y0,2)
   63 continue
!                                                     **************
!                                                     ** Numbers. **
!                                                     **************
!
      hhs= 0.6
      call number (pl(i1)-0.8,ql(j1)-0.45,hhs,xmin,0.,101)
      call number (pr(i1)-1.1,ql(j1)-0.45,hhs,xmax,0.,101)
!
      call number (pl(i1)-1.8,ql(j1)     ,hhs,ymin,0.,101)
      call number (pl(i1)-1.8,qr(j1)-0.30,hhs,ymax,0.,101)
!
!                                                     **************
!                                                     **  Labels. **
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
!                                     **  No plot is made if npt1 < 0 **
!                                     **********************************
   70 if(npt1.lt.0) return
!
      isk= 0
      if(iskip.ne.0) then
        call set_width (2.0)
      end if
!
      if(iskip.eq.1) then
        itot= 5
        ibr = 3 ! 2
      end if
      if(iskip.eq.2) then
        itot= 13
        ibr = 9
      end if
      if(iskip.eq.7) then
        itot= 3
        ibr = 1
      end if
!
      call plotl (u(1),v(1),isc,3)
!**
      if(iplot.eq.1) then
         do 100 i=1,npt
         isk= isk +1
!
         if(iskip.eq.0) then
           call plotl (u(i),v(i),isc,2)
         else
           if(mod(isk,itot).lt.ibr) then
             call plotl (u(i),v(i),isc,2)
           else
             call plotl (u(i),v(i),isc,3)
           end if
         end if
  100    continue
      else
         do 120 i=1,npt-1
         call plotl (u(i+1),v(i)  ,isc,2)
         call plotl (u(i+1),v(i+1),isc,2)
  120    continue
      end if
!**
      call plotl (u(npt),v(npt),isc,3)
      call set_width (1.0)
!
      return
      end subroutine lplot1
!
!
!------------------------------------
      subroutine set_width (awid)
!------------------------------------
      use, intrinsic :: iso_c_binding
      real(C_float) amid
!
      write(77,*) 'stroke'
      write(77,10) awid
   10 format(f4.1,' setlinewidth')
!
      return
      end subroutine set_width 
!
!
!-----------------------------------------------------------------------
      subroutine cplot3 (q,xmax8,ymax8,zmax8,char,nc)
!-----------------------------------------------------------------------
!***********************************************************************
!*   contour plots of scalar quantities.                               *
!***********************************************************************
      use, intrinsic :: iso_c_binding
!
      include    'param-spinRL5.h'
      parameter  (mx1=meshx+1,my1=meshy+1,mz1=meshz+1)
      parameter  (mx=meshx,my=meshy,mz=meshz)
!
      character*8  char,label,cdate*10
      common/headr1/ label,cdate
      common/headr2/ time
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(mx1),pxc(mx1),pxl(mx1),pyr(my1),pyc(my1), &
                     pyl(my1),pzr(mz1),pzc(mz1),pzl(mz1)
      real(C_DOUBLE) q(mx,my,mz),xmax8,ymax8,zmax8
      real(C_float)  a(2048),b(2048),ww(2048),cut(200,4)
!
      j0= my/2 +1
      k0= mz/2 +1
      xmax= xmax8
      ymax= ymax8
      zmax= zmax8
!
!* 1. Plot at k= k0: subscript j first.
!
      npx= 0
      ij= 0
      qc = 1./16.
!***
      do 10 i= 1,mx
      npx= npx +1
      ir= pxr(i)
      il= pxl(i)
!
      npy= 0
      do 10 j= 1,my
      npy= npy +1
      jr= pyr(j)
      jl= pyl(j)
!
      ij= ij+1
      a(ij)= qc*(   q(ir,jr,k0) +2.*q(ir,j,k0)    +q(ir,jl,k0)  &
                +2.*q(i ,jr,k0) +4.*q(i ,j,k0) +2.*q(i ,jl,k0)  &
                +   q(il,jr,k0) +2.*q(il,j,k0)    +q(il,jl,k0) )
   10 continue
!
!* 2. Plot at j= j0: subscript k first.
!
      npx= 0
      ij= 0
      qc = 1./16.
!***
      do 20 i= 1,mx
      npx= npx +1
      ir= pxr(i)
      il= pxl(i)
!
      npz= 0
      do 20 k= 1,mz
      npz= npz +1
      kr= pzr(k)
      kl= pzl(k)
!
      ij= ij+1
      b(ij)= qc*(   q(ir,j0,kr) +2.*q(ir,j0,k)    +q(ir,j0,kl)  &
                +2.*q(i ,j0,kr) +4.*q(i ,j0,k) +2.*q(i ,j0,kl)  &
                +   q(il,j0,kr) +2.*q(il,j0,k)    +q(il,j0,kl) )
   20 continue
!
!
      hh = 0.70
      call symbol (0.1,18.2,hh,label,0.,8)
      call symbol (13.0,0.7,hh,cdate,0.,10)
      call symbol (13.0,0.1,hh,'t =',0.,3)
      call values (999.0,999.0,hh,time,0.,101)
!
      xl1=  1.8
      xr1=  9.3
      xl2= 10.0
      xr2= 17.5
!
      zl=  1.0
      zr=  zl +(xr1 -xl1)*zmax/xmax
      if(zr.gt.10.) zr= 10.
!                          <--- limit elongated y-length.
!
      yl=  zr +1.
      yr=  yl +(xr1 -xl1)*ymax/xmax
      if(yr.gt.25.) yr= 25.
!                          <--- limit elongated y-length.
!
      xc1= 0.5*(xr1+xl1)
      xc2= 0.5*(xr2+xl2)
      yc=  0.5*(yr+yl)
      zc=  0.5*(zr+zl)
!
!---------------------------------------------
!*  **Maximum of the vectors**
!---------------------------------------------
!
      am2= 0.
      am4= 0.
!
      do 100 ij= 1,npx*npy
      am2= amax1(am2,abs(a(ij)))
  100 continue
!
      do 200 ij= 1,npx*npz
      am4= amax1(am4,abs(b(ij)))
  200 continue
!
      ams= amax1(am2,am4)
      if(ams.lt.1.e-10) ams=999.0
!
      call symbol (zl,0.10,hh,'scalar.max=',0.,11)
      call values (999.0,999.0,hh,ams,0.,101)
!
!---------------------------------------------
!*  (1): Contours in (x,z) plane.
!---------------------------------------------
!
      call setscl (0.,0.,zmax,xmax,zl,xl2,zr,xr2,gdz,gdx, &
                   nc,char,6,' (x-z)',0.4, &
                   1,'z',0.4,1,'x',0.4,1)
!
      call values (zl-0.45,xl2-0.5,hh,0.0,0.,101)
      call values (zl-1.3,xr2-0.3, hh,xmax,0.,101)
      call values (zr-1.3,xl2-0.5, hh,zmax,0.,101)
!
      nxz= npx*npz
      call daisho (b,nxz,wamin,wamax)
!
      ncontr= 11
      call eqcntr (b,ww,npz,npx,zl,xl2,zr,xr2,wamin,0.0,wamax, &
                   ncontr,1)
!
!---------------------------------------------
!*  (2): Contours in (x,y) plane.
!---------------------------------------------
!
      call setscl (0.,0.,ymax,xmax,yl,xl2,yr,xr2,gdy,gdx, &
                   1,' ',1,' ',0.4, &  
                   1,'y',0.4,1,'x',0.4,1)
!
      call values (yl-0.45,xl2-0.5,hh,0.0,0.,101)
      call values (yl-1.3,xr2-0.3, hh,xmax,0.,101)
      call values (yr-0.3,xl2-0.5, hh,ymax,0.,101)
!
      nxy= npx*npy
      call daisho (a,nxy,wamin,wamax)
!
      ncontr= 11
      call eqcntr (a,ww,npy,npx,yl,xl2,yr,xr2,wamin,0.0,wamax, &
                   ncontr,1)
!
!---------------------------------------------
!*  (3): Cut plots.
!---------------------------------------------
!
      do 300 jj= 1,4
      j= (my/4)*(jj-1) +1
!
      do 300 i= 1,mx
      cut(i,jj)= q(i,j,k0)
  300 continue
!
!
      amax7= -1.e+10
      amin7=  1.e+10
!
      do 320 jj= 1,4
      do 320 i= 1,mx
      amax7= amax1(cut(i,jj),amax7)
      amin7= amin1(cut(i,jj),amin7)
  320 continue
!
      if(amax7.lt.0.) amax7= 0.
      if(amin7.gt.0.) amin7= 0.
!
!
      dd= amax7 -amin7
      dx= (yr -yl)/6.
      dy= xr1 -xl1
!
      do 340 jj= 1,4
      xo= yl +1.5*(jj-1)*dx
      xu= xo +dx
      call plot (xo,xl1,3)
      call plot (xu,xl1,2)
      call plot (xu,xr1,2)
      call plot (xo,xr1,2)
      call plot (xo,xl1,2)
!
!* zero line.
      x1= xo +dx*(0. -amin7)/dd
      call plot (x1,xl1,3)
      call plot (x1,xr1,2)
!
      x1= xo +dx*(cut(1,jj) -amin7)/dd
      y1= xl1 
      call plot (x1,y1,3)
!
      do 340 i= 1,mx
      x1= xo  +dx*(cut(i,jj) -amin7)/dd
      y1= xl1 +dy*(i-1)/float(mx-1)
!
      call plot (x1,y1,2)
  340 continue
!---------------------
      call chart
!---------------------
      return
      end subroutine cplot3 
!
!
!-----------------------------------------------------------------------
      subroutine eqcntr(u,w,nx,ny,xl,yl,xr,yr,umin,ubund,umax, &
                        lank,iwaku)
!-----------------------------------------------------------------------
!  << eqcntr >>
!         presented by kunihiko.watanabe 14.nov.1989
!         reviced   by hisanori.takamaru 16.mar.1990
!-----------------------------------------------------------------------
!     1. function
!        (1) to draw tokosen
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) u       nx,ny     (i)       world value
!        (2) w       nx,ny     (i)       work array (real*4)
!        (3) xl,xr,yl,yr       (i)       absolute coordinate value
!        (4) umin,umax         (i)       height of max & min
!                                        umin>umax : automatic control
!        (5) ubund             (i)       draw dash line (u < ubund)
!        (6) lank              (i)       number of draw lines
!        (7) iwaku             (i)       =1 : draw frame
!     3. called by
!             (** nothing **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      dimension u(1),w(1)
!
      if (nx.lt.2) return
      if (ny.lt.2) return
      if (xr.lt.xl) return
      if (yr.lt.yl) return
!
      nxy = nx*ny
      nxm1 = nx - 1
      nym1 = ny - 1
!
      dx = (xr-xl)/ float(nxm1)
      dy = (yr-yl)/ float(nym1)
!
      umax1 = umax
      umin1 = umin
!
      if(umax1.gt.(1.000001*umin1)) then
!
        do 10 i = 1 , nxy
          w(i) = u(i) - umin1
          if(u(i).gt.umax1) w(i) = umax1 - umin1
          if(u(i).lt.umin1) w(i) = 0.
   10   continue
!
      else
!
        umax1=-1.e+30
        umin1= 1.e+30
        do 20 i = 1 , nxy
          umax1=amax1(umax1,u(i))
          umin1=amin1(umin1,u(i))
   20   continue
        do 25 i = 1 , nxy
          w(i) = u(i) - umin1
   25   continue
!
      endif
!
!------------------------------------------------
      if(umax1.le.(1.000001*umin1))  return
!------------------------------------------------
!
      if(iwaku.eq.1) then
        call plot(xl,yl,3)
        call plot(xr,yl,2)
        call plot(xr,yr,2)
        call plot(xl,yr,2)
        call plot(xl,yl,2)
        call plot(xl,yl,3)
      endif
!
      uld = float(lank+1) / (umax1-umin1)
      eps = 1.0e-8
!
      nxym1 = nxm1*nym1
      do 9000  ijnxy1 = 1,nxym1
        j = (ijnxy1-1)/nxm1 + 1
        i = ijnxy1 - (j-1)*nxm1
!
          i1 = i + nx * (j - 1)
          i2 = i1 + 1
          i3 = i1 + 1 + nx
          i4 = i1 + nx
!
          u1 =  w(i1) * uld
          u2 =  w(i2) * uld
          u3 =  w(i3) * uld
          u4 =  w(i4) * uld
!
          k1 = ifix(u1)
          k2 = ifix(u2)
          k3 = ifix(u3)
          k4 = ifix(u4)
!
          j1 = iabs(k2-k1)
          j2 = iabs(k3-k2)
          j3 = iabs(k4-k3)
!
          if(j1.ne.0) then
            do 1000 ll = 1 , j1
              u0 = float(ll) + float(min0(k1,k2))
                ujouge = u0/uld + umin1
                if (ujouge.lt.ubund) then
                  jouge = 4
                else
                  jouge = 1
                end if
!
              if(abs(u2-u1).lt.eps)                 go to 1000
!
              x1 = xl + dx * ( (u0-u1)/(u2-u1) + float(i-1) )
              y1 = yl + dy * float(j-1)
!
              if( ((u3-u0)*(u2-u0)).gt.0. )         go to 1100
              if( ( (u0-u2).gt.0. ).and.( (u0-u4).gt.0. ) ) go to 1100
              if( abs(u3-u2).lt.eps )               go to 1100
!
                x2 = xl + dx * float(i)
                y2 = yl + dy * ( (u0-u2)/(u3-u2) + float(j-1) )
!
                call wdash(x1,y1,x2,y2,jouge)
!
 1100         continue
              if( ((u4-u0)*(u3-u0)).gt.0. )         go to 1200
              if( ((u1-u0)*(u3-u0)).gt.0. )         go to 1200
              if( ((u2-u0)*(u4-u0)).gt.0. )         go to 1200
              if( abs(u4-u3).lt.eps )               go to 1200
!
                x2 = xl + dx * ( (u0-u4)/(u3-u4) + float(i-1) )
                y2 = yl + dy * float(j)
                call wdash(x1,y1,x2,y2,jouge)
!
 1200         continue
              if( ((u1-u0)*(u4-u0)).gt.0. )         go to 1300
              if( ( (u0-u1).gt.0. ).and.( (u0-u3).gt.0. ) ) go to 1300
              if( abs(u1-u4).lt.eps )               go to 1300
!
                x2 = xl + dx * float(i-1)
                y2 = yl + dy*((u0-u1)/(u4-u1)+float(j-1))
                call wdash(x1,y1,x2,y2,jouge)
 1300         continue
 1000       continue
!
          endif
!
          if(j2.ne.0) then
!
            do 2000 ll = 1 , j2
              u0 = float(ll) + float(min0(k2,k3))
                ujouge = u0/uld + umin1
                if (ujouge.lt.ubund) then
                  jouge = 4
                else
                  jouge = 1
                end if
              if( abs(u3-u2).lt.eps )               go to 2000
!
              x1 = xl + dx * float(i)
              y1 = yl + dy * ( (u0-u2)/(u3-u2) + float(j-1) )
!
              if( ((u4-u0)*(u3-u0)).gt.0. )         go to 2100
              if( ( (u0-u1).gt.0. ).and.( (u0-u3).gt.0. ) ) go to 2100
              if( abs(u4-u3).lt.eps )               go to 2100
!
                x2 = xl + dx * ( (u0-u4)/(u3-u4) + float(i-1) )
                y2 = yl + dy * float(j)
!
                call wdash(x1,y1,x2,y2,jouge)
!
 2100         continue
              if( ((u1-u0)*(u4-u0)).gt.0. )         go to 2200
              if( ((u1-u0)*(u3-u0)).gt.0. )         go to 2200
              if( ((u2-u0)*(u4-u0)).gt.0. )         go to 2200
              if( abs(u1-u4).lt.eps )               go to 2200
!
                x2 = xl + dx * float(i-1)
                y2 = yl + dy * ( (u0-u1)/(u4-u1)+float(j-1) )
                call wdash(x1,y1,x2,y2,jouge)
 2200         continue
 2000       continue
!
          endif
!
          if(j3.ne.0) then
!
            do 3000 ll = 1 , j3
              u0 = float(ll) + float(min0(k3,k4))
                ujouge = u0/uld + umin1
                if (ujouge.lt.ubund) then
                  jouge = 4
                else
                  jouge = 1
                end if
              if( abs(u4-u3).lt.eps )               go to 3000
!
              x1 = xl + dx * ( (u0-u4)/(u3-u4) + float(i-1) )
              y1 = yl + dy * float(j)
!
              if( ((u1-u0)*(u4-u0)).gt.0. )         go to 3100
              if( ( (u0-u2).gt.0. ).and.( (u0-u4).gt.0. ) ) go to 3100
              if( abs(u1-u4).lt.eps )               go to 3100
!
                x2 = xl + dx * float(i-1)
                y2 = yl + dy * ( (u0-u1)/(u4-u1) + float(j-1) )
                call wdash(x1,y1,x2,y2,jouge)
 3100         continue
 3000       continue
          endif
 9000 continue
!
      return
      end subroutine eqcntr
!
!
!-----------------------------------------------------------------------
      subroutine setscl (wminx,wminy,wmaxx,wmaxy, xl,yl,xr,yr,gdx,gdy, &
                 n1,char1,n2,char2, hight1,                            &
                 nnx,charx,hightx, nny,chary,highty, iwaku)
!-----------------------------------------------------------------------
!  << setscl >>                   /char1/
!                          wmaxy  +--------------------+  (xl,yl)
!                               y |             (xr,yr)|  (xr,yr) on 0
!                               r |                    |
!                               a |                    |
!    (wminx,wminy)              h |                    |
!    (wmaxx,wmaxy) on is        c |                    |
!                                 |(xl,yl)             |
!                          wminy  +--------+--+--------+
!                                 wminx  /charx/       wmaxx
!-----------------------------------------------------------------------
!
!     setscl
!
!     1. function
!        (1) to scale the graphics by calcomp specifications
!     2. arguments            (i/o)     (meaning)
!        (1) wminx,wmaxx,
!            wminy,wmaxy       (i)       world coordinate value
!        (2) xl,xr,yl,yr       (i)       absolute coordinate value
!        (3) gdx,gdy           (o)       scaling factor of coordinate
!                                        from world to absolute
!        (4) char1,charx,cahry (i)       title on graph,x-axis,y-axis
!        (5) iwaku             (i)       draw frame (0:off ; 1:on)
!                                         999 : write out only title,
!                                                    not draw otherwise
!     3. called by
!             (** nothing **)
!     4. calls
!             (** plot   **)
!             (** symbol **)
!             (** number **)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      character*1  char1(1),char2(1),charx(1),chary(1)
!
      if (wmaxx.le.wminy) goto 9999
      if (wmaxx.le.wminy) goto 9999
      if (xr.le.xl)       goto 9999
      if (yr.le.yl)       goto 9999
!
      gdx= (xr-xl)/(wmaxx-wminx)
      gdy= (yr-yl)/(wmaxy-wminy)
!
      xc = 0.5*( xr + xl )
      yc = 0.5*( yr + yl )
!
      if (n1 .gt.0) then
        if (hight1.gt.0) then
          xs1= xc -0.5*n1*hight1
          xs2= xs1 +(n1+1)*hight1
          call symbol(xs1,yr+0.1,hight1,char1(1),0.,n1)
          call symbol(xs2,yr+0.1,hight1,char2(1),0.,n2)
        end if
      end if
!-----------------------------------------------------------------------
      if (iwaku.eq.999) return
!-----------------------------------------------------------------------
!
      if (iwaku.eq.1) then
        call plot (xl,yl,3)
        call plot (xl,yr,2)
        call plot (xr,yr,2)
        call plot (xr,yl,2)
        call plot (xl,yl,2)
        call plot (999.,999.0,3)
      end if
!
      if (nnx.gt.0) then
        if (hightx.gt.0) then
          call symbol(xc-0.5*hightx*nnx,yl-0.5,hightx,charx(1),0.,1)
          do 200 nnx1=2,nnx
  200     call symbol(999.0,999.0,hightx,charx(nnx1),0.,1)
        end if
      end if
      if (nny.gt.0) then
        if (highty.gt.0) then
          call symbol(xl-0.5,yc-0.5*highty*nny,highty,chary(1),0.,1)
          do 300 nny1=2,nny
  300     call symbol(999.0,999.0,highty,chary(nny1),0.,1)
        end if
      else if(nny.lt.0) then
        if (highty.gt.0) then
          call symbol(xc-0.5*highty*nny,yc,highty,chary(1),0.,1)
          do 400 nny1=2,nny
  400     call symbol(999.0,999.0,highty,chary(nny1),0.,1)
        end if
      end if
!
      return
!
!-----------------------------------------------------------------------
!
 9999 continue
      write(6,*) '**********  abnormal world coordinate ********'
      write(6,*) '      '
      write(6,*) '    wmaxx =',wmaxx,' wminx = ',wminx
      write(6,*) '    wmaxy =',wmaxy,' wminy = ',wminy
      write(6,*) '    xl,yl,xr,yr =',xl,yl,xr,yr
      write(6,*) '    fctr  =',fctr
      write(6,*) '      '
      call chart
      call symbol(1.0,10.0,0.2,' abnormal world coordinate call',0.,31)
      call symbol(1.0,09.0,0.2,' wmaxx =',0.,8)
      call number(999.0,999.0,0.2,wmaxx,0.,2)
      call symbol(1.0,08.5,0.2,' wminx =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,08.0,0.2,' wmaxy =',0.,8)
      call number(999.0,999.0,0.2,wmaxy,0.,2)
      call symbol(1.0,07.5,0.2,' wminy =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,07.0,0.2,' fctr  =',0.,8)
      call number(999.0,999.0,0.2,fctr,0.,2)
      call symbol(1.0,06.5,0.2,' xleft =',0.,8)
      call number(999.0,999.0,0.2,xl,0.,2)
      call symbol(1.0,06.0,0.2,' yleft =',0.,8)
      call number(999.0,999.0,0.2,yl,0.,2)
      call symbol(1.0,05.5,0.2,' xright=',0.,8)
      call number(999.0,999.0,0.2,xr,0.,2)
      call symbol(1.0,05.0,0.2,' yright=',0.,8)
      call number(999.0,999.0,0.2,yr,0.,2)
      call exit
      return
      end subroutine setscl 
!
!
!-----------------------------------------------------------------------
      subroutine scalex (xcm,ycm,x00,y00,dx,dy,isc)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
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
      end subroutine scalex 
!
!
!-----------------------------------------------------------------------
      subroutine plotl (x,y,isc,ipl)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      xcm= xl(isc) +dxi(isc)*(x -x0(isc))
      ycm= yl(isc) +dyi(isc)*(y -y0(isc))
!
      call plot (xcm,ycm,ipl)
!
      return
      end subroutine plotl 
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
      use, intrinsic :: iso_c_binding
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
        call number(x,y,height,val,theta,-1)
      else if (n100.eq.1) then
        call number(x,y,height,val,theta,keta)
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
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or. &
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
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or. &
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
      end subroutine values 
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
      use, intrinsic :: iso_c_binding
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
      end subroutine wdash 
!
!
!-----------------------------------------------------------------------
      subroutine daisho(x,nx,xmin1,xmax1)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      dimension x(1)
!
      xmax1= x(1)
      xmin1= x(1)
      do 100 i=2,nx
      xmax1= amax1(xmax1,x(i) )
      xmin1= amin1(xmin1,x(i) )
  100 continue
!
      return
      end subroutine daisho
!
!
!***************************************************************
!*   This program package generates a unix postscript          *
!*   graphic file when called by calcomp-compatible /plot23.f/ *.  
!***************************************************************
!----------------------------------------------------------
!    Postscript header by fortran
!        t. ogino (nagoya university) February 27, 1992
!      modified to conform gsipp commands
!        Motohiko Tanaka (nifs)       November 23, 1993
!
!----------------------------------------------- 5/27/96 -------
!   This ps-adobe-2.0 header allows us full paging features in
!   the ghostview.  to scroll up the page (backward), click the 
!   page number and press two buttons of mouse simultaneously.
!
!   consult: A.Saitou (kyoto u.)  the definition of /@eop  
!   needs stroke for line drawings (not in the tex header).
!---------------------------------------------------------------
       subroutine gopen (nframe)
!----------------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
!*  this is an adobe-2.0 postscript file.
!
       write(77,10)
   10  format('%!ps-adobe-2.0',/      &
              '%%pages: (atend)',/    &
              '%%pageorder: ascend',/ &
              '%%endcomments',/       &
              '%%begindocument')
!
!%%%%%%%%%%%%%%%%%%% procedure defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     write(77,11) 
!  11 format('%%boundingbox: 150. 400. 550. 600.')
!
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/ &
             '/m {moveto} bind def  % x y m -- move to position')
!
      write(77,23) 
   23 format('/tr {/times-roman findfont} bind def',/ &
             '/sf {scalefont} bind def',/  &
             '/se {setfont} bind def',/    &
             '/ro {rotate}  bind def',/    &
             '/tl {translate} bind def',/  &
             '/sc {scale} bind def')
!
      write(77,24) 
   24 format('/@bop          % @bop -- begin the a new page',/ &
             '{erasepage newpath initgraphics',/ &
             '/saveimage save def',/ &
             '} bind def')
!
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/ &
             '{stroke showpage',/   &
             ' saveimage restore',/ &
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
             ' show grestore}',/ &
             ' {show} ifelse',/  &
             '} bind def')
!
      write(77,31)
   31 format('%%enddocument',/ &
             '%%endprolog',/   &
             '%%beginsetup',/  &
             '/resolution 300 def',/ &
             '/#copies 1 def',/ &
             '%%endsetup')
!
!%%%%%%%%%%%%%%%%%%% end of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*  initiate the page one.
!
       nfrm = nframe
!
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%page:',1x,i2,1x,i2)
!
       write(77,30) 
   30  format('%%beginpagesetup',/ &
              '%%endpagesetup',/   &
              '@bop')
!
!
!*  Set magnifying factor (gsipp to sun coordinate).
!   rotate and translate to output on a4-l paper.
!      left corner ...... (  0.,  0.)
!      right corner ..... (600.,780.)
!
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
!
       write(77,*) '90.0 ro'
       write(77,*) '50.0 -550.0 tl'
!
!*  if nfrm=4, four frames in a page (top-left frame).
!
       if(nfrm.eq.1) then
          write(77,*) '1.00 1.00 sc'
       else
          write(77,*) '0.50 0.50 sc'
          write(77,*) '0.0 550.0 tl'
       end if
!
       return
       end subroutine gopen 
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
!*     four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
!
!
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
!
!*  frame 1: open a new page.
!
       if(loc.eq.0) then
          call plote
!
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
!
          write(77,10) 
   10     format('%%pagetrailer    % need for the page count')
!
          write(77,20) lpage,lpage
   20     format('%%page:',1x,i2,1x,i2)
!
          write(77,30) 
   30     format('%%beginpagesetup',/ &
                 '%%endpagesetup',/   &
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
!-----------------------------------------------------
!   First cancel the previous translation, then
!   make a new translation (scale factor alive).
!-----------------------------------------------------
!*   frames 2-4:
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
       end subroutine chart
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
       end subroutine newcolor 
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
!
       return
       end subroutine newpen 
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
!
       return
       end subroutine plot 
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       character isymb*80,ica*80,ich(80)*1
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
       end subroutine symbol 
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
       write(isymb,40) anu
   40  format(1pe9.2)
       write(77,*) '(',isymb,') s'
!
       return 
       end subroutine number 
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
       end subroutine sunscl 
