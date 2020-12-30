c converts binary file to ascii file

        include 'swift.inc'

        real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
        real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

        real*8 mass(NPLMAX),dr,peri
        real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
        real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

        integer istat(NTPMAX,NSTAT)
        real*8 rstat(NTPMAX,NSTATR)
        integer nbod,ntp,ierr,ifol,istep
        integer iflgchk,iu,nleft,i,id
        integer io_read_hdr,io_read_line
        integer io_read_hdr_r,io_read_line_r,n_part,l

        real*8 t0,tstop,dt,dtout,dtdump
        real*8 t,tmax

        real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
        logical*2 lclose
        real*8 a,e,inc,capom,omega,capm,j2rp2,j4rp4
        real*8 elh,elk,elp,elq,apo

        character*80 outfile,inparfile,inplfile,intpfile,fopenstat
	character*20 file1

        real*8 a0,e0,inc0,da,de,dinc

c	degtorad=1.74532925199e-2

c Get data for the run and the test particles
        write(*,*) 'Enter name of parameter data file : '
        read(*,999) inparfile
c       call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
c     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
        write(*,*) ' '
        write(*,*) 'Enter name of planet data file : '
        read(*,999) inplfile
999     format(a)
c       call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
c     &       vxh,vyh,vzh,rplsq,j2rp2,j4rp4)



c Get data for the run and the test particles
        write(*,*) 'Enter name of test particle data file : '
        read(*,999) intpfile
c       call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
c     &               vzht,istat,rstat)

c       inparfile = "param_flora.in"
        call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c       inplfile = "pl.in"
        call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &       vxh,vyh,vzh,rplsq,j2rp2,j4rp4)


c       intpfile = "tp.in"
        call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &               vzht,istat,rstat)


        write(*,*) 'Enter name of output data file : '
        read(*,999) outfile

        n_part=ntp

        do l=1,n_part
           if(l.le.9) then
              write(file1,3000) l
c              write(99,222) file1
           elseif(l.gt.9) then
              if(l.le.99)then
                 write(file1,3001) l
c                 write(99,222) file1
              elseif(l.gt.99) then
                 write(file1,3002) l
c                 write(99,222) file1
              endif
           endif
 222       format(a13)
 3000      format('el_osc_s25_0',i1)
 3001      format('el_osc_s25_',i2)
 3002      format('el_osc_s25_',i3)


           if(l.le.4) then
              call io_open(l,file1,'unknown','formatted',ierr)
c	write(l,*) 1
           elseif(l.gt.4.5) then
              call io_open(l+4,file1,'unknown','formatted',ierr)	   
c	   write(l+3,*) 1
              call io_open(l+4,file1,'unknown','formatted',ierr)	   
           endif	 
c	 write(99,222) file1
c           write(*,*) file1
        enddo

        iu = 99

        dr = 180.0/PI

        if(btest(iflgchk,0)) then
           write(*,*) ' Reading an integer*2 binary file '
        else if(btest(iflgchk,1)) then
           write(*,*) ' Reading an real*4 binary file '
        else
           write(*,*) ' ERROR: no binary file format specified '
           write(*,*) '        in param file '
           stop
        endif

c        write(*,*) ' Input the particle number to follow '
c        read(*,*) ifol
c        write(*,*) ' Following particle ',ifol

        open(unit=iu, file=outfile, status='old',form='unformatted')
        open(unit=7,file='follow_all.dat')
        open(unit=5,file='plan.dat')

        write(*,*) '1 2 3  4    5     6    7    8    9 '
        write(*,*) 't,a,e,inc,capom,omega,capm,peri,apo'

        tmax = t0
 1      continue
             if(btest(iflgchk,0))  then ! bit 0 is set
c                ierr = io_read_hdr(iu,t,nbod,nleft)
             else
                ierr = io_read_hdr_r(iu,t,nbod,nleft)
             endif

             if(ierr.ne.0) goto 2

             istep = 0
             do i=2,nbod
                if(btest(iflgchk,0))  then ! bit 0 is set
c                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm)
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm)
                endif
                if(ierr.ne.0) goto 2
c                if(id.eq.ifol) then
                   istep = 1
                   elh = e*cos(omega+capom)
                   elk = e*sin(omega+capom)
                   elp = sin(inc/2.0)*cos(capom)
                   elq = sin(inc/2.0)*sin(capom)
                   inc = inc*dr
                   capom = capom*dr
                   omega = omega*dr
                   capm = capm*dr
                   peri = a*(1.0d0-e)
                   apo = a*(1.0d0+e)
                   write(5,1010) id,t/365.25,a,e,inc,capom,omega,
     &                  capm,peri,apo
c 1000              format(1x,e15.7,1x,f10.4,1x,f7.5,4(1x,f9.4),
c     &                  2(1x,f10.4))
c                   write(5,1001) id,t,elh,elk,elp,elq
 1010              format(i2,1x,f15.4,1x,1p,e15.7,1x,0p,f10.4,
     &                  1x,f7.5,4(1x,f9.4),2(1x,f10.4))         
 1001              format(1x,i3,1x,e13.5,4(1x,f7.5))
                   tmax = t
c                endif
             enddo

c            close(5)

             do i=1,nleft
                if(btest(iflgchk,0))  then ! bit 0 is set
c                   ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm)
                else
                   ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm)
                endif
                if(ierr.ne.0) goto 2
c                if(id.eq.ifol) then

                   if(t.eq.0.) then
                      a0 = a
                      e0 = e
                      inc0 = inc
                   end if

                   istep = 1

                   elh = e*cos(omega)
                   elk = e*sin(omega)
                   elp = sin(inc/2.0)*cos(capom)
                   elq = sin(inc/2.0)*sin(capom)
                   inc = inc*dr
                   capom = capom*dr
                   omega = omega*dr
                   capm = capm*dr
                   peri = a*(1.0d0-e)
                   apo = a*(1.0d0+e)
                write(7,1000) id,t,a,e,inc,capm,omega,capom
		if(id.le.4) then
		   write(id,1002) t/365.25,a,e,inc,capom,omega,capm
		elseif(id.gt.4.5) then
		   write(id+4,1002)  t/365.25,a,e,inc,capom,omega,capm
		endif


 1000           format(i4,1x,1p,e15.7,1x,0p,f10.4,1x,f7.5,4(1x,f9.4),
     &                  2(1x,f10.4))
 1002           format(e15.7,1x,0p,f10.5,1x,f10.5,4(1x,f10.5),
     &                  2(1x,f10.5))

                   tmax = t
c                   write(8,1001) t,elh,elk,elp,elq

c                endif
             enddo
             if(istep.eq.0) goto 2     ! did not find particle this times step

        goto 1

 2      continue

        write(*,*) ' Tmax = ',tmax

        stop
        end

c*************************************************************************
c                            IO_OPEN.F
c*************************************************************************
c                  THIS FILE MUST BE PRECOMPILED
c*************************************************************************
c open files
c
c             Input:
c                 iu              ==>  unit number (integer scalar)
c                 fname           ==>  file name (character*80)
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c                 format          ==>  format string (character*80)
c             Output:
c                 ierr            ==>  output from iostat
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    3/3/94
c Last revision: 1/30/98

      subroutine io_open(iu,fname,fopenstat,format,ierr)

      include 'swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu
      character*(*) fname,fopenstat,format

c...  Outputs: 
      integer ierr

c----
c...  Executable code 

      if( (fopenstat(1:6).eq.'append') .or. 
     &     (fopenstat(1:6).eq.'APPEND') ) then
         open(unit=iu, file=fname, status='old',
c#ifdef  _OPEN_POSITION
c     &        position='append',
c#else
c     &        access='append',
c#endif
     &        form=format,iostat=ierr)
         if(ierr.ne.0) then
            write(*,*) 'Warning:  Could not open ',fname,' with'
c#ifdef  _OPEN_POSITION
            write(*,*) '          position=append.'
c#else
            write(*,*) '          access=append.'
c#endif
            write(*,*) '          Will open as status=new'
            open(unit=iu, file=fname, status='new',
     &           form=format,iostat=ierr)
         endif
      else
         open(unit=iu, file=fname, status=fopenstat,
     &        form=format,iostat=ierr)

      endif


      return
      end   ! io_open
c************************************************************************
c                              IO_INIT_PL.F
c************************************************************************
c IO_INIT_PL reads in the data for the Sun and planets 
c
c             Input:
c                 infile        ==> File name to read from (character*80)
c                 lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c                 iflgchk        ==>  bit 5 set ==>  include J2 and J4 terms
c
c             Output:
c                 nbod          ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  initial position in Helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  initial position in Helio coord 
c                                    (real arrays)
c                 rplsq         ==>  min distance^2 that a tp can get from pl
c                                    (real array)
c                 j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c                                     (real scalars)
c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision: 9/30/00  HFL

	subroutine io_init_pl(infile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &     vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

	include 'swift.inc'
	include 'io.inc'

c...    Input
	character*(*) infile
	integer iflgchk
        logical*2 lclose

c...    Output
	real*8 mass(NPLMAX),rplsq(NPLMAX),j2rp2,j4rp4
	real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
	real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)
	integer nbod

c...   Internal
	integer j,ierr
        real*8 rpl

c-----
c...  Executable code      

	write(*,*) 'Planet data file is ',infile
        call io_open(7,infile,'old','formatted',ierr)

c Read number of planets
	read(7,*) nbod

        if(nbod.gt.NPLMAX) then
           write(*,*) ' SWIFT ERROR: in io_init_pl: '
           write(*,*) '   The number of massive bodies,',nbod,','
           write(*,*) '   is too large, it must be less than',NPLMAX
           call util_exit(1)
        endif

	write(*,23) nbod
23	format(/,'Number of bodies (incl. the Sun) is ',i3,/,
     &   'For each, list mass ',/,
     &   'Followed by x,y,z,vx,vy,vz : '/)

c For each planet read mass, 
c and helioc. position and vel .
        if(btest(iflgchk,5))  then ! bit 5 is set
           read(7,*) mass(1),j2rp2,j4rp4
        else
           read(7,*) mass(1)
           j2rp2 = 0.0d0
           j4rp4 = 0.0d0
        endif
        read(7,*) xh(1),yh(1),zh(1)
        read(7,*) vxh(1),vyh(1),vzh(1)
        write(*,*) mass(1)
        write(*,*) xh(1),yh(1),zh(1)
        write(*,*) vxh(1),vyh(1),vzh(1)
        rplsq(1) = 0.0d0

        if(  (xh(1).ne.0.0d0) .or.
     &       (yh(1).ne.0.0d0) .or.
     &       (zh(1).ne.0.0d0) .or.
     &       (vxh(1).ne.0.0d0) .or.
     &       (vyh(1).ne.0.0d0) .or.
     &       (vzh(1).ne.0.0d0) ) then
           write(*,*) ' SWIFT ERROR: in io_init_pl: '
           write(*,*) '   Input MUST be in heliocentric coordinates '
           write(*,*) '   Position and Vel. of Massive body 1 .ne. 0'
           call util_exit(1)
        endif

	do j=2,nbod
           if(lclose) then
              read(7,*) mass(j),rpl
              write(*,*) mass(j),rpl
              rplsq(j) = rpl*rpl
           else
              read(7,*) mass(j)
              write(*,*) mass(j)
           endif
	   read(7,*) xh(j),yh(j),zh(j)
	   read(7,*) vxh(j),vyh(j),vzh(j)
	   write(*,*) xh(j),yh(j),zh(j)
	   write(*,*) vxh(j),vyh(j),vzh(j)
	enddo

	close(unit = 7)
	return
	end     ! io_init_pl.f
c--------------------------------------------------------------------------

c**********************************************************************
c			IO_INIT_TP.F
c**********************************************************************
c Read in test particle data
c
c             Input:
c                 infile        ==> File name to read from (character*80)
c
c             Output:
c                 ntp           ==>  number of massive bodies (int scalar)
c                 mass          ==>  mass of bodies (real array)
c              xht,yht,zht      ==>  initial position in Helio coord 
c                                    (real arrays)
c              vxht,vyht,vzht   ==>  initial position in Helio coord 
c                                    (real arrays)
c               istat           ==>  status of the test paricles
c                                      (2d  integer array)
c                                      istat(i,1) = 0  active
c                                      istat(i,1) = 1 not
c               rstat           ==>  status of the test paricles
c                                      (2d  real array)
c
c
c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision:  12/22/95  HFL

	subroutine io_init_tp(infile,ntp,xht,yht,zht,vxht,vyht,
     &     vzht,istat,rstat)

	include 'swift.inc'
	include 'io.inc'

c...    Input
	character*(*) infile

c...    Output
	real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
	real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)
        real*8 rstat(NTPMAX,NSTATR)
	integer istat(NTPMAX,NSTAT)
	integer ntp

c...   Internal
	integer i,j,ierr,ns

c-----
c...  Executable code      

	write(*,*) 'Test particle file called ',infile
        call io_open(7,infile,'old','formatted',ierr)

	read(7,*) ntp

        if(ntp.gt.NTPMAX) then
           write(*,*) ' SWIFT ERROR: in io_init_tp: '
           write(*,*) '   The number of test bodies,',ntp,','
           write(*,*) '   is too large, it must be less than',NTPMAX
           call util_exit(1)
        endif

	write(*,*) ' '
	write(*,*) 'ntp : ',ntp

        if(ntp.eq.0) then
           close(unit = 7)
           write(*,*) ' '
           return
        endif               ! <===== NOTE

c...   Determine the number of istat and rstat variables.  In what follows,
c...   we assume that they are the same.

        call io_getns(7,ns)
        
        if(ns.ne.NSTAT) then
           write(*,*) 'Warning:  The size of istat and rstat arrays is '
           write(*,*) '          not NSTAT=',NSTAT,', but is ',ns
        endif

c Start again:
        rewind(7)
        read(7,*) ntp 

c Read in the x's and v's and istat(*,*)
	  write(*,*) ' '
	  do  i=1,ntp
	    read(7,*) xht(i),yht(i),zht(i)
	    read(7,*) vxht(i),vyht(i),vzht(i)
	    read(7,*) (istat(i,j),j=1,ns)
	    read(7,*) (rstat(i,j),j=1,ns)
            do j=ns+1,NSTAT
               istat(i,j) = 0
            enddo
            do j=ns+1,NSTATR
               rstat(i,j) = 0.0d0
            enddo
	  enddo

	close(unit = 7)
        write(*,*) ' '

	return
	end    ! io_init_tp.f
c-----------------------------------------------------------------

c*************************************************************************
c                            IO_READ_HDR_R
c*************************************************************************
c read in header part of the real*4 file
c
c             Input:
c                 iu            ==> unit number to write to
c             Output:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 nleft         ==>  number of active tp (int scalar)
c
c             Returns:
c               io_read_hdr_r     ==>   =0 read ok
c                                    !=0 read failed is set to iostat variable
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

      integer function io_read_hdr_r(iu,time,nbod,nleft) 

      include 'swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Output
      integer nbod,nleft
      real*8 time

c...  Internals
      real*4 ttmp
      integer*2 nleft2,nbod2
      integer ierr

c----
c...  Executable code 


      read(iu,iostat=ierr) ttmp,nbod2,nleft2
      io_read_hdr_r = ierr
      if(ierr.ne.0) then
         return
      endif

      nbod = nbod2
      nleft = nleft2
      time = ttmp

      return
      end     ! io_read_hdr_r.f
c---------------------------------------------------------------------------

c*************************************************************************
c                            IO_READ_LINE_R
c*************************************************************************
c read one line from real*4 binary file.
c
c      Input:
c            iu       ==> unit number to write to
c      Output:
C	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C	     omega    ==> argument of perihelion (real scalar)
C	     capm     ==> mean anomoly(real scalar)
c       Returns:
c      io_read_line_r    ==>   =0 read ok
c                           !=0 read failed is set to iostat variable
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

      integer function io_read_line_r(iu,id,a,e,inc,capom,omega,capm) 

      include 'swift.inc'
      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Output: 
      integer id
      real*8 a,e,inc,capom,omega,capm

c...  Internals
      integer*2 id2
      real*4 a4,e4,inc4,capom4,omega4,capm4
      integer ierr

c----
c...  Executable code 

      read(iu,iostat=ierr) id2,a4,e4,inc4,capom4,omega4,capm4
      io_read_line_r = ierr
      if(ierr.ne.0) then
         return
      endif

      id = id2

      a = a4
      e = e4
      inc = inc4
      capom = capom4
      capm = capm4
      omega = omega4

      return
      end      ! io_read_line_r
c--------------------------------------------------------------------------

c*************************************************************************
c                            UTIL_EXIT.F
c*************************************************************************
c Exits program
c
c             Input:
c                 iflg          ==>  status of exit
c                                       = 0 if normal exit
c                                       = 1 if exit because error
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    8/6/93
c Last revision: MD : change calc. of rhil Apr. 25

      subroutine util_exit(iflg)

      include 'swift.inc'

c...  Inputs: 
      integer iflg


c-----
c...  Executable code 

      write(*,*) ' '

      if(iflg.eq.0) then
        write(*,1000) VER_NUM 
 1000   format('Normal termination of SWIFT (version ',f3.1,')')
      else
        write(*,2000) VER_NUM 
 2000   format('Terminating SWIFT (version',f3.1,') due to ERROR!!! ')
      endif

      write(*,*) '----------------------------------------------------'

      stop
      end  ! util_exit

c---------------------------------------------------


c************************************************************************
c                          IO_INIT_PARAM.F
c************************************************************************
c INIT_PARAM reads in the parameters for the integration. 
c
c      Input:
c            infile   ==> File name to read from (character*80)
c
c      Output:
c            t0       ==> Initial time (real scalar)
c            tstop    ==> final time (real scalar)
c            dt       ==> time step  (real scalar)
c            dtout    ==> time between binary outputs (real scalar)
c            dtdump   ==> time between dumps  (real scalar)
c            iflgchk  ==>  =0 don't run diagnostic routines
c                          bit 0 set ==>  write int*2 binary data file
c                          bit 1 set ==>  write real*4 binary file 
c                          bit 2 set ==>  calc energy of system wrt time
c                          bit 3 set ==>  calc jacobi of the test particles
c                          bit 4 set ==>  check if particles are removed
c                          bit 5 set ==>  include J2 and J4 terms
c      rmin,rmax      ==>  maximum and min distance from Sun
c                                if <0  then don't check
c                                    (real scalar)
c      rmaxu          ==>  maximum distance from Sun in not bound
c                                 if <0  then don't check
c                                      (real scalar)
c       qmin          ==> Smallest perihelion distance
c                                 if <0  then don't check
c                                      (real scalar)
c       lclose        ==> .true. --> discard particle if it gets 
c                                    too close to a planet. Read in that 
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c       outfile       ==>  Name of binary output file (character*80)
c       fopenstat     ==>  The status flag for the open statements of the
c                          output files.  Must be one of the following:
c                                 new      (die if the file exists)
c                                 append   (add to what is there)
c                                 unknown  (just write over what is there)
c                                 (character*80)
c
c
c Remarks: 
c Authors:  Martin Duncan
c Date:    3/2/93 
c Last revision:  5/10/94  HFL

        subroutine io_init_param(infile,t0,tstop,dt,dtout,dtdump,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

	include 'swift.inc'
	include 'io.inc'

c...    Input
	character*(*) infile

c...  Outputs: 
	integer iflgchk
	real*8 t0,tstop,dt
	real*8 dtout,dtdump
	real*8 rmin,rmax,rmaxu,qmin
        logical*2 lclose
	character*80 outfile,fopenstat

c...  Internals
        logical*1 lflg(0:IO_NBITS-1)
        integer i,ierr

c-----
c...  Executable code 

	write(*,*) 'Parameter data file is ',infile
        call io_open(7,infile,'old','formatted',ierr)

	read(7,*) t0,tstop,dt
	write(*,*) 't0,tstop,dt : ',t0,tstop,dt
	read(7,*) dtout,dtdump
	write(*,*) 'dtout,dtdump : ',dtout,dtdump
        read(7,*) (lflg(i),i=IO_NBITS-1,0,-1)

        iflgchk=0
        do i=0,IO_NBITS-1
           if(lflg(i)) then
              iflgchk = ibset(iflgchk,i)
           endif
        enddo

        write(*,*) (lflg(i),i=IO_NBITS-1,0,-1),' = ',iflgchk

        if(btest(iflgchk,0) .and. btest(iflgchk,1))  then 
           write(*,*) ' SWIFT ERROR: in io_init_param:'
           write(*,*) '    Invalid logical flags '
           write(*,*) '    You cannot request that both a real and ',
     &                '       an integer binary file be written '
           call util_exit(1)
        endif

        if(btest(iflgchk,4))  then ! bit 4 is set
           read(7,*) rmin,rmax,rmaxu,qmin,lclose
           write(*,*) 'rmin,rmax,rmaxu,qmin,lclose :',
     &          rmin,rmax,rmaxu,qmin,lclose
        else
           rmin = -1.0
           rmax = -1.0
           rmaxu = -1.0
           qmin = -1.0
           lclose = .false.
        endif

        if(btest(iflgchk,0) .or. btest(iflgchk,1))  then 
           read(7,999) outfile
 999       format(a)
           write(*,*) 'outfile : ', outfile
           write(*,*) ' '
        endif

        read(7,999) fopenstat
        if(  (fopenstat(1:3).ne.'new') .and. 
     &       (fopenstat(1:3).ne.'NEW') .and.
     &       (fopenstat(1:7).ne.'unknown') .and. 
     &       (fopenstat(1:7).ne.'UNKNOWN') .and.
     &       (fopenstat(1:6).ne.'append') .and. 
     &       (fopenstat(1:6).ne.'APPEND') ) then
           write(*,*) ' SWIFT ERROR: in io_init_param:'
           write(*,*) '    Invalid status flag:',fopenstat(1:7),':'
           call util_exit(1)
        endif
        
	close(unit = 7)

	return
	end     ! io_init_param
c____________________________________________________________________________
c
c

c**********************************************************************
c			IO_GETNS.F
c**********************************************************************
c Determines that number of istat variables there are in a tp.in file
c
c             Input:
c                 iu            ==> unit number (integer)
c
c             Output:
c                 ns            ==>  number of istat variable (int scalar)
c
c
c Remarks: 
c Authors:  Hal Levison
c Date:    10/1/96
c Last revision:  3/18/97

      subroutine io_getns(iu,ns)

      include 'swift.inc'
      include 'io.inc'

c...  Input
      integer iu

c...  Output
      integer ns
      
c...  Internal
      character*1024 line
      integer i,i1,ib
      real*8 xht,yht,zht    
      real*8 vxht,vyht,vzht

c-----
c...  Executable code      

c...  get the irrelavant stuff
      read(7,*) xht,yht,zht    
      read(7,*) vxht,vyht,vzht

      ns = 0
      do while(.true.)
         read(7,fmt='(a)') line

c...     if there are `.' then it is not a istat line. Therefore leave
         do i = 1,1024
            if(line(i:i).eq.'.') goto 99
         enddo


c...     Find the first non-blank character
         i1 = 0
         do i=1,1024
            if( (i1.eq.0) .and. (line(i:i).ne.' ') ) then
               i1 = i
            endif
         enddo

         ib = 1
         do i=i1+1,1024
            if( (ib.eq.1) .and. (line(i:i).eq.' ') ) then
               ns = ns + 1
            endif
            if (line(i:i).eq.' ') then
               ib = 0
            else
               ib = 1
            endif
         enddo

      enddo
 99   continue

      rewind(7)

      return
      end     ! io_getns
c----------------------------------------------------


