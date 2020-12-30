      program res_arg_nu1c_mult

      double precision time,a,e,sini,comega,omega,M, varpi
      double precision ac,ec,sinic,comegac,omegac,MC, varpic
      double precision res_ang_sin,res_ang_cos,res_ang, sigma
      double precision degtorad, pi
      character*16 file_el,file_res
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      npart=50
      nint=1123
      
      degtorad=1.74532925199e-2
      pi=3.14159265358979
      
      open(1,file='mars')

      do l=1,npart
         write(file_el,3000) l
         write(file_res,4000) l
         write(*,*) file_el,file_res
         rewind(1)
         close(2)
         open(2,file=file_el)
         close(3)
         open(3,file=file_res)
         do i=1,nint
            read(1,*) time,ac,ec,sinic,comegac,omegac,MC
            read(2,*) time,a,e,sini,comega,omega,M
            varpi=comega+omega
            varpic=comegac+omegac
c            sigma=2*(M+varpi)-(MC+varpic)-varpi
            sigma=2*(M+varpi)-(MC+varpic)-varpic
            res_ang_sin=sin((sigma)*pi/180.)
            res_ang_cos=cos((sigma)*pi/180.)
            res_ang=atan2(res_ang_sin,res_ang_cos)*180./pi
            if(res_ang.le.0) res_ang=res_ang+360.
            write(3,777) time/1d6,res_ang_cos,res_ang_sin,res_ang
         enddo
      enddo
 777  format(4(f11.4,2x))
 3000 format('el_osc_s25_',i2.2)
 4000 format('res_arg_',i2.2)
 999  end
