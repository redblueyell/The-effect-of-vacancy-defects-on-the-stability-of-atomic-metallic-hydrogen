!*************************************************************************
       program ws   ! writed by cgzhang@theory.issp.ac.cn @ 2025.06.18
**************************************************************************
      implicit none
      real*8       autotemp1,autotemp2,autotemp3
     >            ,autotemp7,autotemp8,autotemp9
      integer*4    autotemp4,autotemp5,autotemp6
      real*8       b011,b022,b033,alatt
     >            ,deltax,deltay
     >            ,t2_tmp_big,t2_tmp_sml
     >            ,xs,ys,zs,xt2,yt2,zt2,t2
     >            ,distance_tmp, rc

      integer*4    V_n,V_n_1,V_n_0 
     >            ,J_n
     >            ,I_n,I_n_1,I_n_0
     >            ,n2,n3,n4,n5,n6,n7,n8,n9,n10
     >            ,n11,n12,n13,n14
     >            ,P_n, IV
     >            ,id_tmp_big, id_tmp_sml, id_tmp_big_j, id_tmp_big_k
     >            ,k1,k2,k3
     >            ,id_tmp

! note the array size
! at most 5000 interst and vacancies
! and J_n should big more than million
      real*8       V_x(1:10000000),V_y(1:10000000),V_z(1:10000000)
     >            ,I_x(1:10000000),I_y(1:10000000),I_z(1:10000000)
     >            ,b0(3,3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8    hmeps,dtemp
      real*8    fnlcx, fnlcy, fnlcz
      integer*4 ix, iy, iz,jx,jy,jz,kc,error
      integer*4 ip,jg,jat,jat_temp,jc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer*4 i,j,ii,m,k
     >         ,findex,fnum,itcum,step  !fnum is the number of cascade files need be analysed
      integer*4  nlcx,nlcy,nlcz,nlc
      integer*4, allocatable ::  ltop(:),linkmp(:)
     >                          ,nws(:),wsid(:,:)
      integer*4  nn1,mstart,max1,nnbrs,icount,inum
      real*8  ,allocatable :: xx0(:),yy0(:),zz0(:)  !xx0(i),yy0(i),zz0(i) are base configeration possitions
     >                        ,x0(:),y0(:),z0(:)
     >                        ,wsx0(:,:),wsy0(:,:),wsz0(:,:)
     >                        ,distance(:)

      integer*4 ,allocatable :: cid(:),id_distance(:)
      character*80 filename,chartemp,string1,string2,string
      real*8     sx,sy,sz,s,x0t,y0t,z0t
      integer*4  nix(27),niy(27),niz(27) 
      logical*4  pbcx,pbcy,pbcz
      real*8     rijsq,dx12t,dy12t,dz12t,rijsq_temp
      real*8     xl, yl, zl, xh, yh, zh
!**********************************************************************
      parameter  ( fnum=200              ) ! cascade file number
      parameter  ( step=1000             ) ! every n steps output config
      parameter  ( alatt=1.11979*2.0     )
      parameter  ( b011=13.4375          )
      parameter  ( b022=13.4375          )
      parameter  ( b033=12.952           )

      parameter  ( nlcx = int(b011/alatt) )  !b011/alatt
      parameter  ( nlcy = int(b022/alatt) )  !b022/alatt
      parameter  ( nlcz = int(b033/alatt) )  !b033/alatt
      parameter  ( nlc =nlcx*nlcy*nlcz    )          


      parameter  ( pbcx=.true.         )
      parameter  ( pbcy=.true.         )
      parameter  ( pbcz=.true.         )
      parameter  ( nnbrs=500           )
      parameter  ( hmeps = -1e-9       )
      parameter  ( s=1.0d0)


      real*8     x0v(nnbrs),y0v(nnbrs),z0v(nnbrs)
      integer*4  idv(nnbrs)
!
! link cell map
!
        data nix/ 0,-1,-1,-1,0,0,-1,1,-1, 0, 1,-1,0,1, 0,1,1, 1, 0,
     &            0, 0, 1, 1, 1,-1,-1,-1/
        data niy/ 0, 0,-1, 1,1,0, 0,0,-1,-1,-1, 1,1,1,-1,0,1,-1, 1,
     &            0,-1, 0,-1, 1, 0,-1, 1/
        data niz/ 0, 0, 0, 0,0,1, 1,1, 1, 1, 1, 1,1,1, 0,0,0, 0,-1,
     &           -1,-1,-1,-1,-1,-1,-1,-1/

      allocate(ltop(1:nlc),stat=error)  
      b0=0.0
      b0(1,1)=b011
      b0(2,2)=b022
      b0(3,3)=b033
!
!   initialize (clear) link cells
!
      do 10 i = 1, nlc
        ltop(i) = 0
 10   continue
      fnlcx = dfloat (nlcx)
      fnlcy = dfloat (nlcy)
      fnlcz = dfloat (nlcz)
c
c    read base configeration into xx0 yy0 zz0
c
      filename="bulk.cfg"
      open(unit=101,file=filename,form="FORMATTED",status="OLD")

        read(101,'(a21,i10)') string,nn1   !Number of particles = 364470
        read(101,'(a35)')    string       !A = 1 Angstrom (basic length-scale)
        do i=1,3
          do j=1,3
           read(101,*) string
          enddo
        enddo
        read(101,'(a13)')  string         !.NO_VELOCITY. 
        read(101,'(a13,i2)') string,inum  !entry_count = 9

        do i=1,inum-3
          read(101,*) string 
        enddo
      allocate(linkmp(1:nn1),stat=error)
! base isite
      allocate(xx0(1:nn1),stat=error)
      allocate(yy0(1:nn1),stat=error)
      allocate(zz0(1:nn1),stat=error)
      allocate(nws(1:nn1),stat=error)
!
! cascade ws atoms in each ws cell at most 4 atoms
      allocate(wsid(1:nn1,1:15),stat=error )
      allocate(wsx0(1:nn1,1:15),stat=error )
      allocate(wsy0(1:nn1,1:15),stat=error )
      allocate(wsz0(1:nn1,1:15),stat=error )
! configer analysed       
      allocate(x0(1:nn1),stat=error)
      allocate(y0(1:nn1),stat=error)
      allocate(z0(1:nn1),stat=error)
      allocate(cid(1:nn1),stat=error)

      if(error.ne.0) then
         stop "acclocate wrong!"
      endif
c
c     initialize parameters about xx0(i),yy0(i),zz0(i)
c
      do 11 i=1,nn1
         linkmp(i)=0
         xx0(i)=0.0d0
         yy0(i)=0.0d0
         zz0(i)=0.0d0
         x0(i)=0.0d0
         y0(i)=0.0d0
         z0(i)=0.0d0
11    continue
      do i = 1, nn1
          read(101,*) string 
          read(101,*) string 

         read(101,*) autotemp1,autotemp2,autotemp3
          xx0( i )=autotemp1  
          yy0( i )=autotemp2   ! autotemp1-3 is xx0,yy0,zz0
          zz0( i )=autotemp3
      enddo
      close(101)
! remap  base site configeration  xx0,yy0,zz0 in (0.0,1.0)
      do i=1,nn1
        if(xx0(i).gt.1.0d0) xx0(i)=xx0(i)-1.0d0
        if(xx0(i).lt.0.0d0) xx0(i)=xx0(i)+1.0d0

        if(yy0(i).gt.1.0d0) yy0(i)=yy0(i)-1.0d0
        if(yy0(i).lt.0.0d0) yy0(i)=yy0(i)+1.0d0

        if(zz0(i).gt.1.0d0) zz0(i)=zz0(i)-1.0d0
        if(zz0(i).lt.0.0d0) zz0(i)=zz0(i)+1.0d0
      enddo

! read base configer right
      do 12 i = 1, nn1
c
c    xx0(i),yy0(i),zz0(i) in (0.0,1.0),ix,iy,iz  begin from 0~~~~~(nlcx-1,nlxy-1,nlcz-1)
c    hmeps=-1.0e-9
c
          ix = int( ( xx0(i) + hmeps ) * fnlcx )
          iy = int( ( yy0(i) + hmeps ) * fnlcy )
          iz = int( ( zz0(i) + hmeps ) * fnlcz )
c          print*,'ix=',ix,'iy=',iy,'iz=',iz,i

          ip = 1 + ix + nlcx * (iy + nlcy*iz)       ! ip begin from 1~~~nlc
          if(ip.gt.nlc.or.ip.lt.1) then
            print*,'Error in ip',ip,nlc,i
            print*,'xx0(i)=',xx0(i),'yy0(i)=',yy0(i),'zz0(i)=',zz0(i)
          endif
c
c     assign atom i to link cell ip
c
          j = ltop(ip)
          ltop(ip) = i
          linkmp(i) = j
12    continue


!      write(*,*) "hello 0"

! open DefectVsTime_W.dat
      open(unit=77,file="DefectVsTime_W.dat",form="FORMATTED"
     >                                                ,status="REPLACE")
! open WSVsTime_W.dat
      open(unit=88,file="WSVsTime_W.dat",form="FORMATTED"
     >                                                ,status="REPLACE")
      do 23  findex=1,fnum                   ! do loop over every cascade file
         itcum=step*(findex-1)
!         write(*,*) "findex=?",findex

         nws=0
c        
c    open file named auto~itcum.cfg      
c     
         write(chartemp,'(I9)') itcum
c
c   autoitcum.cfg
c
         filename="voids-relax."//trim(adjustl(chartemp))//".cfg"
        open(102+findex,file=filename,form="FORMATTED",status="OLD")
c
c     read data in front of the cascade configeration
c
        read(102+findex,'(a21,i10)') string,icount   !
!        write(*,*) "icount",icount
        read(102+findex,'(a37)') string
        do i=1,3
          do j=1,3
            read(102+findex,*) string
          enddo
        enddo
        read(102+findex,'(a13)')  string         !.NO_VELOCITY. 
        read(102+findex,'(a13,i2)') string,inum  !entry_count = 7

        do i=1,inum-3
          read(102+findex,*) string
        enddo
!        write(*,*) "icount=",icount, "inum:",inum
c
c     read cascade configeration into-->x0(ii),y0(ii),z0(ii)
c
       do ii=1,icount
          read(102+findex,*) string
          read(102+findex,*) string

       read(102+findex,*)   autotemp1,autotemp2,autotemp3
          x0( ii )=autotemp1
          y0( ii )=autotemp2
          z0( ii )=autotemp3
      enddo
           close(102+findex)     ! close cascade file


!      write(*,*) "hello 1"

! remap cascade configeration 
      do i=1,icount
        if(x0(ii).gt.1.0d0) x0(ii)=x0(ii)-1.0d0
        if(x0(ii).lt.0.0d0) x0(ii)=x0(ii)+1.0d0

        if(y0(ii).gt.1.0d0) y0(ii)=y0(ii)-1.0d0
        if(y0(ii).lt.0.0d0) y0(ii)=y0(ii)+1.0d0

        if(z0(ii).gt.1.0d0) z0(ii)=z0(ii)-1.0d0
        if(z0(ii).lt.0.0d0) z0(ii)=z0(ii)+1.0d0
      enddo

!      write(*,*) "hello 2"
      do 13 ii=1,icount   ! ii is cascade configer id,x0(ii),y0(ii),z0(ii)
c       
c    x0(i),y0(i),z0(i) in which linkcell,ix,iy,iz begin from 0 to (nlcx-1),(nlcy-1),(nlcz-1)
c      
          ix = int( ( x0(ii) + hmeps ) * fnlcx )
          iy = int( ( y0(ii) + hmeps ) * fnlcy )
          iz = int( ( z0(ii) + hmeps ) * fnlcz )

          ip = 1 + ix + nlcx * (iy + nlcy*iz)  ! ip begin from 1~~nlc        
          if(ip.gt.nlc) then
            print*,'Error in ip',ip,nlc,ii
            print*,'x0(ii)=',x0(ii),'y0(ii)=',y0(ii),'z0(ii)=',z0(ii)
          endif
          
         i=ltop(ip)                 ! i=ltop(ip)=0 indicate have no atom, we must think about neigbour linkcell
         if (i.eq.0) then 
             m=0
            goto 26       ! i will not be 0
         endif
         m = 0
14       m = m+1
         idv(m) = i       ! store base id number
         x0v(m)  = xx0(i) ! store base x,y,z
         y0v(m)  = yy0(i)
         z0v(m)  = zz0(i)

         i        = linkmp(i)
         if (i.gt.0) goto 14
26         mstart = m
!
!**** Now select all particles from neighboring cells
!
!      write(*,*) "hello 3"
         do 15 kc = 2,27

           sx = 0.d0
           sy = 0.d0
           sz = 0.d0
           
          jx = ix + nix(kc)
          jy = iy + niy(kc)
          jz = iz + niz(kc)
!
!**** Minimum image conversion   ! if no periodicity,not count this (jx,jy,jz)
!
          if ((ix.eq.nlcx-1).and.(jx.gt.ix)) then
          jx = 0
          sx = s
          if(pbcx.eqv..false.)  goto 15
          elseif ((ix.eq.0).and.(jx.lt.ix)) then
          jx = nlcx-1
          sx = -s
          if(pbcx.eqv..false.)  goto 15
          endif

          if ((iy.eq.nlcy-1).and.(jy.gt.iy)) then
          jy = 0
          sy = s
          if(pbcy.eqv..false.)  goto 15
          elseif ((iy.eq.0).and.(jy.lt.iy)) then
          jy = nlcy-1
          sy = -s
          if(pbcy.eqv..false.)  goto 15
          endif
          if ((iz.eq.nlcz-1).and.(jz.gt.iz)) then
          jz = 0
          sz = s
          if(pbcz.eqv..false.)  goto 15
          elseif ((iz.eq.0).and.(jz.lt.iz)) then
          jz = nlcz-1
          sz = -s
          if(pbcz.eqv..false.) goto 15
          endif
c
c**** Index of neighbouring cell
c
             jc = 1+jx + nlcx*( jy + nlcy*jz )   !j
             j  = ltop(jc)
c
c**** Bypass this neighbouring cell if it is empty
c
             if (j.eq.0) goto 15            ! if no atoms in this linkcell jc, then next linkcell 2~~27
17          m = m+1
             idv(m) = j
             x0v(m) = xx0(j) + sx
             y0v(m) = yy0(j) + sy
             z0v(m) = zz0(j) + sz
             j = linkmp(j)
             if (j.gt.0) goto 17
!
!**** Save number of particles in first half of link cells
!
15     continue
!      write(*,*) "hello 4"
!
!**** We have now found all the neighbouring particles of cell ic.
         max1= m
         if (max1.gt.nnbrs) then
             write(*,*) 'max1.gt.nnbrs something wrong'
             stop
         endif

         do 18 jg=1,max1          ! x0v (0.0,1.0)
          x0v(jg) = b011*x0v(jg)
          y0v(jg) = b022*y0v(jg)
          z0v(jg) = b033*z0v(jg)
18      continue
!      write(*,*) "hello 5"
c
c      atom ii  x0(ii),y0(ii),z0(ii) in cascade configeration,convert into real coordinatation
c           
          x0t = b011*x0(ii)  ! x0(ii) (0.0,1.0)
          y0t = b022*y0(ii)
          z0t = b033*z0(ii)
         
           rijsq_temp=10000.0d0     !  as large as possible
c      
c    link cell the ideal config,and find which ideal configeration site is closest to this atom ii 
c
           do 20 m=1,max1
              jat=idv(m)
              dx12t=x0v(m)-x0t
              dy12t=y0v(m)-y0t
              dz12t=z0v(m)-z0t
              rijsq=dx12t*dx12t+dy12t*dy12t+dz12t*dz12t
             if(rijsq_temp.gt.rijsq) then
                jat_temp=jat
                rijsq_temp=rijsq
             endif
20       continue
      
!      write(*,*) "hello 6"
!      write(*,*) "ii:",ii,"jat_temp:", jat_temp
c
c   WS methods onely if the closest then have
c
        nws(jat_temp)=nws(jat_temp)+1 

        wsid( jat_temp,nws(jat_temp) )=ii
        wsx0(jat_temp,nws(jat_temp))=x0(ii)
        wsy0(jat_temp,nws(jat_temp))=y0(ii)
        wsz0(jat_temp,nws(jat_temp))=z0(ii)
!      write(*,*) "hello 7"
13    continue


!      write(*,*) "hello 8"
! nws(1:nn1),wsid(1:nn1,1:nws(1:nn1)),wsx0,wsy0,wsz0(1:nn1,nws(1:nn1))
       V_n=0
       I_n=0
       J_n=0

       n2=0
       n3=0
       n4=0

       do i=1,nn1
         if( nws(i).eq.0 ) then
               V_n=V_n+1
               V_x(V_n)=xx0(i)  ! V_x,V_y,V_z in (0,1)
               V_y(V_n)=yy0(i)
               V_z(V_n)=zz0(i)
         elseif(nws(i).eq.1) then
               J_n=J_n+1
         elseif(nws(i).eq.2) then
               n2=n2+1
               I_n=I_n+1
               I_x(I_n)=wsx0(i,1) 
               I_y(I_n)=wsy0(i,1)
               I_z(I_n)=wsz0(i,1)
         elseif(nws(i).eq.3) then
               n3=n3+1   ! only add one number
               do j=1,2
                 I_n=I_n+1
                 I_x(I_n)=wsx0(i,j)
                 I_y(I_n)=wsy0(i,j)
                 I_z(I_n)=wsz0(i,j)
               enddo
         elseif(nws(i).eq.4) then
               n4=n4+1   ! only add one number
               do j=1,3
                 I_n=I_n+1
                 I_x(I_n)=wsx0(i,j)
                 I_y(I_n)=wsy0(i,j)
                 I_z(I_n)=wsz0(i,j)
               enddo
         else
               write(*,*) "warning !!!!!! nws(i)",nws(i)
         endif
       enddo
         if( findex.eq.1 ) then
          write(77,*) "--(1)step--(2)V_n--(3) I_n--"

        write(88,*) "-(1)step---(2)V_n---(3)J_n-----"
     >             , "(4)n2----(5)n3---(6)n4"
         endif
         write(77,777) itcum,V_n,I_n,V_n-I_n
777      format (i7,1x,i5,5x,i5,5x,i5)

         write(88,888) itcum,V_n,J_n,n2,n3,n4
888      format (i7,1x,i5,5x,i7,i5,2x,i5,2x,i5)

! detect print
         if(V_n+J_n+n2+n3+n4.ne.nn1) then
           write(*,*) "warning !!!! V_n+J_n+n2+n3+n4.ne.nn1"
         endif
!         if( V_n-4.ne.I_n) then
!            write(*,*) "warning !!!!!V_n-1.ne.I_n"
!         endif

!**************************************************************************
!
!  open vac.${step}.cfg
!
        filename="vac."//trim(adjustl(chartemp))//".cfg"
       open(unit=1000+findex,file=filename,form="FORMATTED"
     >      ,status="REPLACE")
        inum=3 
        write(1000+findex,'(a,i7)')"Number of particles =",V_n
        write(1000+findex,'(a)')"A = 1.0 Angstrom (basic length-scale)"
        do i=1,3
          do j=1,3
      write(1000+findex,'(a,i1,a,i1,a,f16.8,1x,a)')"H0(",i,",",j,") =",
     >            b0(j,i),"A"
          enddo
        enddo
        write(1000+findex,'(a)')  ".NO_VELOCITY."
        write(1000+findex,'(a,i2)') "entry_count = ",inum
        write(1000+findex,'(a)')"1.0"
        write(1000+findex,'(a)')"H"
        do i=1,V_n
          write(1000+findex,*) V_x(i),V_y(i),V_z(i)
        enddo
        close(1000+findex)
!
! open int.${step}.cfg
!
       filename="int."//trim(adjustl(chartemp))//".cfg"
       open(unit=2000+findex,file=filename,form="FORMATTED"
     >      ,status="REPLACE")
        inum=3
        write(2000+findex,'(a,i7)')"Number of particles =",I_n
        write(2000+findex,'(a)')"A = 1.0 Angstrom (basic length-scale)"
        do i=1,3
          do j=1,3
      write(2000+findex,'(a,i1,a,i1,a,f16.8,1x,a)')"H0(",i,",",j,") =",
     >            b0(j,i),"A"
          enddo
        enddo
        write(2000+findex,'(a)')  ".NO_VELOCITY."
        write(2000+findex,'(a,i2)') "entry_count = ",inum
        write(2000+findex,'(a)')"1.0"
        write(2000+findex,'(a)')"H"
        do i=1,I_n
          write(2000+findex,*) I_x(i),I_y(i),I_z(i)
        enddo
        close(2000+findex)
!
! open defect.${Step}.cfg 
!
       filename="defect."//trim(adjustl(chartemp))//".cfg"
       open(unit=3000+findex,file=filename,form="FORMATTED"
     >      ,status="REPLACE")
        inum=4
        write(3000+findex,'(a,i7)')"Number of particles =",I_n+V_n
        write(3000+findex,'(a)')"A = 1.0 Angstrom (basic length-scale)"
        do i=1,3
          do j=1,3
      write(3000+findex,'(a,i1,a,i1,a,f16.8,1x,a)')"H0(",i,",",j,") =",
     >            b0(j,i),"A"
          enddo
        enddo
        write(3000+findex,'(a)')  ".NO_VELOCITY."
        write(3000+findex,'(a,i2)') "entry_count = ",inum
        write(3000+findex,'(a)') "auxiliary[0] = IV"
        write(3000+findex,'(a)')"1.0"
        write(3000+findex,'(a)')"H"

        IV=1
        do i=1,I_n
          write(3000+findex,*) I_x(i),I_y(i),I_z(i),IV
        enddo

        IV=0
        do i=1,V_n
          write(3000+findex,*) V_x(i),V_y(i),V_z(i),IV
        enddo

        close(3000+findex)

!        write(*,*) "findex end good!", findex
!
!****************************************************************************
! get the pure v
! V_x(1:V_n), V_y(1:V_n), V_z(1:V_n)
! I_x(1:I_n), I_y(1:I_n), I_z(1:I_n)

      allocate(distance(1:V_n),stat=error)
      allocate(id_distance(1:V_n),stat=error)
      if(error.ne.0) then
         stop "acclocate wrong!"
      endif

       do i=1,V_n
         xs=V_x(i)
         ys=V_y(i)
         zs=V_z(i)
         t2_tmp_sml=1000000.0
         id_tmp_sml=0
         do j=1,I_n
           do k1=-1,1
             xt2=(xs-I_x(j)+k1)*(xs-I_x(j)+k1)*b011*b011
           do k2=-1,1
             yt2=(ys-I_y(j)+k2)*(ys-I_y(j)+k2)*b022*b022
           do k3=-1,1
             zt2=(zs-I_z(j)+k3)*(zs-I_z(j)+k3)*b033*b033
             t2=xt2+yt2+zt2
             if (t2.lt.t2_tmp_sml) then
               t2_tmp_sml=t2
               id_tmp_sml=j
             endif
           enddo
           enddo
           enddo
         enddo
         distance(i)=t2_tmp_sml
         id_distance(i)=i
       enddo

        do i = 1, V_n-1
            do j = 1, V_n-i
                if (distance(j) < distance(j + 1)) then
                    distance_tmp = distance(j)
                    distance(j) = distance(j + 1)
                    distance(j + 1) = distance_tmp

                    id_tmp=id_distance(j)
                    id_distance(j)=id_distance(j+1)
                    id_distance(j+1)=id_tmp
                endif
            enddo
        enddo

!**************************************************************************
!
!  open vac-pure.${step}.cfg
!
        filename="vac-pure."//trim(adjustl(chartemp))//".cfg"
       open(unit=4000+findex,file=filename,form="FORMATTED"
     >      ,status="REPLACE")
        inum=3
        write(4000+findex,'(a,i7)')"Number of particles =",V_n-I_n
        write(4000+findex,'(a)')"A = 1.0 Angstrom (basic length-scale)"
        do i=1,3
          do j=1,3
      write(4000+findex,'(a,i1,a,i1,a,f16.8,1x,a)')"H0(",i,",",j,") =",
     >            b0(j,i),"A"
          enddo
        enddo
        write(4000+findex,'(a)')  ".NO_VELOCITY."
        write(4000+findex,'(a,i2)') "entry_count = ",inum
        write(4000+findex,'(a)')"1.0"
        write(4000+findex,'(a)')"H"

       do j=1,V_n-I_n
        id_tmp_big=id_distance(j)
        write(4000+findex,*) V_x(id_tmp_big),
     >                       V_y(id_tmp_big),V_z(id_tmp_big)
       enddo
        close(4000+findex)

! get the dissociation time, when the closest distance among vacancies is greater than alatt*2, alatt=1.11979 A
       t2_tmp_sml=1000000.0
       do j=1,V_n-I_n-1
         id_tmp_big_j=id_distance(j)
         xs=V_x(id_tmp_big_j)
         ys=V_y(id_tmp_big_j)
         zs=V_z(id_tmp_big_j)
         do k=j+1,V_n-I_n
           id_tmp_big_k=id_distance(k)
           do k1=-1,1
       xt2=(xs-V_x(id_tmp_big_k)+k1)*(xs-V_x(id_tmp_big_k)+k1)*b011*b011
           do k2=-1,1
       yt2=(ys-V_y(id_tmp_big_k)+k2)*(ys-V_y(id_tmp_big_k)+k2)*b022*b022
           do k3=-1,1
       zt2=(zs-V_z(id_tmp_big_k)+k3)*(zs-V_z(id_tmp_big_k)+k3)*b033*b033
             t2=xt2+yt2+zt2
             if (t2.lt.t2_tmp_sml) then
               t2_tmp_sml=t2
             endif
           enddo
           enddo
           enddo
         enddo 
       enddo

! output the dissociation time
       rc=1.11979*3.0
       if (t2_tmp_sml.gt.rc*rc) then
         write(*,*) "step: ", itcum,
     > "dissociation_time: ", itcum*0.2, " fs,", 
     >  " distance: ", sqrt(t2_tmp_sml)
       endif

       deallocate(distance)
       deallocate(id_distance)
!****************************************************************************

23      continue  ! findex loop
        close(66)
        close(77)
        close(88)

      deallocate(xx0)
      deallocate(yy0)
      deallocate(zz0)

      deallocate(x0)
      deallocate(y0)
      deallocate(z0)

      deallocate(nws) 
      deallocate(wsid)
      deallocate(wsx0)
      deallocate(wsy0)
      deallocate(wsz0)


      deallocate(ltop)
      deallocate(linkmp)

      end program
