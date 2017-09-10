!c================================c
        subroutine ipd0fm
!c================================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.gd2'
        include 'table.gd3'
        include 'table.gd4'
        include 'table.gd5'
        include 'table.gd6'
        include 'table.gd7'
!c----------------------------------------------------------------------c
        strh=1.0;reno=200;dt=0.01
        m11=258;jb1=66;je1=194;n11=258;kb1=66;ke1=194     
        m12=6;jb2=3;je2=5;n12=6;kb2=3;ke2=5     
        m13=10;jb3=4;je3=8;n13=10;kb3=4;ke3=8     
        m14=18;jb4=6;je4=14;n14=18;kb4=6;ke4=14     
        m15=34;jb5=10;je5=26;n15=34;kb5=10;ke5=26     
        m16=66;jb6=18;je6=50;n16=66;kb6=18;ke6=50
        m17=130;jb7=34;je7=98;n17=130;kb7=34;ke7=98          
        return
        end
!c================================c
        subroutine ipd0uf
!c================================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.gd2'
        include 'table.gd3'
        include 'table.gd4'
        include 'table.gd5'
        include 'table.gd6'
        include 'table.gd7'
!c----------------------------------------------------------------------c
        open(unit=ioc,file=dt0ufm,form='unformatted',status='old')
        read(ioc) strh,reno
        read(ioc) dt
        read(ioc) l11
        read(ioc) m11,n11,jb1,je1,kb1,ke1
        if((levgrd-1).eq.0) goto 1001
        read(ioc) m12,n12,jb2,je2,kb2,ke2
        if((levgrd-2).eq.0) goto 1001
        read(ioc) m13,n13,jb3,je3,kb3,ke3
        if((levgrd-3).eq.0) goto 1001
        read(ioc) m14,n14,jb4,je4,kb4,ke4
        if((levgrd-4).eq.0) goto 1001
        read(ioc) m15,n15,jb5,je5,kb5,ke5
        if((levgrd-5).eq.0) goto 1001
        read(ioc) m16,n16,jb6,je6,kb6,ke6
        if((levgrd-6).eq.0) goto 1001
        read(ioc) m17,n17,jb7,je7,kb7,ke7
        if((levgrd-7).eq.0) goto 1001
1001    close(ioc)
        return
        end

!c================================c
        subroutine ipd1fm
!c================================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c--------------------------------c
        open(unit=ioc,file=dt1fmt,form='formatted',status='old')
        read(ioc,101)
101     format('*',25(1h-),'momentum grid on level one',25(1h-),'*')
        call rdcfmt(yv1,' yv1= ',ioc,2,m11,nj1)
        call rdcfmt(zw1,' zw1= ',ioc,2,n11,nk1)
        close(unit=ioc)
        return
        end
!c===============================c
        subroutine ipd1uf
!c===============================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
!c--------------------------c
        open(19,file='..\input\xy.dat')
        read(19,*)yv1
        close(19)
        open(19,file='..\input\xy.dat')
        read(19,*)zw1
        close(19)
        return
        end



!c===============================c
        subroutine ipdff0
!c===============================c
        use varalc 
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
        include 'table.alc'
!c-------------------------------c
        do j=1,nj1
            do k=1,nk1
                vi1(j,k)=0
                wi1(j,k)=0
            enddo
        enddo
        
        open(19,file='..\input\vav.dat')
        read(19,*)vi1
        close(19)
        open(19,file='..\input\wav.dat')
        read(19,*)wi1
        close(19)
!        open(19,file='..\input\pav.dat')
!        read(19,*)pre
!        close(19)
        
        return
        end

!c================================c
        subroutine opd0fm
!c================================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.gd2'
        include 'table.gd3'
        include 'table.gd4'
        include 'table.gd5'
        include 'table.gd6'
        include 'table.gd7'
!c----------------------------------------------------------------------c
20      format(' *',6x,'stralh    number    strh=',1pd10.3,6x,'*')
21      format(' *',6x,'renold    number    reno=',1pd10.3,6x,'*')

40      format(' *',4x,'the time marching step   dt: dt=',1pd9.2,2x,'*')
102     format(' *',4x,'level  1  point number in X:l11=',i6,5x,'*')
103     format(' *',4x,'level 1: m11',i4,4x,'jb1',i4,4x,'je1',i4,5x,'*')
104     format(' *',4x,'level 1: n11',i4,4x,'kb1',i4,4x,'ke1',i4,5x,'*')
108     format(' *',4x,'level 2: m12',i4,4x,'jb2',i4,4x,'je2',i4,5x,'*')
109     format(' *',4x,'level 2: n12',i4,4x,'kb2',i4,4x,'ke2',i4,5x,'*')
111     format(' *',4x,'level 3: m13',i4,4x,'jb3',i4,4x,'je3',i4,5x,'*')
112     format(' *',4x,'level 3: n13',i4,4x,'kb3',i4,4x,'ke3',i4,5x,'*')
114     format(' *',4x,'level 4: m14',i4,4x,'jb4',i4,4x,'je4',i4,5x,'*')
115     format(' *',4x,'level 4: n14',i4,4x,'kb4',i4,4x,'ke4',i4,5x,'*')
117     format(' *',4x,'level 5: m15',i4,4x,'jb5',i4,4x,'je5',i4,5x,'*')
118     format(' *',4x,'level 5: n15',i4,4x,'kb5',i4,4x,'ke5',i4,5x,'*')
120     format(' *',4x,'level 6: m16',i4,4x,'jb6',i4,4x,'je6',i4,5x,'*')
121     format(' *',4x,'level 6: n16',i4,4x,'kb6',i4,4x,'ke6',i4,5x,'*')
50      format(' *                                               *')
52      format(' *             The Criterion Numbers             *')
53      format(' *           The Coordinate Parameters           *')
54      format(' *-----------------------------------------------*')
!c----------------------------------------------------------------------c
        open(unit=ioc,file='dat0.fmt',form='formatted')
        write(ioc,54)
        write(ioc,50)
        write(ioc,52)
        write(ioc,50)
        write(ioc,20) strh
        write(ioc,21) reno
        write(ioc,50)
        write(ioc,53)
        write(ioc,50)
        write(ioc,40) dt
        write(ioc,102) l11
        write(ioc,103) m11,jb1,je1
        write(ioc,104) n11,kb1,ke1
        if((levgrd-1).eq.0) goto 1001
        write(ioc,108) m12,jb2,je2
        write(ioc,109) n12,kb2,ke2
        if((levgrd-2).eq.0) goto 1001
        write(ioc,111) m13,jb3,je3
        write(ioc,112) n13,kb3,ke3
        if((levgrd-3).eq.0) goto 1001
        write(ioc,114) m14,jb4,je4
        write(ioc,115) n14,kb4,ke4
        if((levgrd-4).eq.0) goto 1001
        write(ioc,117) m15,jb5,je5
        write(ioc,118) n15,kb5,ke5
        if((levgrd-5).eq.0) goto 1001
        write(ioc,120) m16,jb6,je6
        write(ioc,121) n16,kb6,ke6
        if((levgrd-6).eq.0) goto 1001
1001    write(ioc,50)
        write(ioc,54)
        close(unit=ioc)
        return
        end

!c================================c
        subroutine opd1fm
!c================================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c--------------------------------c
        open(unit=ioc,file='dat1.fmt',form='formatted')
        write(ioc,101)
101     format('*',25(1h-),'momentum grid on level one',25(1h-),'*')
        call wtcfmt(yv1,' yv1= ',ioc,2,m11,nj1)
        call wtcfmt(zw1,' zw1= ',ioc,2,n11,nk1)
        close(unit=ioc)
        return
        end

!c================================c
        subroutine opdff0
!c================================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c--------------------------------c
        open(unit=ioc,file='data.ff0',form='unformatted')
        call wtfufm(uu1,ioc,1,1,m11,n11,nj1,nk1)
        call wtfufm(vv1,ioc,2,1,m11,n11,nj1,nk1)
        call wtfufm(ww1,ioc,1,2,m11,n11,nj1,nk1)
        call wtfufm(pre,ioc,1,1,m11,n11,nj1,nk1)
        close(unit=ioc)
        return
        end
!c================================c
        subroutine opdhf0
!c================================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c--------------------------------c
        open(unit=ioc,file='data.hf0',form='unformatted')
        call wtfufm(hu0,ioc,2,2,m11-1,n11-1,nj1,nk1)
        call wtfufm(hv0,ioc,3,2,m11-1,n11-1,nj1,nk1)
        call wtfufm(hw0,ioc,2,3,m11-1,n11-1,nj1,nk1)
        close(unit=ioc)
        return
        end

        subroutine wtcfmt(x,cx,ioc,i0,i1,ni)
        include 'table.prc'
        dimension x(ni)
        character cx*6
!c------------------------------------------------c
100     format()
101     format('*',76(1h-),'*')
111     format(2x,' i =',4x,5(i4,8x),i4)
112     format(a6,1p6d12.4)
!c------------------------------------------------c
        write(ioc,100)
        write(ioc,101)
        iend=i0-1
10      if(iend.eq.i1) goto 11
        ibeg=iend+1
        iend=iend+6
        iend=min0(iend,i1)
        write(ioc,100)
        write(ioc,111)(i,i=ibeg,iend)
        write(ioc,112) cx,(x(i),i=ibeg,iend)
        goto 10
11      write(ioc,101)
        return
        end
        subroutine rdcfmt(x,cx,ioc,i0,i1,ni)
        include 'table.prc'
        dimension x(ni)
        character cx*6
!c------------------------------------------------c
100     format()
101     format('*',76(1h-),'*')
111     format(2x,'i =',4x,5(i4,8x),i4)
112     format(a6,1p6d12.4)
!c------------------------------------------------c
        read(ioc,100)
        read(ioc,101)
        iend=i0-1
10      if(iend.eq.i1) goto 11
        ibeg=iend+1
        iend=iend+6
        iend=min0(iend,i1)
        read(ioc,100)
        read(ioc,111)(iwrk,i=ibeg,iend)
        read(ioc,112) cx,(x(i),i=ibeg,iend)
        goto 10
11      read(ioc,101)
        return
        end

        subroutine wtcufm(x,ioc,i0,i1,ni)
        include 'table.prc'
        dimension x(ni)
        write(ioc)(x(i),i=i0,i1)
        return
        end
        subroutine rdcufm(x,ioc,i0,i1,ni)
        include 'table.prc'
        dimension x(ni)
        read(ioc)(x(i),i=i0,i1)
        return
        end

        subroutine wtffmt(f,mode,ioc,i0,j0,k0,i1,j1,k1,ni,nj,nk)
        include 'table.prc'
        dimension f(ni,nj,nk)
!c------------------------------------------c
100     format()
121     format(7x,'i =',i4)
122     format(7h    j =,i5,1x,6i10)
123     format(4h   k)
124     format(i4,1x,1p7d10.2)

131     format(7x,'j =',i4)
132     format(7h    k =,i5,1x,6i10)
133     format(4h   i)
134     format(i4,1x,1p7d10.2)

141     format(7x,'k =',i4)
142     format(7h    i =,i5,1x,6i10)
143     format(4h   j)
144     format(i4,1x,1p7d10.2)
!c------------------------------------------c
        ifl=i0+i1
        jfl=j0+j1
        kfl=k0+k1
        if(mode.eq.1) then
        do 11 i=i0,i1
        write(ioc,121) i
        write(ioc,100)
        jbeg=j0-7
12      jbeg=jbeg+7
        jend=jbeg+6
        jend=min0(jend,j1)
        do 13 kk=k0,k1
        k=kfl-kk
        write(ioc,124) k,(f(i,j,k),j=jbeg,jend)
13      continue
        write(ioc,123)
        write(ioc,122)(j,j=jbeg,jend)
        write(ioc,100)
        if(jend.lt.j1) goto 12
11      continue
        else if(mode.eq.2) then
        do 14 j=j0,j1
        write(ioc,131) j
        write(ioc,100)
        kbeg=k0-7
15      kbeg=kbeg+7
        kend=kbeg+6
        kend=min0(kend,k1)
        do 16 ii=i0,i1
        i=ifl-ii
        write(ioc,134) i,(f(i,j,k),k=kbeg,kend)
16      continue
        write(ioc,133)
        write(ioc,132)(k,k=kbeg,kend)
        write(ioc,100)
        if(kend.lt.k1) goto 15
14      continue
        else if(mode.eq.3) then
        do 17 k=k0,k1
        write(ioc,141) k
        write(ioc,100)
        ibeg=i0-7
18      ibeg=ibeg+7
        iend=ibeg+6
        iend=min0(iend,i1)
        do 19 jj=j0,j1
        j=jfl-jj
        write(ioc,144) j,(f(i,j,k),i=ibeg,iend)
19      continue
        write(ioc,143)
        write(ioc,142)(i,i=ibeg,iend)
        write(ioc,100)
        if(iend.lt.i1) goto 18
17      continue
        end if
        return
        end
        subroutine rdffmt(f,mode,ioc,i0,j0,k0,i1,j1,k1,ni,nj,nk)
        include 'table.prc'
        dimension f(ni,nj,nk)
!c------------------------------------------c
100     format()
121     format(7x,'i =',i4)
122     format(7h    j =,i5,1x,6i10)
123     format(4h   k)
124     format(i4,1x,1p7d10.2)

131     format(7x,'j =',i4)
132     format(7h    k =,i5,1x,6i10)
133     format(4h   i)
134     format(i4,1x,1p7d10.2)

141     format(7x,'k =',i4)
142     format(7h    i =,i5,1x,6i10)
143     format(4h   j)
144     format(i4,1x,1p7d10.2)
!c------------------------------------------c
        ifl=i0+i1
        jfl=j0+j1
        kfl=k0+k1
        if(mode.eq.1) then
        do 11 ii=i0,i1
        read(ioc,121) i
        read(ioc,100)
        jbeg=j0-7
12      jbeg=jbeg+7
        jend=jbeg+6
        jend=min0(jend,j1)
        do 13 kk=k0,k1
        read(ioc,124) k,(f(i,j,k),j=jbeg,jend)
13      continue
        read(ioc,123)
        read(ioc,122)(jwrk,j=jbeg,jend)
        read(ioc,100)
        if(jend.lt.j1) goto 12
11      continue
        end if
        if(mode.eq.2) then
        do 14 jj=j0,j1
        read(ioc,131) j
        read(ioc,100)
        kbeg=k0-7
15      kbeg=kbeg+7
        kend=kbeg+6
        kend=min0(kend,k1)
        do 16 ii=i0,i1
        read(ioc,134) i,(f(i,j,k),k=kbeg,kend)
16      continue
        read(ioc,133)
        read(ioc,132)(kwrk,k=kbeg,kend)
        read(ioc,100)
        if(kend.lt.k1) goto 15
14      continue
        end if
        if(mode.eq.3) then
        do 17 kk=k0,k1
        read(ioc,141) k
        read(ioc,100)
        ibeg=i0-7
18      ibeg=ibeg+7
        iend=ibeg+6
        iend=min0(iend,i1)
        do 19 jj=j0,j1
        read(ioc,144) j,(f(i,j,k),i=ibeg,iend)
19      continue
        read(ioc,143)
        read(ioc,142)(iwrk,i=ibeg,iend)
        read(ioc,100)
        if(iend.lt.i1) goto 18
17      continue
        end if
        return
        end

        subroutine wtfufm(f,ioc,j0,k0,j1,k1,nj,nk)
        include 'table.prc'
        dimension f(nj,nk)
        write(ioc) ((f(j,k),j=j0,j1),k=k0,k1)
        return
        end
        subroutine rdfufm(f,ioc,i0,j0,k0,i1,j1,k1,nj,nk)
        include 'table.prc'
        dimension f(nj,nk)
        read(ioc) ((f(j,k),j=j0,j1),k=k0,k1)
        return
        end

        subroutine chk3d1(ichk,f,fnm,nch,j0,k0,j1,k1,nj,nk)
        include 'table.prc'
        dimension f(nj,nk)
        character fnm*8
!c------------------------------------------c
100     format()
121     format(7x,'i =',i4)
122     format(7h    j =,i5,1x,6i10)
123     format(4h   k)
124     format(i4,1x,1p7d10.2)
!c------------------------------------------c
        jfl=j0+j1
        kfl=k0+k1
        if(i.eq.ichk) then
        open(unit=nch,file=fnm)
        write(nch,121) i
        write(nch,100)
        jbeg=j0-7
12      jbeg=jbeg+7
        jend=jbeg+6
        jend=min0(jend,j1)
        do 13 kk=k0,k1
        k=kfl-kk
        write(nch,124) k,(f(j,k),j=jbeg,jend)
13      continue
        write(nch,123)
        write(nch,122)(j,j=jbeg,jend)
        write(nch,100)
        if(jend.lt.j1) goto 12
        end if
        return
        end
        subroutine chk3d2(ichk,f,i,j0,k0,j1,k1,ni,nj,nk)
        include 'table.prc'
        dimension f(ni,nj,nk)
!c------------------------------------------c
100     format()
121     format(7x,'i =',i4)
122     format(7h    j =,i5,1x,6i10)
123     format(4h   k)
124     format(i4,1x,1p7d10.2)
!c------------------------------------------c
        jfl=j0+j1
        kfl=k0+k1
        if(i.eq.ichk) then
        write(*,121) i
        write(*,100)
        jbeg=j0-7
12      jbeg=jbeg+7
        jend=jbeg+6
        jend=min0(jend,j1)
        do 13 kk=k0,k1
        k=kfl-kk
        write(*,124) k,(f(i,j,k),j=jbeg,jend)
13      continue
        write(*,123)
        write(*,122)(j,j=jbeg,jend)
        write(*,100)
        if(jend.lt.j1) goto 12
        end if
        close(10)
        return
        end
!c*********************************c
!c*********************************c
        subroutine shuchu
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.lei'
	
	open(91,file='..\output\shuchu.dat',form='formatted')
	write(91,*) 'TITLE    ="Plot3D DataSet"'
        write(91,*) 'VARIABLES = "y" "z" "p" "v" "w"'
        write(91,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
        write(91,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
        write(91,*) 'ZONE T="Zone-original grid"'
        write(91,*) 'STRANDID=0, SOLUTIONTIME=0'
        write(91,*) 'I=258, J=258, K=1, ZONETYPE=Ordered'
        write(91,*) 'DATAPACKING=POINT'
        write(91,*) 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE )'
	do k=1,258
	    do j=1,258
	        write(91,1010)ym1(j),zm1(k),pre(j,k),vv1(j,k),ww1(j,k)
1010	        format(5e17.9)
            enddo
        enddo
	close(91)
	
	open(92,file='..\output\grid.p3d',form='unformatted')
        write(92) nj1,nk1
        write(92)((sngl(ym1(j)),j=1,nj1),k=1,nk1),((sngl(zm1(k)),j=1,nj1),k=1,nk1)
        close(92)
	
	open(93,file='..\output\vw.p3d',form='unformatted')
        write(93) nj1,nk1,2
        write(93)((sngl(0.5*(vv1(j,k)+vv1(j,k))),j=1,nj1),k=1,nk1),((sngl(0.5*(ww1(j,k)+ww1(j,k))),j=1,nj1),k=1,nk1)
        close(93)
	
        open(20,file='..\output\leinuo_dat.dat')
!        write(20,*) 'VARIABLES = "y" "z" "vpvp" "vpwp" "wpwp"'
!        write(20,*)'I=',258,', J=',258,', K=',1
        do k=1,258
	    do j=1,258
	        write(20,1111)ym1(j),zm1(k),vpvp(j,k),vpwp(j,k),wpwp(j,k)
1111	        format(5e17.9)
            enddo
        enddo
        close(20)
        
        return
        end
