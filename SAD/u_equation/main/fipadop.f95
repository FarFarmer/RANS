!================================c
        subroutine ipd0fm
!================================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.gd2'
        include 'table.gd3'
        include 'table.gd4'
        include 'table.gd5'
        include 'table.gd6'
        include 'table.gd7'
        
!----------------------------------------------------------------------c
        strh=1.0;reno=200;dt=0.05
        m11=258;jb1=66;je1=194;n11=258;kb1=66;ke1=194     
        m12=6;jb2=3;je2=5;n12=6;kb2=3;ke2=5     
        m13=10;jb3=4;je3=8;n13=10;kb3=4;ke3=8     
        m14=18;jb4=6;je4=14;n14=18;kb4=6;ke4=14     
        m15=34;jb5=10;je5=26;n15=34;kb5=10;ke5=26     
        m16=66;jb6=18;je6=50;n16=66;kb6=18;ke6=50
        m17=130;jb7=34;je7=98;n17=130;kb7=34;ke7=98          

        return
        end

!===============================!
        subroutine ipd1uf
!===============================!
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!-------------------------------!
        open(199,file='..\input\xy.dat')
        read(199,*)yv1
        close(199)
        open(190,file='..\input\xy.dat')
        read(190,*)zw1
        close(190)
        return
        end

!===============================!
        subroutine ipdff0u
!===============================!
        use varalc 
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!-------------------------------!
        include 'table.alc'
!-------------------------------!
        open(19,file='..\input\vav.dat')
        read(19,*) vi1
        close(19)
        open(19,file='..\input\wav.dat')
        read(19,*) wi1
        close(19)
        do j=1,nj1
            do k=1,nk1
                ui1(j,k)=0.0
            enddo
        enddo

        open(191,file='..\input\uvw.dat',form='formatted')
	do j=1,258
	    do k=1,258
	        read(191,1010)ym1(j),zm1(k),ui1(j,k),vi1(j,k),wi1(j,k)
1010	        format(5e17.9)
            enddo
        enddo
        close(191)
        
        upri(:,:)=ui1(:,:)
        
        return
        end
	  
	  

!*********************************c
	subroutine shuchu2(ncompt)
!================================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'	
	  
	open(91,file='..\output\uav.dat',form='formatted')
	write(91,*) 'TITLE    ="Plot3D DataSet"'
	write(91,*) 'VARIABLES = "y" "z" "u" "v" "w" "miu"'
	write(91,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
        write(91,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
        write(91,*) 'ZONE T="Zone-original grid"'
        write(91,*) 'STRANDID=0, SOLUTIONTIME=0'
        write(91,*) 'I=258, J=258, K=1, ZONETYPE=Ordered'
        write(91,*) 'DATAPACKING=POINT'
        write(91,*) 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE single )'
	do j=1,258
	    do k=1,258
	        write(91,1010) ym1(j),zm1(k),uu1(j,k),vv1(j,k),ww1(j,k),dmiu1212(j,k)
1010	        format(5e17.9)
	    enddo
        enddo
        close(91)
        
        open(191,file='..\input\uvw.dat',form='formatted')
	do j=1,258
	    do k=1,258
	        read(191,1011) yv1(j),zw1(k),upri(j,k),vi1(j,k),wi1(j,k)
1011	        format(5e17.9)
            enddo
        enddo
        close(191)
        
        open(111,file='..\output\ucomp.dat',form='formatted')
        write(111,*) 'TITLE    ="Plot3D DataSet"'
	write(111,*) 'VARIABLES = "y" "z" "umod" "upri"'
	write(111,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
        write(111,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
        write(111,*) 'ZONE T="Zone-original grid"'
        write(111,*) 'STRANDID=0, SOLUTIONTIME=0'
        write(111,*) 'I=258, J=258, K=1, ZONETYPE=Ordered'
        write(111,*) 'DATAPACKING=POINT'
        write(111,*) 'DT=(SINGLE SINGLE SINGLE SINGLE )'
	do j=1,258
	    do k=1,258
	        write(111,1111) ym1(j),zm1(k),uu1(j,k),upri(j,k)
1111	        format(5e17.9)
            enddo
        enddo
        close(111)
        
        open(121,file='..\output\upri.dat',form='formatted')
        write(121,*) 'TITLE    ="Plot3D DataSet"'
	write(121,*) 'VARIABLES = "y" "z" "upri"'
        write(121,*) 'I=258, J=258, K=1, ZONETYPE=Ordered'
        write(121,*) 'DATAPACKING=POINT'
        write(121,*) 'DT=(SINGLE SINGLE SINGLE )'
	do j=1,258
	    do k=1,258
	        write(121,1211) ym1(j),zm1(k),upri(j,k)
1211	        format(5e17.9)
            enddo
        enddo
        close(121)
        
        open(122,file='..\output\umod.dat',form='formatted')
        write(122,*) 'TITLE    ="Plot3D DataSet"'
	write(122,*) 'VARIABLES = "y" "z" "umod" '
        write(122,*) 'I=258, J=258, K=1, ZONETYPE=Ordered'
        write(122,*) 'DATAPACKING=POINT'
        write(122,*) 'DT=(SINGLE SINGLE SINGLE )'
	do j=1,258
	    do k=1,258
	        write(122,1212) ym1(j),zm1(k),uu1(j,k)
1212	        format(5e17.9)
            enddo
        enddo
        close(122)
        
        
        
        open(112,file='..\output\resume_initial.dat',form='formatted')
        do j=1,ncompt
            write(112,1112) resume_initial(j)
        enddo
1112	format(5e17.9)
        close(112)

	end
