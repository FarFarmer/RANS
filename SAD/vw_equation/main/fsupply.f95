!c===============================c
        subroutine augrid
!c===============================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
        ibg=1
        jbg=1
        kbg=1
!c-------grid level one----------c
        call mgrdyz(jbg,kbg,m11,n11,nj1,nk1,ym1,yv1,dym1,dyv1,fy11,fy21,zm1,zw1,dzm1,dzw1,fz11,fz21)
        return
        end


!c===============================c
        subroutine mgrdyz(j1,k1,m1,n1,nj,nk,ym,yv,dym,dyv,fy1,fy2,zm,zw,dzm,dzw,fz1,fz2)
!c===============================c
        include 'table.prc'
        dimension ym(nj),yv(nj),fy1(nj),fy2(nj),dym(nj),dyv(nj),zm(nk),zw(nk),fz1(nk),fz2(nk),dzm(nk),dzw(nk)
!c-------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        do 21 j=j2,m2
21      ym(j)=5.0d-01*(yv(j+1)+yv(j))
!        ym(j1)=yv(j2)
!        ym(m1)=yv(m1)
        ym(j1)=yv(j2)-(ym(j2)-yv(j2))
        ym(m1)=yv(m1)+(yv(m1)-ym(m2))
        do 22 j=j2,m1
22      dyv(j)=ym(j)-ym(j-1)
        do 23 j=j2,m2
23      dym(j)=yv(j+1)-yv(j)
        do 24 j=j2,m1
        fy1(j)=(ym(j)-yv(j))/dyv(j)
24      fy2(j)=(yv(j)-ym(j-1))/dyv(j)

        do 31 k=k2,n2
31      zm(k)=5.0d-01*(zw(k+1)+zw(k))
!        zm(k1)=zw(k2)
!        zm(n1)=zw(n1)
        zm(k1)=zw(k2)-(zm(k2)-zw(k2))
        zm(n1)=zw(n1)+(zw(n1)-zm(n2))
        do 32 k=k2,n1
32      dzw(k)=zm(k)-zm(k-1)
        do 33 k=k2,n2
33      dzm(k)=zw(k+1)-zw(k)
        do 34 k=k2,n1
        fz1(k)=(zm(k)-zw(k))/dzw(k)
34      fz2(k)=(zw(k)-zm(k-1))/dzw(k)
        return
        end

!c===============================c
        subroutine begins
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
        j11=2
        k11=1
        do 112 j=j11,m11
        do 112 k=k11,n11
112     vv1(j,k)=vi1(j,k)
        j11=1
        k11=2
        do 113 j=j11,m11
        do 113 k=k11,n11
113     ww1(j,k)=wi1(j,k)
        return
        end
!c===============================c
        subroutine update
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c

        i11=1
        j11=2
        k11=1
        do 112 j=j11,m11
        do 112 k=k11,n11
112     vi1(j,k)=vv1(j,k)
        j11=1
        k11=2
        do 113 j=j11,m11
        do 113 k=k11,n11
113     wi1(j,k)=ww1(j,k)
        return
        end
