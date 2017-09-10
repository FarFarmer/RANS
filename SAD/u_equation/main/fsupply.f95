!===============================!
        subroutine augrid
!===============================!
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!-------------------------------!
        ibg=1
        jbg=1
        kbg=1
!-------grid level one----------c
        call mgrdyz(jbg,kbg,m11,n11,nj1,nk1,ym1,yv1,dym1,dyv1,fy11,fy21,zm1,zw1,dzm1,dzw1,fz11,fz21)
        return
        end


!===============================!
        subroutine mgrdyz(j1,k1,m1,n1,nj,nk,ym,yv,dym,dyv,fy1,fy2,zm,zw,dzm,dzw,fz1,fz2)
!===============================!
        include 'table.prc'
        dimension ym(nj),yv(nj),fy1(nj),fy2(nj),dym(nj),dyv(nj),zm(nk),zw(nk),fz1(nk),fz2(nk),dzm(nk),dzw(nk)
!-------------------------------!
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        do j=j2,m2
            ym(j)=5.0d-01*(yv(j+1)+yv(j))
        enddo 
        ym(j1)=yv(j2)
        ym(m1)=yv(m1)
!        ym(j1)=yv(j2)-(ym(j2)-yv(j2))
!        ym(m1)=yv(m1)+(yv(m1)-ym(m2))
        do j=j2,m1
            dyv(j)=ym(j)-ym(j-1)
        enddo
        do j=j2,m2
            dym(j)=yv(j+1)-yv(j)
        enddo
        do j=j2,m1
            fy1(j)=(ym(j)-yv(j))/dyv(j)
            fy2(j)=(yv(j)-ym(j-1))/dyv(j)
        enddo

        do k=k2,n2
            zm(k)=5.0d-01*(zw(k+1)+zw(k))
        enddo
        zm(k1)=zw(k2)
        zm(n1)=zw(n1)
!        zm(k1)=zw(k2)-(zm(k2)-zw(k2))
!        zm(n1)=zw(n1)+(zw(n1)-zm(n2))
        do k=k2,n1
            dzw(k)=zm(k)-zm(k-1)
        enddo
        do k=k2,n2
            dzm(k)=zw(k+1)-zw(k)
        enddo
        do k=k2,n1
            fz1(k)=(zm(k)-zw(k))/dzw(k)
            fz2(k)=(zw(k)-zm(k-1))/dzw(k)
        enddo
        return
        end

!===============================!
        subroutine begins
!===============================!
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!-------------------------------!
        j11=1
        k11=1
        do j=j11,m11
            do k=k11,n11
                uu1(j,k)=ui1(j,k)
                j11=2;k11=1
            enddo
        enddo
        do j=j11,m11
            do k=k11,n11
                vv1(j,k)=vi1(j,k)
                j11=1;k11=2
            enddo
        enddo
        do j=j11,m11
            do k=k11,n11
                ww1(j,k)=wi1(j,k)
            enddo
        enddo
        return
        end
!===============================!
        subroutine update
!===============================!
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!-------------------------------!
        j11=1
        k11=1
        do j=j11,m11
            do k=k11,n11
                ui1(j,k)=uu1(j,k)
                j11=2;k11=1
            enddo
        enddo
        do j=j11,m11
            do k=k11,n11
                vi1(j,k)=vv1(j,k)
            enddo
        enddo
        j11=1
        k11=2
        do j=j11,m11
            do k=k11,n11
                wi1(j,k)=ww1(j,k)
            enddo
        enddo            
        return
        end
