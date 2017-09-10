!c===============================c
        subroutine equatp
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
        j11=1;k11=1
        j21=j11+1;k21=k11+1
        m21=m11-1;n21=n11-1
!c-------------------------------c
        divsum0=0.0d+00
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                aiw1(j,k)=0.0d+00
                aie1(j,k)=0.0d+00
                ajs1(j,k)=dzm1(k)/dyv1(j)
                ajn1(j,k)=dzm1(k)/dyv1(j+1)
                akb1(j,k)=dym1(j)/dzw1(k)
                akt1(j,k)=dym1(j)/dzw1(k+1)
                if(j.eq.m21.or.(j.eq.(jb1-1).and.(k.ge.kb1.and.k.lt.ke1))) ajn1(j,k)=0.0d+00
                if(k.eq.n21.or.(k.eq.(kb1-1).and.(j.ge.jb1.and.j.lt.je1))) akt1(j,k)=0.0d+00
                if(j.eq.j21.or.(j.eq.je1.and.(k.ge.kb1.and.k.lt.ke1))) ajs1(j,k)=0.0d+00
                if(k.eq.k21.or.(k.eq.ke1.and.(j.ge.jb1.and.j.lt.je1))) akb1(j,k)=0.0d+00
                app1(j,k)=aiw1(j,k)+aie1(j,k)+ajs1(j,k)+ajn1(j,k)+akb1(j,k)+akt1(j,k)
        
                dvy=(vv1(j+1,k)-vv1(j,k))*dzm1(k)
                dwz=(ww1(j,k+1)-ww1(j,k))*dym1(j)
                bpp1(j,k)=-(dvy+dwz)
                divsum0=divsum0+dabs(bpp1(j,k))
            enddo
        enddo
                
        app1(jr1,kr1)=1.0d+00
        aiw1(jr1,kr1)=0.0d+00
        aie1(jr1,kr1)=0.0d+00
        ajs1(jr1,kr1)=0.0d+00
        ajn1(jr1,kr1)=0.0d+00
        akb1(jr1,kr1)=0.0d+00
        akt1(jr1,kr1)=0.0d+00
        bpp1(jr1,kr1)=0.0d+00
        
        mode=2
!        mode=1
        if(levgrd.eq.1) then
            call solve1(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre,res1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.2) then
            call solve2(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre,res1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.3) then
            call solve3(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre,res1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.4) then
            call solve4(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre,res1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.5) then
            call solve5(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre,res1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.6) then
            call solve6(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre,res1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.7) then
            call solve7(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        pre,res1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
!c-------pressure boundary value-------c

        do j=j21,m21
            pre(j,k11)=pre(j,k21)+dzw1(k21)/dzw1(k21+1)*(pre(j,k21)-pre(j,k21+1))
            pre(j,n11)=pre(j,n21)+dzw1(n21+1)/dzw1(n21)*(pre(j,n21)-pre(j,n21-1))
        enddo

        do k=k21,n21
            pre(j11,k)=pre(j21,k)+dyv1(j21)/dyv1(j21+1)*(pre(j21,k)-pre(j21+1,k))
            pre(m11,k)=pre(m21,k)+dyv1(m21+1)/dyv1(m21)*(pre(m21,k)-pre(m21-1,k))
        enddo

        prej=pre(j21,k11)+dyv1(j21)/dyv1(j21+1)*(pre(j21,k11)-pre(j21+1,k11))
        prek=pre(j11,k21)+dzw1(k21)/dzw1(k21+1)*(pre(j11,k21)-pre(j11,k21+1))
        pre(j11,k11)=5.0d-01*(prej+prek)

        prej=pre(j21,n11)+dyv1(j21)/dyv1(j21+1)*(pre(j21,n11)-pre(j21+1,n11))
        prek=pre(j11,n21)+dzw1(n21+1)/dzw1(n21)*(pre(j11,n21)-pre(j11,n21-1))
        pre(j11,n11)=5.0d-01*(prej+prek)

        prej=pre(m21,k11)+dyv1(m21+1)/dyv1(m21)*(pre(m21,k11)-pre(m21-1,k11))
        prek=pre(m11,k21)+dzw1(k21)/dzw1(k21+1)*(pre(m11,k21)-pre(m11,k21+1))
        pre(m11,k11)=5.0d-01*(prej+prek)

        prej=pre(m21,n11)+dyv1(m21+1)/dyv1(m21)*(pre(m21,n11)-pre(m21-1,n11))
        prek=pre(m11,n21)+dzw1(n21+1)/dzw1(n21)*(pre(m11,n21)-pre(m11,n21-1))
        pre(m11,n11)=5.0d-01*(prej+prek)

        do j=jb1+1,je1-2
            pre(j,ke1-1)=pre(j,ke1)+dzw1(ke1)/dzw1(ke1+1)*(pre(j,ke1)-pre(j,ke1+1))
            pre(j,kb1)=pre(j,kb1-1)+dzw1(kb1)/dzw1(kb1-1)*(pre(j,kb1-1)-pre(j,kb1-2))
        enddo

        do k=kb1+1,ke1-2
            pre(je1-1,k)=pre(je1,k)+dyv1(je1)/dyv1(je1+1)*(pre(je1,k)-pre(je1+1,k))
            pre(jb1,k)=pre(jb1-1,k)+dyv1(jb1)/dyv1(jb1-1)*(pre(jb1-1,k)-pre(jb1-2,k))
        enddo

        prej=pre(je1,ke1-1)+dyv1(je1)/dyv1(je1+1)*(pre(je1,ke1-1)-pre(je1+1,ke1-1))
        prek=pre(je1-1,ke1)+dzw1(ke1)/dzw1(ke1+1)*(pre(je1-1,ke1)-pre(je1-1,ke1+1))
        pre(je1-1,ke1-1)=5.0d-01*(prej+prek)

        prej=pre(je1,kb1)+dyv1(je1)/dyv1(je1+1)*(pre(je1,kb1)-pre(je1+1,kb1))
        prek=pre(je1,kb1-1)+dzw1(kb1)/dzw1(kb1-1)*(pre(je1,kb1-1)-pre(je1,kb1-2))
        pre(je1-1,kb1)=5.0d-01*(prej+prek)

        prej=pre(jb1-1,ke1)+dyv1(jb1)/dyv1(jb1-1)*(pre(jb1-1,ke1)-pre(jb1-2,ke1))
        prek=pre(jb1,ke1)+dzw1(ke1)/dzw1(ke1+1)*(pre(jb1,ke1)-pre(jb1,ke1+1))
        pre(jb1,ke1-1)=5.0d-01*(prej+prek)

        prej=pre(jb1-1,kb1)+dyv1(jb1)/dyv1(jb1-1)*(pre(jb1-1,kb1)-pre(jb1-2,kb1))
        prek=pre(jb1,kb1-1)+dzw1(kb1)/dzw1(kb1-1)*(pre(jb1,kb1-1)-pre(jb1,kb1-2))
        pre(jb1,kb1)=5.0d-01*(prej+prek)


!c-------correction of v by p-------c
        j11=2
        k11=1
        j21=j11+1
        k21=k11+1
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.le.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle

                vv1(j,k)=vv1(j,k)-(pre(j,k)-pre(j-1,k))/dyv1(j)
            enddo
        enddo
!c-------correction of w by p-------c
        j11=1
        k11=2
        j21=j11+1
        k21=k11+1
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.le.ke1)) cycle
                ww1(j,k)=ww1(j,k)-(pre(j,k)-pre(j,k-1))/dzw1(k)
            enddo
        enddo

!c-------continuity satisfaction level check-------c
        divsum=0.0d+00
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                dvy=(vv1(j+1,k)-vv1(j,k))*dzm1(k)
                dwz=(ww1(j,k+1)-ww1(j,k))*dym1(j)
                divsum=divsum+dabs(dvy+dwz)
            enddo
        enddo
        
        write(*,601) divsum0,divsum
601     format(' *',8x,'divergence drop',1pe9.2,'/',1pe9.2,8x,'*')
        return
        end
