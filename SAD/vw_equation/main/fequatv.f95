!c===============================c
        subroutine equatv
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
        j11=2;k11=1
        j21=j11+1;k21=k11+1
        m21=m11-1;n21=n11-1
        coef=5.0d-01/reno
!c-------------------------------c
!c-------doing vyy,ajs,ajn-------c
        do k=k21,n21
            fact=coef*dzm1(k)
            do j=j21,m21
                if((j.ge.jb1.and.j.le.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                ajn1(j,k)=fact/dym1(j)
                ajs1(j,k)=fact/dym1(j-1)
                viyy=fact*((vi1(j+1,k)-vi1(j,k))/dym1(j)-(vi1(j,k)-vi1(j-1,k))/dym1(j-1))
                bpp1(j,k)=bpp1(j,k)+viyy
            enddo
        enddo
!c-------doing vzz,akb,akt-------c
        do j=j21,m21
            fact=coef*dyv1(j)
            do k=k21,n21
                if((j.ge.jb1.and.j.le.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                akb1(j,k)=fact/dzw1(k)
                akt1(j,k)=fact/dzw1(k+1)
                vizz=fact*((vi1(j,k+1)-vi1(j,k))/dzw1(k+1)-(vi1(j,k)-vi1(j,k-1))/dzw1(k))
                bpp1(j,k)=bpp1(j,k)+vizz
            enddo
        enddo
!c-------doing vxx,aiw,aie-------c
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.le.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                aiw1(j,k)=0
                aie1(j,k)=0
            enddo
        enddo
!c-------------------------------c
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.le.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                app1(j,k)=strh/dt*dyv1(j)*dzm1(k)+aie1(j,k)+aiw1(j,k)+ajn1(j,k)+ajs1(j,k)+akt1(j,k)+akb1(j,k)
!                app1(j,k)=aie1(j,k)+aiw1(j,k)+ajn1(j,k)+ajs1(j,k)+akt1(j,k)+akb1(j,k)

            enddo
        enddo
!c********************************c
        call bondv2(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,vv1,app1,bpp1,aiw1,ajs1,akb1,aie1,ajn1,akt1)
!c********************************c
        call residum(j21,k21,m21,n21,jb1,je1+1,kb1,ke1,ni1,nj1,nk1,ressum0,vv1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1,res1)
        ressum1=ressum0
        write(*,200) ressum0
        iter=0
        mode=2
!c--------------------------------c
100     iter=iter+1
!c--------------------------------c

!c--------------------------------c
        do k=k21,n21
            if(k.ge.kb1.and.k.lt.ke1) then
                call trdgmj(j11,jb1,k,nj1,nk1,vv1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
                call trdgmj(je1,m11,k,nj1,nk1,vv1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
            else
                call trdgmj(j11,m11,k,nj1,nk1,vv1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
            end if
        enddo
!c--------------------------------c
        do j=j21,m21
            if(j.ge.jb1.and.j.le.je1) then
                call trdgmk(k11,kb1,j,nj1,nk1,vv1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
                call trdgmk(ke1-1,n11,j,nj1,nk1,vv1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
            else
                call trdgmk(k11,n11,j,nj1,nk1,vv1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
            end if
        enddo
!c--------------------------------c
154     call residum(j21,k21,m21,n21,jb1,je1+1,kb1,ke1,ni1,nj1,nk1,ressum,vv1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1,res1)
        write(*,201) ressum

        if(ressum.ge.ressum1) then
            write(*,202)
            goto 1000
        end if

!c mode=1: control residual precision
        if(mode.eq.1) then
            if(ressum.le.epsm21) then
                write(*,203)
                goto 1000
            else
                ressum1=ressum
                write(*,204)
                goto 100
            end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
            if(ressum.le.epsm21) then
                write(*,203)
                goto 1000
            end if
            declev=ressum/ressum0
            if(declev.gt.epsm31) then
                ressum1=ressum
                write(*,205)
                goto 100
            else
                write(*,206)
                goto 1000
            end if
        end if

!c mode=3: control residual decreasing rate
        if(mode.eq.3) then
            if(ressum.le.epsm21) then
                write(*,203)
                goto 1000
            end if
            decrat=ressum/ressum1
            if(decrat.le.epsm41) then
                ressum1=ressum
                write(*,207)
                goto 100
            else
                write(*,208)
                goto 1000
            end if
        end if

200     format(' *',8x,'initial residual ressum=',1pe10.3,8x,'*')
201     format(' *',8x,'total  residual  ressum=',1pe10.3,8x,'*')
202     format(' *        residual not improving,     return        *')
203     format(' *        total residual satisfied,   return        *')
204     format(' *        residual not satisfied, continuing        *')
205     format(' *        declev not satisfied,   continuing        *')
206     format(' *        declev satisfied,  go out of L-B-L        *')
207     format(' *        decrat small, go on next iteration        *')
208     format(' *        decrat large, stop the LBL process        *')

1000    return
        end
