!c===============================c
        subroutine equatu
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.lei'
        dimension wmv1212(258,258)
!c-------------------------------c
        open(19,file='..\input\vav.dat')
        read(19,*)vv1
        close(19)
        open(19,file='..\input\wav.dat')
        read(19,*)ww1
        close(19)
        open(19,file='..\input\mv1212.dat')
        read(19,*)wmv1212
        close(19)
        
        j11=1;k11=1
        j21=j11+1;k21=k11+1
        m21=m11-1;n21=n11-1
        coef=1./reno
        do k=k11,n11
            do j=j11,m11
                bpp1(j,k)=2*dym1(j)*dzm1(k)
                uu1(j,k)=0
            enddo
        enddo
          
!c-------doing ajs,ajn-------c
        do k=k21,n21
            do j=j21,m21
                coef=(1./reno-wmv1212(j,k))
                fact=coef*dzm1(k)
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                ajs1(j,k)=fact/dyv1(j)
                ajn1(j,k)=fact/dyv1(j+1)
            enddo
        enddo
!c-------doing akb,akt-------c
        do j=j21,m21
            do k=k21,n21
                coef=(1./reno-wmv1212(j,k))
                fact=coef*dym1(j)
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                akb1(j,k)=fact/dzw1(k)
                akt1(j,k)=fact/dzw1(k+1)
            enddo
        enddo
!c-------------------------------c
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                app1(j,k)=ajn1(j,k)+ajs1(j,k)+akt1(j,k)+akb1(j,k)+0.5*(vv1(j+1,k)-vv1(j,k))*dzm1(k) &
                          +0.5*(ww1(j,k+1)-ww1(j,k))*dym1(j)
                ajs1(j,k)=ajs1(j,k)+0.5*vv1(j,k)*dzm1(k)
                ajn1(j,k)=ajn1(j,k)-0.5*vv1(j+1,k)*dzm1(k)
                akb1(j,k)=akb1(j,k)+0.5*ww1(j,k)*dym1(k)
                akt1(j,k)=akt1(j,k)-0.5*ww1(j,k+1)*dym1(k)
                aiw1(j,k)=0
                aie1(j,k)=0
            enddo
        enddo
        call bondu2(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,uu1,app1,bpp1,aiw1,ajs1,akb1,aie1,ajn1,akt1)
!c********************************c
        call residum(j21,k21,m21,n21,jb1,je1,kb1,ke1,ni1,nj1,nk1,ressum0,uu1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1,res1)
        ressum1=ressum0
        write(*,200) ressum0
        iter=0
        mode=2
        
!c--------------------------------c
100     iter=iter+1

!c--------------------------------c
        do k=k21,n21
            if(k.ge.kb1.and.k.lt.ke1) then
                call trdgmj(j11,jb1,k,nj1,nk1,uu1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
                call trdgmj(je1-1,m11,k,nj1,nk1,uu1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
            else
                call trdgmj(j11,m11,k,nj1,nk1,uu1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
            end if
        enddo
!c--------------------------------c
        do j=j21,m21
            if(j.ge.jb1.and.j.lt.je1) then
                call trdgmk(k11,kb1,j,nj1,nk1,uu1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
                call trdgmk(ke1-1,n11,j,nj1,nk1,uu1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
            else
                call trdgmk(k11,n11,j,nj1,nk1,uu1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
            end if
        enddo

!c--------------------------------c
154     call residum(j21,k21,m21,n21,jb1,je1,kb1,ke1,ni1,nj1,nk1,ressum,uu1,app1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1,res1)
        write(*,201) ressum

        if(ressum.ge.ressum1) then
            write(*,202)
!c        goto 1000
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
    
1000      open(91,file='..\output\uav.dat',form='formatted')
	  do 157 j=1,258
	  do 157 k=1,258
	  write(91,1010)ym1(j),zm1(k), uu1(j,k),vv1(j,k),ww1(j,k),wmv1212(j,k)
1010	  format(6e17.9)
157	  continue
        return
        end
