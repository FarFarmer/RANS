
!===============================!
        subroutine initiu(kkk)
!===============================!
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.lei'     
!-------------------------------!
        j11=1;k11=1
        j21=j11+1;k21=k11+1
        m21=m11-1;n21=n11-1
!-------------------------------!
        do j=j11,m11
            do k=k11,n11
                aiw1(j,k)=0.0d+00;aie1(j,k)=0.0d+00
                ajs1(j,k)=0.0d+00;ajn1(j,k)=0.0d+00
                akb1(j,k)=0.0d+00;akt1(j,k)=0.0d+00
                app1(j,k)=0.0d+00;bpp1(j,k)=0.0d+00
            enddo
        enddo
!-------------------------------!
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) then
                    cycle
                endif
!-------doing velocity-------c
                uiw=0
                uie=0

                if(j.eq.j21.or.(j.eq.je1.and.(k.ge.kb1.and.k.lt.ke1))) then
                    uis=0.0d+00
                else
                    uis=ui1(j-1,k)*fy11(j)+ui1(j,k)*fy21(j)
                end if
                if(j.eq.m21.or.(j.eq.(jb1-1).and.(k.ge.kb1.and.k.lt.ke1))) then
                    uin=0.0d+00
                else
                    uin=ui1(j,k)*fy11(j+1)+ui1(j+1,k)*fy21(j+1)
                end if

                if(k.eq.k21.or.(k.eq.ke1.and.(j.ge.jb1.and.j.lt.je1))) then
                    uib=0.0d+00
                else
                    uib=ui1(j,k-1)*fz11(k)+ui1(j,k)*fz21(k)
                end if
                if(k.eq.n21.or.(k.eq.(kb1-1).and.(j.ge.jb1.and.j.lt.je1))) then
                    uit=0.0d+00
                else
                    uit=ui1(j,k)*fz11(k+1)+ui1(j,k+1)*fz21(k+1)
                end if

                viw=0;vie=0
                vin=vi1(j+1,k);vis=vi1(j,k)

                wiw=0;wie=0
                wit=wi1(j,k+1);wib=wi1(j,k)
!-------convective term-------c
                uux=0.0   
                uvy=(uin*vin-uis*vis)*dzm1(k)
                uwz=(uit*wit-uib*wib)*dym1(j)
!-----------------------------c
                volu=dym1(j)*dzm1(k)
                hu0(j,k)=-(uux+uvy+uwz)
                uupvp=(UPVP(j+1,k)-UPVP(j-1,k))/2./DYM1(j)
                uupwp=(UPWP(j,k+1)-UPWP(j,k-1))/2./DZM1(k)
                ulei=-(uupvp+uupwp)*volu
                bpp1(j,k)=(strh/dt*ui1(j,k)-prx)*volu+hu0(j,k)+ulei
!                bpp1(j,k)=(strh/dt*ui1(j,k)-prx)*volu+hu0(j,k)
!                bpp1(j,k)=-prx*volu+hu0(j,k)
            enddo
        enddo

        write(*,*)'*--------------------------------------------------*'
        write(*,*) '         kkk=',kkk
        write(*,*)'*--------------------------------------------------*'
        return
        end

