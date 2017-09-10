!c===============================c
        subroutine initiv
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.lei'
!c-------------------------------c
        m21=m11-1
        n21=n11-1
!c-------------------------------c
        j11=1
        k11=1
        do j=j11,m11
            do k=k11,n11
                aiw1(j,k)=0.0d+00
                aie1(j,k)=0.0d+00
                ajs1(j,k)=0.0d+00
                ajn1(j,k)=0.0d+00
                akb1(j,k)=0.0d+00
                akt1(j,k)=0.0d+00
                app1(j,k)=0.0d+00
                bpp1(j,k)=0.0d+00
            enddo
        enddo
!c-------------------------------c
        j11=2
        k11=1
        j21=j11+1
        k21=k11+1
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.le.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
!c-------doing velocity-------c
                viw=0
	        vie=0
                vis=5.0d-01*(vi1(j-1,k)+vi1(j,k))
                vin=5.0d-01*(vi1(j+1,k)+vi1(j,k))
                vib=(vi1(j,k-1)*fz11(k)+vi1(j,k)*fz21(k))
                vit=(vi1(j,k)*fz11(k+1)+vi1(j,k+1)*fz21(k+1))

                wit=wi1(j-1,k+1)*fy11(j)+wi1(j,k+1)*fy21(j)
                wib=wi1(j-1,k)*fy11(j)+wi1(j,k)*fy21(j)
                wis=5.0d-01*(wi1(j-1,k+1)+wi1(j-1,k))
                win=5.0d-01*(wi1(j,k+1)+wi1(j,k))
!c-------convective term-------c
                vvy=(vin*vin-vis*vis)*dzm1(k)
                vwz=(vit*wit-vib*wib)*dyv1(j)
	        vux=0
!c-----------------------------c
                volv=dyv1(j)*dzm1(k)
                hv0(j,k)=-(vvy+vwz+vux)
                
!                vvpvp=(VPVP(j,k)-VPVP(j-1,k))/DYV1(j)
!                vvpwp=(0.25*(VPWP(j-1,k)+VPWP(j-1,k+1)+VPWP(j,k)+VPWP(j,k+1)) &
!                       -0.25*(VPWP(j-1,k-1)+VPWP(j-1,k)+VPWP(j,k-1)+VPWP(j,k)))/DZM1(k)
!                vlei=-(vvpvp+vvpwp)*volv

                vvpvp=(vpvp(j,k)-vpvp(j-1,k))*dzm1(k)
                
                if (k==n21.or.k==kb1) then
                    vpwpt=0.0d+00
                else
                    vpwpt=(vpwp(j-1,k)*fy11(j)+vpwp(j,k)*fy21(j))*fz11(k+1)+(vpwp(j-1,k+1)*fy11(j)+vpwp(j,k+1)*fy21(j))*fz21(k+1)
                endif
                
                if (k==k21.or.k==ke1) then
                    vpwpb=0.0d+00
                else
                    vpwpb=(vpwp(j-1,k-1)*fy11(j)+vpwp(j,k-1)*fy21(j))*fz11(k)+(vpwp(j-1,k)*fy11(j)+vpwp(j,k)*fy21(j))*fz21(k)
                endif
                
                vvpwp=(vpwpt-vpwpb)*dyv1(j)
                vlei=-(vvpvp+vvpwp)
                
                bpp1(j,k)=(strh/dt*vi1(j,k))*volv+hv0(j,k)+vlei
!                bpp1(j,k)=(strh/dt*vi1(j,k)-pry)*volv+hv0(j,k)+vlei
!                bpp1(j,k)=(-pry)*volv+hv0(j,k)+vlei

            enddo
        enddo
        return
        end

!c===============================c
        subroutine initiw
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.lei'
!c-------------------------------c
        m21=m11-1
        n21=n11-1
!c-------------------------------c
        j11=1
        k11=1
        do j=j11,m11
            do k=k11,n11
                aiw1(j,k)=0.0d+00
                aie1(j,k)=0.0d+00
                ajs1(j,k)=0.0d+00
                ajn1(j,k)=0.0d+00
                akb1(j,k)=0.0d+00
                akt1(j,k)=0.0d+00
                app1(j,k)=0.0d+00
                bpp1(j,k)=0.0d+00
            enddo
        enddo
!c-------------------------------c
        j11=1
        k11=2
        j21=j11+1
        k21=k11+1
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.le.ke1)) cycle
!c-------doing velocity-------c
                wiw=0
	        wie=0
                wis=(wi1(j-1,k)*fy11(j)+wi1(j,k)*fy21(j))
                win=(wi1(j,k)*fy11(j+1)+wi1(j+1,k)*fy21(j+1))
                wib=5.0d-01*(wi1(j,k-1)+wi1(j,k))
                wit=5.0d-01*(wi1(j,k+1)+wi1(j,k))

                vin=vi1(j+1,k-1)*fz11(k)+vi1(j+1,k)*fz21(k)
                vis=vi1(j,k-1)*fz11(k)+vi1(j,k)*fz21(k)
                vib=5.0d-01*(vi1(j,k-1)+vi1(j+1,k-1))
                vit=5.0d-01*(vi1(j,k)+vi1(j+1,k))
!c-------convective term-------c
                wwz=(wit*wit-wib*wib)*dym1(j)
	        wux=0
                wvy=(win*vin-wis*vis)*dzw1(k)
!c-----------------------------c
                volw=dzw1(k)*dym1(j)
                hw0(j,k)=-(wwz+wux+wvy)
                
!                wvpwp=(0.25*(VPWP(j,k)+VPWP(j,k-1)+VPWP(j+1,k)+VPWP(j+1,k-1)) &
!                       -0.25*(VPWP(j-1,k)+VPWP(j-1,k-1)+VPWP(j,k)+VPWP(j,k-1)))/DYM1(j)
!                wwpwp=(WPWP(j,k)-WPWP(j,k-1))/DZW1(k)
!                wlei=-(wvpwp+wwpwp)*volw
                
                wwpwp=(wpwp(j,k)-wpwp(j,k-1))*dym1(j)
                
                if (j==m21.or.j==jb1) then
                    vpwpn=0.0d+00
                else
                    vpwpn=(vpwp(j,k-1)*fz11(k)+vpwp(j,k)*fz21(k))*fy11(j+1)+(vpwp(j+1,k-1)*fz11(k)+vpwp(j+1,k)*fz21(k))*fy21(j+1)
                endif
                
                if (j==j21.or.j==je1) then
                    vpwps=0.0d+00
                else
                    vpwps=(vpwp(j-1,k-1)*fz11(k)+vpwp(j-1,k)*fz21(k))*fy11(j)+(vpwp(j,k-1)*fz11(k)+vpwp(j,k)*fz21(k))*fy21(j)
                endif

                wvpwp=(vpwpn-vpwps)*dzw1(k)
                wlei=-(wwpwp+wvpwp)

                bpp1(j,k)=(strh/dt*wi1(j,k))*volw+hw0(j,k)+wlei
!                bpp1(j,k)=(strh/dt*wi1(j,k)-prz)*volw+hw0(j,k)+wlei
!                bpp1(j,k)=(-prz)*volw+hw0(j,k)+wlei
            enddo
        enddo
        return
        end

!c===============================c
        subroutine initip
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
        j11=1
        k11=1
        do j=j11,m11
            do k=k11,n11
                pre(j,k)=0.0d+00
                bpp1(j,k)=0.0d+00
            enddo
        enddo
        do j=j11,m11
            do k=k11,n11
                aiw1(j,k)=0.0d+00
                aie1(j,k)=0.0d+00
                ajs1(j,k)=0.0d+00
                ajn1(j,k)=0.0d+00
                akb1(j,k)=0.0d+00
                akt1(j,k)=0.0d+00
                app1(j,k)=0.0d+00
            enddo
        enddo
        return
        end
