
!===============================!
        subroutine residum(j2,k2,m2,n2,jb,je,kb,ke,ni,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!===============================!
        include 'table.prc'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!*******************************!
        ressum=0.0d+00
        do 101 j=j2,m2
        do 101 k=k2,n2
        if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) goto 101
        res(j,k)=bpp(j,k)-app(j,k)*ff(j,k)+aiw(j,k)*ff(j,k)+aie(j,k)*ff(j,k)+ajs(j,k)*ff(j-1,k)+ajn(j,k)*ff(j+1,k)+ &
                 akb(j,k)*ff(j,k-1)+akt(j,k)*ff(j,k+1)
        ressum=ressum+dabs(res(j,k))
101     continue
        return
        end

!********************************c
        subroutine trdgmj(j1,m1,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!********************************c
        include 'table.prc'
        include 'table.dmw'
        dimension p(njw),q(njw)
        dimension ff(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
!--------------------------------c
        j2=j1+1
        j3=j2+1
        m2=m1-1
        m3=m2-1
        denom=1.0d+00/app(j2,k)
        p(j2)=ajn(j2,k)*denom
        temp=bpp(j2,k)+akb(j2,k)*ff(j2,k-1)+akt(j2,k)*ff(j2,k+1)
        q(j2)=temp*denom
        do j=j3,m2
            denom=1.0d+00/(app(j,k)-p(j-1)*ajs(j,k))
            p(j)=ajn(j,k)*denom
            temp=bpp(j,k)+akb(j,k)*ff(j,k-1)+akt(j,k)*ff(j,k+1)
            q(j)=(temp+ajs(j,k)*q(j-1))*denom
        enddo
        ff(m2,k)=q(m2)
        do j=m3,j2,-1
            ff(j,k)=ff(j+1,k)*p(j)+q(j)
        enddo
        return
        end

!********************************c
        subroutine trdgmk(k1,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!********************************c
        include 'table.prc'               
        include 'table.dmw'
        dimension p(nkw),q(nkw)
        dimension ff(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
!--------------------------------c
        k2=k1+1
        k3=k2+1
        n2=n1-1
        n3=n2-1
        denom=1.0d+00/app(j,k2)
        p(k2)=akt(j,k2)*denom
        temp=bpp(j,k2)+ajs(j,k2)*ff(j-1,k2)+ajn(j,k2)*ff(j+1,k2)
        q(k2)=temp*denom
        do k=k3,n2
            denom=1.0d+00/(app(j,k)-p(k-1)*akb(j,k))
            p(k)=akt(j,k)*denom
            temp=bpp(j,k)+ajs(j,k)*ff(j-1,k)+ajn(j,k)*ff(j+1,k)
            q(k)=(temp+akb(j,k)*q(k-1))*denom
        enddo
        ff(j,n2)=q(n2)
        do  k=n3,k2,-1
             ff(j,k)=ff(j,k+1)*p(k)+q(k)
        enddo
        return
        end


!********************************c
        subroutine trdgpj(j1,m1,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!********************************c
        include 'table.prc'
        include 'table.dmw'
        dimension p(njw),q(njw)
        dimension ff(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
!--------------------------------c 
        j2=j1+1
        j3=j2+1
        m2=m1-1
        m3=m2-1


        denom=1.0d+00/app(j2,k)
        p(j2)=ajn(j2,k)*denom
        temp=bpp(j2,k)+akb(j2,k)*ff(j2,k-1)+akt(j2,k)*ff(j2,k+1)
        q(j2)=temp*denom
        do 100 j=j3,m2
        denom=1.0d+00/(app(j,k)-p(j-1)*ajs(j,k))
        p(j)=ajn(j,k)*denom
        temp=bpp(j,k)+akb(j,k)*ff(j,k-1)+akt(j,k)*ff(j,k+1)
        q(j)=(temp+ajs(j,k)*q(j-1))*denom
100     continue

102     ff(m2,k)=q(m2)
        do 101 j=m3,j2,-1
101     ff(j,k)=ff(j+1,k)*p(j)+q(j)
        return
        end

!********************************c
        subroutine trdgpk(k1,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!********************************c
        include 'table.prc'
        include 'table.dmw'
        dimension p(nkw),q(nkw)
        dimension ff(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
!--------------------------------c
        k2=k1+1
        k3=k2+1
        n2=n1-1
        n3=n2-1


        denom=1.0d+00/app(j,k2)
        p(k2)=akt(j,k2)*denom
        temp=bpp(j,k2)+ajs(j,k2)*ff(j-1,k2)+ajn(j,k2)*ff(j+1,k2)
        q(k2)=temp*denom
        do 100 k=k3,n2
        denom=1.0d+00/(app(j,k)-p(k-1)*akb(j,k))
        p(k)=akt(j,k)*denom
        temp=bpp(j,k)+ajs(j,k)*ff(j-1,k)+ajn(j,k)*ff(j+1,k)
        q(k)=(temp+akb(j,k)*q(k-1))*denom
100     continue

102     ff(j,n2)=q(n2)
        do 101 k=n3,k2,-1
101     ff(j,k)=ff(j,k+1)*p(k)+q(k)
        return
        end
