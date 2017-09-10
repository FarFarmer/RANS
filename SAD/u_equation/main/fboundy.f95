!*******************************!
        subroutine bondu2(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,uu1,app,con,aim,ajm,akm,aip,ajp,akp)
!*******************************!
        include 'table.prc'
        dimension app(nj,nk),con(nj,nk),aim(nj,nk),aip(nj,nk),ajm(nj,nk),ajp(nj,nk),akm(nj,nk),akp(nj,nk)
        dimension uu1(nj,nk)
!*******************************!
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
!-------------------------------!
        do 101 k=k2,n2
        j=j2
        con(j,k)=con(j,k)+ajm(j,k)*uu1(j-1,k)
        ajm(j,k)=0.0d+00
        j=m2
        con(j,k)=con(j,k)+ajp(j,k)*uu1(j+1,k)
        ajp(j,k)=0.0d+00
101     continue
!-------------------------------!
        do 102 j=j2,m2
        k=k2
        con(j,k)=con(j,k)+akm(j,k)*uu1(j,k-1)
        akm(j,k)=0.0d+00
        k=n2
        con(j,k)=con(j,k)+akp(j,k)*uu1(j,k+1)
        akp(j,k)=0.0d+00
102     continue
!-------------------------------! 
        do 103 k=kb,ke-1
        j=je
        con(j,k)=con(j,k)+ajm(j,k)*uu1(j-1,k)
        ajm(j,k)=0.0d+00
        j=jb-1
        con(j,k)=con(j,k)+ajp(j,k)*uu1(j+1,k)
        ajp(j,k)=0.0d+00
103     continue
!-------------------------------!
        do 104 j=jb,je-1
        k=ke
        con(j,k)=con(j,k)+akm(j,k)*uu1(j,k-1)
        akm(j,k)=0.0d+00
        k=kb-1
        con(j,k)=con(j,k)+akp(j,k)*uu1(j,k+1)
        akp(j,k)=0.0d+00
104     continue
!-------------------------------!
        return
        end
