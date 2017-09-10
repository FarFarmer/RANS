!c**********************c
        subroutine leinuo
!c**********************c
        include 'table.prc'
        include 'table.cre'
        include 'table.lei'
        
        open(20,file='..\input\vpvp.dat')
        read(20,*) vpvp
        close(20)

        open(21,file='..\input\vpwp.dat')
        read(21,*) vpwp
        close(21)
        
        open(22,file='..\input\wpwp.dat')
        read(22,*) wpwp
        close(22)
!        open(20,file='..\input\leinuo.dat')
!        do k=1,258
!	    do j=1,258
!	        read(20,1011) vpvp(j,k),vpwp(j,k),wpwp(j,k)
!1011	        format(3e17.9)
!            enddo
!        enddo
!        close(20)
        end
