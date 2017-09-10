!**********************c
        subroutine leinuoU
!**********************c
        include 'table.prc'
        include 'table.cre'
        include 'table.lei'
        open(20,file='..\input\upvp.dat')
        read(20,*)upvp
        close(20)
        
        open(21,file='..\input\upwp.dat')
        read(21,*)upwp
        close(21)
        end
