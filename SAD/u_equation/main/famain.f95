!===============================!
        include 'table.mdu'
!===============================!
        program lesannsda
!-------------------------------!

        call ipd0fm
        call ipd1uf
        call ipdff0u
        call augrid
	
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        write(*,*)'*    The program is to calculate 2D-rans  u        *'
        write(*,*)'*                                                  *'
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        call begins
        call nseqts
        write(*,*)'*         The----End----Of---The---Program         *'
        write(*,*)'*                                                  *'
        write(*,*)'*--------------------------------------------------*'
        stop
        end

!======================================================================c
        subroutine nseqts 
!======================================================================c
        use varalc
        include 'table.prc'
!----------------------------------------------------------------------c        
        call cpu_time(time0)
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        write(*,*)'*        solving  process for momentum Eqs.        *'
        write(*,*)'*                                                  *'
        call leinuoU
        
        ncompt=10000
        allocate(resume_initial(ncompt))
        do kkk=1,ncompt
            write(*,*)'*        solving  momentum  equation  of  u        *'
            call initiu(kkk)
            call equatu(kkk)
            write(*,*)'*                                                  *'
            write(*,*)'*--------end of handling momentum equations--------*'
            call update
        enddo

        call shuchu2(ncompt) 
        !----------------------------------------------------------------------c
        call cpu_time(time1)
        write(*,*)'         total computation time',time1-time0
        write(*,*)'*                                                  *'
        write(*,*)'*--------------------------------------------------*'
        return
        end

