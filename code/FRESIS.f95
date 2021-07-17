MODULE MODULE1
IMPLICIT NONE
!----------------INTEGER
INTEGER,PARAMETER :: IM=1025
INTEGER,PARAMETER :: JM=257
INTEGER,PARAMETER :: KMAX=2
INTEGER,PARAMETER :: NPRINT1=2
INTEGER,PARAMETER :: NPRINT2=10
INTEGER,PARAMETER :: ITERM=1000000
INTEGER,PARAMETER :: KSTART=0
!----------------REAL
REAL(8),PARAMETER :: PI=3.1415926535897932384626433832795D0
REAL(8),PARAMETER :: D=1.2D0
REAL(8),PARAMETER :: RC=0.5D0*D
REAL(8),PARAMETER :: THICKNESS=0.2D0
REAL(8),PARAMETER :: RO=0.5
REAL(8),PARAMETER :: RI=RO-THICKNESS
REAL(8),PARAMETER :: H=2.D0
REAL(8),PARAMETER :: L=4.D0*H
REAL(8),PARAMETER :: LY=H
REAL(8),PARAMETER :: LX=L
REAL(8),PARAMETER :: XC=-0.25D0*L
REAL(8),PARAMETER :: YC=0.D0
REAL(8),PARAMETER :: DX=L/(IM-1)
REAL(8),PARAMETER :: DY=H/(JM-1)
REAL(8),PARAMETER :: DX2=DX*DX
REAL(8),PARAMETER :: DY2=DY*DY
REAL(8),PARAMETER :: FDXDY=4.D0*DX*DY
REAL(8),PARAMETER :: BETA=DX/DY
REAL(8),PARAMETER :: BETA2=BETA*BETA
REAL(8),PARAMETER :: DT=5.D-6
REAL(8),PARAMETER :: UMAX=10.D0
REAL(8),PARAMETER :: NUF=0.1D0
REAL(8),PARAMETER :: NUS=0.1D0
REAL(8),PARAMETER :: RHO=1.D0
REAL(8),PARAMETER :: INV_RHO=1.D0/RHO
REAL(8),PARAMETER :: C1=2000.D0
REAL(8),PARAMETER :: C2=0.D0
REAL(8),PARAMETER :: C3=0.D0
REAL(8),PARAMETER :: EPS=1.D-6
REAL(8),PARAMETER :: CX=DT/DX
REAL(8),PARAMETER :: CY=DT/DY
REAL(8),PARAMETER :: CX2=CX*CX
REAL(8),PARAMETER :: CY2=CY*CY
REAL(8),PARAMETER :: PHIS_MIN=0.1D0
!------------------------------------------------------------------------------
LOGICAL,PARAMETER :: ICFLAG=.FALSE.                !TRUE if exists the initial file as initial condition or else FALSE
LOGICAL,PARAMETER :: BCFLAG=.TRUE.                !TRUE for CBC and FALSE for NBC
LOGICAL,PARAMETER :: WRITE_POINT_STATUS=.FALSE.   !TRUE for writing the point variables or else FALSE
LOGICAL,PARAMETER :: WRITE_FIELD_STATUS=.TRUE.    !TRUE for writing the field variables or else FALSE
!------------------------------------------------------------------------------
END MODULE
!==============================================================================
MODULE MULTIGRID
IMPLICIT NONE
!------------------INTEGER
INTEGER,PARAMETER            :: NG_X=10                                 !NUMBER OF COARSE GRID IN X DIRECTION
INTEGER,PARAMETER            :: NG_Y=8                                       !NUMBER OF COARSE GRID IN Y DIRECTION
INTEGER,PARAMETER            :: NG_M=(NG_X+NG_Y)/2+2                             !AVERAGE OF NG_X,NG_Y
INTEGER,PARAMETER            :: NG_D=ABS(NG_X-NG_Y)                              !DIFFRENCE OF NG_X,NG_Y
INTEGER,PARAMETER            :: NG_L=MIN(NG_X,NG_Y)                              !MINIMOM OF NG_X,NG_Y
INTEGER,PARAMETER            :: NCYCLE=2                                         !NUMBER OF MG V-CYCLE
INTEGER,PARAMETER            :: NPRE=1                                           !NUMBER OF PRESMOOTHING LOOP
INTEGER,PARAMETER            :: NPOST=1                                          !NUMBER OF POST SMOOTHING LOOP
INTEGER,PARAMETER            :: NFINAL=1                                         !NUMBER OF CORSEST GRID
INTEGER,PARAMETER            :: MEMLEN=13*2**(2*NG_M)/3+14*2**NG_M+8*NG_M-100/3  !MEMORY LENGHT FOR ARRAY
!------------------REAL
!------------------CHAR
END MODULE MULTIGRID
!==============================================================================
PROGRAM MAIN
USE MODULE1
IMPLICIT NONE
INTEGER                      :: K,I,J
REAL(8),DIMENSION(IM,JM)     :: X,Y,U,V,OMEGA,RHS,VEL_DIV,BXX,BYY,BXY,BYX,PHIS,SI
!------------------------------------------------------------------------------
CALL GRID (X,Y)
CALL INITIAL_CONDITIONS (X,Y,U,V,SI,OMEGA,PHIS,BXX,BXY,BYX,BYY)
CALL UV_BC (Y,U,V)
CALL SI_BC (U,SI)
!------------------------------------------------------------------------------
DO K=KSTART+1,KMAX
CALL OUTLET_BOUNDARY_CONDITIONS (U,V,SI)
!CALL EULER_VV (K,X,Y,U,V,OMEGA,PHIS,BXX,BXY,BYX,BYY)
CALL EULER_VS (K,X,Y,U,V,SI,OMEGA,PHIS,BXX,BXY,BYX,BYY)
!PRINT*,K
END DO
!------------------------------------------------------------------------------
END PROGRAM MAIN
!==============================================================================
SUBROUTINE GRID (X,Y)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: X,Y
!------------------------------------------------------------------------------
Y(:,1)=-0.5D0*H
DO I=1,IM
DO J=1,JM-1
Y(I,J+1)=Y(I,J)+DY
END DO
END DO
!------------------------------------------------------------------------------
X(1,:)=-0.5D0*L
DO I=1,IM-1
DO J=1,JM
X(I+1,J)=X(I,J)+DX
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE GRID
!==============================================================================
SUBROUTINE INITIAL_CONDITIONS (X,Y,U,V,SI,OMEGA,PHIS,BXX,BXY,BYX,BYY)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: OMEGA,BXX,BXY,BYX,BYY,PHIS,U,V,X,Y,SI
REAL(8)                      :: DUMMY
!------------------------------------------------------------------------------
IF (ICFLAG) THEN
OPEN(11230,FILE='FLD   9000.DAT',POSITION='REWIND')
READ(11230,*)
READ(11230,*)
DO J=1,JM
DO I=1,IM
READ(11230,*)    X(I,J),Y(I,J),U(I,J),V(I,J),SI(I,J),OMEGA(I,J),PHIS(I,J),BXX(I,J),BXY(I,J),BYX(I,J),BYY(I,J),DUMMY
END DO
END DO
PAUSE "==File is read=="
CLOSE(11230)
ELSE
DO I=1,IM
DO J=1,JM
OMEGA(I,J)=0.D0
!IF (((X(I,J)-XC)**2.D0+(Y(I,J)-YC)**2.D0).LE.R**2.D0) THEN
IF ((((X(I,J)-XC)**2.D0+(Y(I,J)-YC)**2.D0).LE.RO**2.D0).AND.(((X(I,J)-XC)**2.D0+(Y(I,J)-YC)**2.D0).GE.RI**2.D0)) THEN!
U(I,J)=0.D0
V(I,J)=0.D0
PHIS(I,J)=1.D0
BXX(I,J)=1.D0*PHIS(I,J)**0.5D0
BXY(I,J)=0.D0*PHIS(I,J)**0.5D0
BYX(I,J)=0.D0*PHIS(I,J)**0.5D0
BYY(I,J)=1.D0*PHIS(I,J)**0.5D0
END IF
END DO
END DO
END IF
!------------------------------------------------------------------------------
END SUBROUTINE INITIAL_CONDITIONS
!==============================================================================
SUBROUTINE UV_BC (Y,U,V)
USE MODULE1
IMPLICIT NONE
REAL(8),DIMENSION(IM,JM)     :: U,V,Y
!------------------------------------------------------------------------------
V=0.D0
U(:,1)=0.D0
U(:,JM)=0.D0
U(:,2:JM-1)=UMAX*(1.D0-(2.D0*Y(:,2:JM-1)/H)**2.D0)
!------------------------------------------------------------------------------
END SUBROUTINE UV_BC
!==============================================================================
SUBROUTINE SI_BC (U,SI)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: SI,U
!------------------------------------------------------------------------------
SI(:,1)=0.D0
DO J=1,JM-1
SI(:,J+1)=SI(:,J)+U(:,J)*DY
END DO
SI(:,JM)=SI(1,JM)
!------------------------------------------------------------------------------
END SUBROUTINE SI_BC
!==============================================================================
SUBROUTINE WRITE_TOTAL_SOL_VS (X,Y,U,V,SI,OMEGA,PHIS,BXX,BXY,BYX,BYY,VEL_DIV,K)
USE MODULE1
IMPLICIT NONE
INTEGER                             :: I,J,K,COUNTER
REAL(8),DIMENSION(IM,JM)            :: U,V,SI,X,Y,VEL_DIV,OMEGA,PHIS,BXX,BXY,BYX,BYY
CHARACTER*7                         :: EXT
CHARACTER*12                        :: FN1
CHARACTER*27                        :: FNAME
!------------------------------------------------------------------------------
FN1='/results/FLD'
COUNTER=0
!------------------------------------------------------------------------------
COUNTER=COUNTER+10
WRITE(EXT,'(I7)') K
FNAME=FN1//EXT//'.DAT'
!------------------------------------------------------------------------------
OPEN(COUNTER,FILE=FNAME,POSITION='REWIND')
WRITE(COUNTER,*)'VARIABLES= "X", "Y", "U", "V", "SI", "OMEGA","PHIS","BXX","BXY","BYX","BYY","VEL_DIV"'
WRITE(COUNTER,*)'ZONE, I=',IM,',J=',JM,',F=POINT'
DO J=1,JM
DO I=1,IM
WRITE(COUNTER,*) X(I,J),Y(I,J),U(I,J),V(I,J),SI(I,J),OMEGA(I,J),PHIS(I,J),BXX(I,J),BXY(I,J),BYX(I,J),BYY(I,J),VEL_DIV(I,J)
END DO
END DO
CLOSE(COUNTER)
!------------------------------------------------------------------------------
!WRITE(*,*) '========================'
!WRITE(*,*) 'PRINTING ON ',FNAME
!WRITE(*,*) '========================'
!------------------------------------------------------------------------------
END SUBROUTINE WRITE_TOTAL_SOL_VS
!==============================================================================
SUBROUTINE WRITE_INITIAL (X,Y,U,V,SI,OMEGA,PHIS,BXX,BXY,BYX,BYY,VEL_DIV)
USE MODULE1
IMPLICIT NONE
INTEGER                             :: I,J
REAL(8),DIMENSION(IM,JM)            :: U,V,SI,X,Y,VEL_DIV,OMEGA,PHIS,BXX,BXY,BYX,BYY
!------------------------------------------------------------------------------
OPEN(1,FILE="INITIAL.DAT",POSITION='REWIND')
WRITE(1,*)'VARIABLES= "X", "Y", "U", "V", "SI", "OMEGA","PHIS","BXX","BXY","BYX","BYY","VEL_DIV"'
WRITE(1,*)'ZONE, I=',IM,',J=',JM,',F=POINT'
DO J=1,JM
DO I=1,IM
WRITE(1,*) X(I,J),Y(I,J),U(I,J),V(I,J),SI(I,J),OMEGA(I,J),PHIS(I,J),BXX(I,J),BXY(I,J),BYX(I,J),BYY(I,J),VEL_DIV(I,J)
END DO
END DO
CLOSE(1)
!------------------------------------------------------------------------------
END SUBROUTINE WRITE_INITIAL
!==============================================================================
SUBROUTINE WRITE_TOTAL_SOL_VV (X,Y,U,V,OMEGA,PHIS,BXX,BXY,BYX,BYY,VEL_DIV,K)
USE MODULE1
IMPLICIT NONE
INTEGER                             :: I,J,K,COUNTER
REAL(8),DIMENSION(IM,JM)            :: U,V,X,Y,VEL_DIV,OMEGA,PHIS,BXX,BXY,BYX,BYY
CHARACTER*7                         :: EXT
CHARACTER*3                         :: FN1
CHARACTER*18                        :: FNAME
!------------------------------------------------------------------------------
FN1='FLD'
COUNTER=0
!------------------------------------------------------------------------------
COUNTER=COUNTER+10
WRITE(EXT,'(I7)') K
FNAME=FN1//EXT//'.DAT'
!------------------------------------------------------------------------------
OPEN(COUNTER,FILE=FNAME,POSITION='REWIND')
WRITE(COUNTER,*)'VARIABLES= "X", "Y", "U", "V","OMEGA","PHIS","BXX","BXY","BYX","BYY","VEL_DIV"'
WRITE(COUNTER,*)'ZONE, I=',IM,',J=',JM,',F=POINT'
DO J=1,JM
DO I=1,IM
WRITE(COUNTER,*) X(I,J),Y(I,J),U(I,J),V(I,J),OMEGA(I,J),PHIS(I,J),BXX(I,J),BXY(I,J),BYX(I,J),BYY(I,J),VEL_DIV(I,J)
END DO
END DO
CLOSE(COUNTER)
!------------------------------------------------------------------------------
WRITE(*,*) '========================'
WRITE(*,*) 'PRINTING ON ',FNAME
WRITE(*,*) '========================'
!------------------------------------------------------------------------------
END SUBROUTINE WRITE_TOTAL_SOL_VV
!==============================================================================
SUBROUTINE VOR_DEFINITION (ORDER_FLAG,U,V,OMEGA)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J,ORDER_FLAG
REAL(8),DIMENSION(IM,JM)     :: U,V,OMEGA
!------------------------------------------------------------------------------
IF (ORDER_FLAG==2) THEN
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
OMEGA(I,J)=(V(I+1,J)-V(I-1,J))/(2.D0*DX)-(U(I,J+1)-U(I,J-1))/(2.D0*DY)
END DO
END DO
!------------------------------------------------------------------------------
ELSE IF (ORDER_FLAG==4) THEN
!------------------------------------------------------------------------------
OMEGA(2,2)=(-3.D0*V(1,2)-10.D0*V(2,2)+18.D0*V(3,2)-6.D0*V(4,2)+V(5,2))/(12.D0*DX)&
          -(-3.D0*U(2,1)-10.D0*U(2,2)+18.D0*U(2,3)-6.D0*U(2,4)+U(2,5))/(12.D0*DY)
!------------------------------------------------------------------------------
OMEGA(IM-1,2)=(-V(IM-4,2)+6.D0*V(IM-3,2)-18.D0*V(IM-2,2)+10.D0*V(IM-1,2)+3.D0*V(IM,2))/(12.D0*DX)&
             -(-3.D0*U(IM-1,1)-10.D0*U(IM-1,2)+18.D0*U(IM-1,3)-6.D0*U(IM-1,4)+U(IM-1,5))/(12.D0*DY)
!------------------------------------------------------------------------------
OMEGA(2,JM-1)=(-3.D0*V(1,JM-1)-10.D0*V(2,JM-1)+18.D0*V(3,JM-1)-6.D0*V(4,JM-1)+V(5,JM-1))/(12.D0*DX)&
             -(-U(2,JM-4)+6.D0*U(2,JM-3)-18.D0*U(2,JM-2)+10.D0*U(2,JM-1)+3.D0*U(2,JM))/(12.D0*DY)
!------------------------------------------------------------------------------
OMEGA(IM-1,JM-1)=(-V(IM-4,JM-1)+6.D0*V(IM-3,JM-1)-18.D0*V(IM-2,JM-1)+10.D0*V(IM-1,JM-1)+3.D0*V(IM,JM-1))/(12.D0*DX)&
                -(-U(IM-1,JM-4)+6.D0*U(IM-1,JM-3)-18.D0*U(IM-1,JM-2)+10.D0*U(IM-1,JM-1)+3.D0*U(IM-1,JM))/(12.D0*DY)
!------------------------------------------------------------------------------
DO J=3,JM-2
OMEGA(2,J)=(-3.D0*V(1,J)-10.D0*V(2,J)+18.D0*V(3,J)-6.D0*V(4,J)+V(5,J))/(12.D0*DX)&
          -(U(2,J-2)-8.D0*U(2,J-1)+8.D0*U(2,J+1)-U(2,J+2))/(12.D0*DY)
END DO
!------------------------------------------------------------------------------
DO J=3,JM-2
OMEGA(IM-1,J)=(-V(IM-4,J)+6.D0*V(IM-3,J)-18.D0*V(IM-2,J)+10.D0*V(IM-1,J)+3.D0*V(IM,J))/(12.D0*DX)&
             -(U(IM-1,J-2)-8.D0*U(IM-1,J-1)+8.D0*U(IM-1,J+1)-U(IM-1,J+2))/(12.D0*DY)
END DO
!------------------------------------------------------------------------------
DO I=3,IM-2
OMEGA(I,JM-1)=(V(I-2,JM-1)-8.D0*V(I-1,JM-1)+8.D0*V(I+1,JM-1)-V(I+2,JM-1))/(12.D0*DX)&
             -(-U(I,JM-4)+6.D0*U(I,JM-3)-18.D0*U(I,JM-2)+10.D0*U(I,JM-1)+3.D0*U(I,JM))/(12.D0*DY)
END DO
!------------------------------------------------------------------------------
DO I=3,IM-2
OMEGA(I,2)=(V(I-2,2)-8.D0*V(I-1,2)+8.D0*V(I+1,2)-V(I+2,2))/(12.D0*DX)&
          -(-3.D0*U(I,1)-10.D0*U(I,2)+18.D0*U(I,3)-6.D0*U(I,4)+U(I,5))/(12.D0*DY)
END DO
!------------------------------------------------------------------------------
DO I=3,IM-2
DO J=3,JM-2
OMEGA(I,J)=(V(I-2,J)-8.D0*V(I-1,J)+8.D0*V(I+1,J)-V(I+2,J))/(12.D0*DX)&
          -(U(I,J-2)-8.D0*U(I,J-1)+8.D0*U(I,J+1)-U(I,J+2))/(12.D0*DY)
END DO
END DO
!------------------------------------------------------------------------------
ELSE
!------------------------------------------------------------------------------
PRINT*,'!!!Invalid input for order of accuracy please revise the ORDER_FLAG'
!------------------------------------------------------------------------------
END IF
!------------------------------------------------------------------------------
END SUBROUTINE VOR_DEFINITION
!==============================================================================
SUBROUTINE OMEGA_BOUNDARY_CONDITIONS (ORDER_FLAG,U,V,OMEGA)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: ORDER_FLAG,I,J
REAL(8),DIMENSION(IM,JM)     :: OMEGA,U,V
!------------------------------------------------------------------------------
IF (ORDER_FLAG==2) THEN
!------------------------------------------------------------------------------
DO I=2,IM-1
OMEGA(I,1)=(V(I-1,1)-2.D0*V(I,1)+V(I+1,1))/(2.D0*DX)-(-3.D0*U(I,1)+4.D0*U(I,2)-U(I,3))/(2.D0*DY)
OMEGA(I,JM)=(V(I-1,JM)-2.D0*V(I,JM)+V(I+1,JM))/(2.D0*DX)-(U(I,JM-2)-4.D0*U(I,JM-1)+3.D0*U(I,JM))/(2.D0*DY)
END DO
DO J=2,JM-1
OMEGA(1,J)=(-3.D0*V(1,J)+4.D0*V(2,J)-V(3,J))/(2.D0*DX)-(U(1,J-1)-2.D0*U(1,J)+U(1,J+1))/(2.D0*DY)
OMEGA(IM,J)=(V(IM-2,J)-4.D0*V(IM-1,J)+3.D0*V(IM,J))/(2.D0*DX)-(U(IM,J-1)-2.D0*U(IM,J)+U(IM,J+1))/(2.D0*DY)
END DO
!------------------------------------------------------------------------------
ELSE IF (ORDER_FLAG==4) THEN
!------------------------------------------------------------------------------
OMEGA(1,1)=(-25.D0*V(1,1)+48.D0*V(2,1)-36.D0*V(3,1)+16.D0*V(4,1)-3.D0*V(5,1))/(12.D0*DX)&
          -(-25.D0*U(1,1)+48.D0*U(1,2)-36.D0*U(1,3)+16.D0*U(1,4)-3.D0*U(1,5))/(12.D0*DY)
OMEGA(IM,1)=(3.D0*V(IM-4,1)-16.D0*V(IM-3,1)+36.D0*V(IM-2,1)-48.D0*V(IM-1,1)+25.D0*V(IM,1))/(12.D0*DX)&
           -(-25.D0*U(IM,1)+48.D0*U(IM,2)-36.D0*U(IM,3)+16.D0*U(IM,4)-3.D0*U(IM,5))/(12.D0*DY)
OMEGA(1,JM)=(-25.D0*V(1,JM)+48.D0*V(2,JM)-36.D0*V(3,JM)+16.D0*V(4,JM)-3.D0*V(5,JM))/(12.D0*DX)&
           -(3.D0*U(1,JM-4)-16.D0*U(1,JM-3)+36.D0*U(1,JM-2)-48.D0*U(1,JM-1)+25.D0*U(1,JM))/(12.D0*DY)
OMEGA(IM,JM)=(3.D0*V(IM-4,JM)-16.D0*V(IM-3,JM)+36.D0*V(IM-2,JM)-48.D0*V(IM-1,JM)+25.D0*V(IM,JM))/(12.D0*DX)&
            -(3.D0*U(IM,JM-4)-16.D0*U(IM,JM-3)+36.D0*U(IM,JM-2)-48.D0*U(IM,JM-1)+25.D0*U(IM,JM))/(12.D0*DY)
!------------------------------------------------------------------------------
OMEGA(2,1)=(-3.D0*V(1,1)-10.D0*V(2,1)+18.D0*V(3,1)-6.D0*V(4,1)+V(5,1))/(12.D0*DX)&
          -(-25.D0*U(2,1)+48.D0*U(2,2)-36.D0*U(2,3)+16.D0*U(2,4)-3.D0*U(2,5))/(12.D0*DY)
OMEGA(IM-1,1)=(-V(IM-4,1)+6.D0*V(IM-3,1)-18.D0*V(IM-2,1)+10.D0*V(IM-1,1)+3.D0*V(IM,1))/(12.D0*DX)&
             -(-25.D0*U(IM-1,1)+48.D0*U(IM-1,2)-36.D0*U(IM-1,3)+16.D0*U(IM-1,4)-3.D0*U(IM-1,5))/(12.D0*DY)
OMEGA(1,2)=(-25.D0*V(1,2)+48.D0*V(2,2)-36.D0*V(3,2)+16.D0*V(4,2)-3.D0*V(5,2))/(12.D0*DX)&
          -(-3.D0*U(1,1)-10.D0*U(1,2)+18.D0*U(1,3)-6.D0*U(1,4)+U(1,5))/(12.D0*DY)
OMEGA(IM,2)=(3.D0*V(IM-4,2)-16.D0*V(IM-3,2)+36.D0*V(IM-2,2)-48.D0*V(IM-1,2)+25.D0*V(IM,2))/(12.D0*DX)&
           -(-3.D0*U(IM,1)-10.D0*U(IM,2)+18.D0*U(IM,3)-6.D0*U(IM,4)+U(IM,5))/(12.D0*DY)
OMEGA(1,JM-1)=(-25.D0*V(1,JM-1)+48.D0*V(2,JM-1)-36.D0*V(3,JM-1)+16.D0*V(4,JM-1)-3.D0*V(5,JM-1))/(12.D0*DX)&
             -(-U(1,JM-4)+6.D0*U(1,JM-3)-18.D0*U(1,JM-2)+10.D0*U(1,JM-1)+3.D0*U(1,JM))/(12.D0*DY)
OMEGA(IM,JM-1)=(3.D0*V(IM-4,JM-1)-16.D0*V(IM-3,JM-1)+36.D0*V(IM-2,JM-1)-48.D0*V(IM-1,JM-1)+25.D0*V(IM,JM-1))/(12.D0*DX)&
              -(-U(IM,JM-4)+6.D0*U(IM,JM-3)-18.D0*U(IM,JM-2)+10.D0*U(IM,JM-1)+3.D0*U(IM,JM))/(12.D0*DY)
OMEGA(IM-1,JM)=(-V(IM-4,JM)+6.D0*V(IM-3,JM)-18.D0*V(IM-2,JM)+10.D0*V(IM-1,JM)+3.D0*V(IM,JM))/(12.D0*DX)&
              -(3.D0*U(IM-1,JM-4)-16.D0*U(IM-1,JM-3)+36.D0*U(IM-1,JM-2)-48.D0*U(IM-1,JM-1)+25.D0*U(IM-1,JM))/(12.D0*DY)
OMEGA(2,JM)=(-3.D0*V(1,JM)-10.D0*V(2,JM)+18.D0*V(3,JM)-6.D0*V(4,JM)+V(5,JM))/(12.D0*DX)&
           -(3.D0*U(2,JM-4)-16.D0*U(2,JM-3)+36.D0*U(2,JM-2)-48.D0*U(2,JM-1)+25.D0*U(2,JM))/(12.D0*DY)
!------------------------------------------------------------------------------
DO I=3,IM-2
OMEGA(I,1)=(V(I-2,1)-8.D0*V(I-1,1)+8.D0*V(I+1,1)-V(I+2,1))/(12.D0*DX)&
          -(-25.D0*U(I,1)+48.D0*U(I,2)-36.D0*U(I,3)+16.D0*U(I,4)-3.D0*U(I,5))/(12.D0*DY)
OMEGA(I,JM)=(V(I-2,JM)-8.D0*V(I-1,JM)+8.D0*V(I+1,JM)-V(I+2,JM))/(12.D0*DX)&
           -(3.D0*U(I,JM-4)-16.D0*U(I,JM-3)+36.D0*U(I,JM-2)-48.D0*U(I,JM-1)+25.D0*U(I,JM))/(12.D0*DY)
END DO
DO J=3,JM-2
OMEGA(1,J)=(-25.D0*V(1,J)+48.D0*V(2,J)-36.D0*V(3,J)+16.D0*V(4,J)-3.D0*V(5,J))/(12.D0*DX)&
          -(U(1,J-2)-8.D0*U(1,J-1)+8.D0*U(1,J+1)-U(1,J+2))/(12.D0*DY)
OMEGA(IM,J)=(3.D0*V(IM-4,J)-16.D0*V(IM-3,J)+36.D0*V(IM-2,J)-48.D0*V(IM-1,J)+25.D0*V(IM,J))/(12.D0*DX)&
           -(U(IM,J-2)-8.D0*U(IM,J-1)+8.D0*U(IM,J+1)-U(IM,J+2))/(12.D0*DY)
END DO
!------------------------------------------------------------------------------
ELSE
!------------------------------------------------------------------------------
PRINT*,'!!!Invalid input for order of accuracy please revise the ORDER_FLAG'
!------------------------------------------------------------------------------
END IF
!------------------------------------------------------------------------------
END SUBROUTINE OMEGA_BOUNDARY_CONDITIONS
!==============================================================================
SUBROUTINE VELOCITY_DIVERGENCE (U,V,VEL_DIV)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: U,V,VEL_DIV
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
VEL_DIV(I,J)=(U(I+1,J)-U(I-1,J))/(2.D0*DX)+(V(I,J+1)-V(I,J-1))/(2.D0*DY)
END DO
END DO
!------------------------------------------------------------------------------
!Lower Boundary
DO I=2,IM-1
VEL_DIV(I,1)=(U(I+1,1)-U(I-1,1))/(2.D0*DX)+(-3.D0*V(I,1)+4.D0*V(I,2)-V(I,3))/(2.D0*DY)
END DO
!------------------------------------------------------------------------------
!Upper Boundary
DO I=2,IM-1
VEL_DIV(I,JM)=(U(I+1,JM)-U(I-1,JM))/(2.D0*DX)+(3.D0*V(I,JM)-4.D0*V(I,JM-1)+V(I,JM-2))/(2.D0*DY)
END DO
!------------------------------------------------------------------------------
!Left Boundary
DO J=2,JM-1
VEL_DIV(1,J)=(-3.D0*U(1,J)+4.D0*U(2,J)-U(3,J))/(2.D0*DX)+(V(1,J+1)-V(1,J-1))/(2.D0*DY)
END DO
!------------------------------------------------------------------------------
!Right Boundary
DO J=2,JM-1
VEL_DIV(IM,J)=(3.D0*U(IM,J)-4.D0*U(IM-1,J)+U(IM-2,J))/(2.D0*DX)+(V(IM,J+1)-V(IM,J-1))/(2.D0*DY)
END DO
!------------------------------------------------------------------------------
!Upper Right
VEL_DIV(IM,JM)=(3.D0*U(IM,JM)-4.D0*U(IM-1,JM)+U(IM-2,JM))/(2.D0*DX)+(3.D0*V(IM,JM)-4.D0*V(IM,JM-1)+V(IM,JM-2))/(2.D0*DY)
!------------------------------------------------------------------------------
!Lower Right
VEL_DIV(IM,1)=(3.D0*U(IM,1)-4.D0*U(IM-1,1)+U(IM-2,1))/(2.D0*DX)+(-3.D0*V(IM,1)+4.D0*V(IM,2)-V(IM,3))/(2.D0*DY)
!------------------------------------------------------------------------------
!Upper Left
VEL_DIV(1,JM)=(-3.D0*U(1,JM)+4.D0*U(2,JM)-U(3,JM))/(2.D0*DX)+(3.D0*V(1,JM)-4.D0*V(1,JM-1)+V(1,JM-2))/(2.D0*DY)
!------------------------------------------------------------------------------
!Lower Left
VEL_DIV(1,1)=(-3.D0*V(1,1)+4.D0*V(1,2)-V(1,3))/(2.D0*DY)+(-3.D0*U(1,1)+4.D0*U(2,1)-U(3,1))/(2.D0*DX)
!------------------------------------------------------------------------------
END SUBROUTINE VELOCITY_DIVERGENCE
!==============================================================================
SUBROUTINE GS_FOURTH (U,F)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J,ITER
REAL(8),DIMENSION(IM,JM)     :: F,U,UP
REAL(8)                      :: ERROR
!------------------------------------------------------------------------------
ITER=0
ERROR=1.D0
DO WHILE (ERROR.GT.EPS.AND.ITER.LE.ITERM)
ITER=ITER+1
ERROR=0.D0
UP=U
!------------------------------------------------------------------------------
U(2,2)=(10.D0*U(1,2)-4.D0*U(3,2)+14.D0*U(4,2)-6.D0*U(5,2)+U(6,2)&
      +BETA2*(10.D0*U(2,1)-4.D0*U(2,3)+14.D0*U(2,4)-6.D0*U(2,5)+U(2,6))&
	  -12.D0*DX2*F(2,2))&
	  /(15.D0*(1.D0+BETA2))
!------------------------------------------------------------------------------
U(IM-1,2)=(U(IM-5,2)-6.D0*U(IM-4,2)+14.D0*U(IM-3,2)-4.D0*U(IM-2,2)+10.D0*U(IM,2)&
      +BETA2*(10.D0*U(IM-1,1)-4.D0*U(IM-1,3)+14.D0*U(IM-1,4)-6.D0*U(IM-1,5)+U(IM-1,6))&
	  -12.D0*DX2*F(IM-1,2))&
	  /(15.D0*(1.D0+BETA2))
!------------------------------------------------------------------------------
U(IM-1,JM-1)=(U(IM-5,JM-1)-6.D0*U(IM-4,JM-1)+14.D0*U(IM-3,JM-1)-4.D0*U(IM-2,JM-1)+10.D0*U(IM,JM-1)&
      +BETA2*(U(IM-1,JM-5)-6.D0*U(IM-1,JM-4)+14.D0*U(IM-1,JM-3)-4.D0*U(IM-1,JM-2)+10.D0*U(IM-1,JM))&
	  -12.D0*DX2*F(IM-1,JM-1))&
	  /(15.D0*(1.D0+BETA2))
!------------------------------------------------------------------------------
U(2,JM-1)=(10.D0*U(1,JM-1)-4.D0*U(3,JM-1)+14.D0*U(4,JM-1)-6.D0*U(5,JM-1)+U(6,JM-1)&
      +BETA2*(U(2,JM-5)-6.D0*U(2,JM-4)+14.D0*U(2,JM-3)-4.D0*U(2,JM-2)+10.D0*U(2,JM))&
	  -12.D0*DX2*F(2,JM-1))&
	  /(15.D0*(1.D0+BETA2))
!------------------------------------------------------------------------------
DO J=3,JM-2
U(2,J)=(10.D0*U(1,J)-4.D0*U(3,J)+14.D0*U(4,J)-6.D0*U(5,J)+U(6,J)&
      +BETA2*(-U(2,J-2)+16.D0*U(2,J-1)+16.D0*U(2,J+1)-U(2,J+2))&
	  -12.D0*DX2*F(2,J))&
	  /(15.D0*(1.D0+2.D0*BETA2))
END DO
!------------------------------------------------------------------------------
DO J=3,JM-2
U(IM-1,J)=(U(IM-5,J)-6.D0*U(IM-4,J)+14.D0*U(IM-3,J)-4.D0*U(IM-2,J)+10.D0*U(IM,J)&
      +BETA2*(-U(IM-1,J-2)+16.D0*U(IM-1,J-1)+16.D0*U(IM-1,J+1)-U(IM-1,J+2))&
	  -12.D0*DX2*F(IM-1,J))&
	  /(15.D0*(1.D0+2.D0*BETA2))
END DO
!------------------------------------------------------------------------------
DO I=3,IM-2
U(I,2)=(-U(I-2,2)+16.D0*U(I-1,2)+16.D0*U(I+1,2)-U(I+2,2)&
      +BETA2*(10.D0*U(I,1)-4.D0*U(I,3)+14.D0*U(I,4)-6.D0*U(I,5)+U(I,6))&
	  -12.D0*DX2*F(I,2))&
	  /(15.D0*(2.D0+BETA2))
END DO
!------------------------------------------------------------------------------
DO I=3,IM-2
U(I,JM-1)=(-U(I-2,JM-1)+16.D0*U(I-1,JM-1)+16.D0*U(I+1,JM-1)-U(I+2,JM-1)&
         +BETA2*(U(I,JM-5)-6.D0*U(I,JM-4)+14.D0*U(I,JM-3)-4.D0*U(I,JM-2)+10.D0*U(I,JM))&
		 -12.D0*DX2*F(I,JM-1))&
		 /(15.D0*(2.D0+BETA2))
END DO
!------------------------------------------------------------------------------
DO I=3,IM-2
DO J=3,JM-2
U(I,J)=(-U(I-2,J)+16.D0*U(I-1,J)+16.D0*U(I+1,J)-U(I+2,J)&
      +BETA2*(-U(I,J-2)+16.D0*U(I,J-1)+16.D0*U(I,J+1)-U(I,J+2))&
	  -12.D0*DX2*F(I,J))&
	  /(30.D0*(1.D0+BETA2))
END DO
END DO
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
ERROR=MAX(ERROR,ABS(U(I,J)-UP(I,J)))
END DO
END DO
!------------------------------------------------------------------------------
END DO
!------------------------------------------------------------------------------
IF (ERROR.GT.EPS) THEN
PRINT*, "Convergence can not be met with",ITERM,"iterations."
END IF
!------------------------------------------------------------------------------
END SUBROUTINE GS_FOURTH
!==============================================================================
SUBROUTINE GS_SECOND (U,F)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J,ITER
REAL(8),DIMENSION(IM,JM)     :: F,U,UP
REAL(8)                      :: ERROR
!------------------------------------------------------------------------------
ITER=0
ERROR=1.D0
DO WHILE (ERROR.GT.EPS.AND.ITER.LT.ITERM)
ITER=ITER+1
ERROR=0.D0
UP=U
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
U(I,J)=(U(I-1,J)+U(I+1,J)+BETA2*(U(I,J-1)+U(I,J+1))-DX2*F(I,J))&
	  /(2.D0*(1.D0+BETA2))
END DO
END DO
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
ERROR=MAX(ERROR,ABS(U(I,J)-UP(I,J)))
END DO
END DO
!------------------------------------------------------------------------------
END DO
!------------------------------------------------------------------------------
IF (ITER.GT.ITERM.AND.ERROR.GT.EPS) THEN
PRINT*, "Convergence can not be met with",ITERM-1,"iterations."
END IF
!------------------------------------------------------------------------------
END SUBROUTINE GS_SECOND
!==============================================================================
SUBROUTINE EULER_VS (K,X,Y,U,V,SI,OMEGA,PHIS,BXX,BXY,BYX,BYY)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: K,I,J
REAL(8),DIMENSION(IM,JM)     :: U,V,OMEGA_RHS,OMEGA,F,F1,F2,X,Y,VEL_DIV,BXX,BXY,BYX,BYY
REAL(8),DIMENSION(IM,JM)     :: BXX_RHS,BXY_RHS,BYX_RHS,BYY_RHS,PHIS,SI,SIGMA_SH_XX,SIGMA_SH_XY,SIGMA_SH_YX,SIGMA_SH_YY,CDHST
REAL(8)                      :: UNORM
!------------------------------------------------------------------------------
CALL FORCING_FUNCTION_VS (OMEGA,F)
CALL SOLVER (SI,F)
CALL U_V (SI,U,V)
CALL ZERO_FIELD_VELOCITY (X,Y,PHIS,U,V)
CALL VOR_DEFINITION (2,U,V,OMEGA)
CALL FORCING_FUNCTION_VS (OMEGA,F)
CALL SOLVER (SI,F)
CALL U_V (SI,U,V)
CALL OMEGA_BOUNDARY_CONDITIONS (2,U,V,OMEGA)
CALL DEVIATORIC_HYPERELASTIC_TENSOR (PHIS,BXX,BXY,BYX,BYY,SIGMA_SH_XX,SIGMA_SH_XY,SIGMA_SH_YX,SIGMA_SH_YY)
CALL CURL_DIVERGENCE_HYPERELASTIC_STRESS_TENSOR (SIGMA_SH_XX,SIGMA_SH_XY,SIGMA_SH_YX,SIGMA_SH_YY,CDHST)
CALL VOR_TRANSPORT_RHS (PHIS,U,V,OMEGA,CDHST,OMEGA_RHS)
CALL LEFT_CAUCHY_GREEN_DEFORMATION (U,V,PHIS,BXX,BXY,BYX,BYY)
CALL VELOCITY_DIVERGENCE (U,V,VEL_DIV)
IF (WRITE_FIELD_STATUS) THEN
IF (((K/NPRINT1)*NPRINT1)==K) CALL WRITE_TOTAL_SOL_VS (X,Y,U,V,SI,OMEGA,PHIS,BXX,BXY,BYX,BYY,VEL_DIV,K)
END IF
CALL LEVEL_SET_FUNCTION_MCCORMACK (U,V,PHIS)
!CALL LEVEL_SET_FUNCTION_LAXWENDROFF (U,V,PHIS)
CALL CONVERGENCE_TEST (K,PHIS,U,V,UNORM)
OMEGA=OMEGA+DT*OMEGA_RHS
!------------------------------------------------------------------------------
END SUBROUTINE EULER_VS
!==============================================================================
SUBROUTINE EULER_VV (K,X,Y,U,V,OMEGA,PHIS,BXX,BXY,BYX,BYY)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: K,I,J
REAL(8),DIMENSION(IM,JM)     :: U,V,OMEGA_RHS,OMEGA,F,F1,F2,X,Y,VEL_DIV,UP,VP,BXX,BXY,BYX,BYY
REAL(8),DIMENSION(IM,JM)     :: BXX_RHS,BXY_RHS,BYX_RHS,BYY_RHS,PHIS,SI,SIGMA_SH_XX,SIGMA_SH_XY,SIGMA_SH_YX,SIGMA_SH_YY,CDHST
REAL(8)                      :: UNORM
!------------------------------------------------------------------------------
UP=U
VP=V
CALL FORCING_FUNCTION_VV (4,OMEGA,F1,F2)
CALL SOLVER (U,F1)
CALL SOLVER (V,F2)
CALL ZERO_FIELD_VELOCITY (X,Y,PHIS,UP,VP,U,V)
CALL VOR_DEFINITION (2,U,V,OMEGA)
CALL FORCING_FUNCTION_VV (4,OMEGA,F1,F2)
CALL SOLVER (U,F1)
CALL SOLVER (V,F2)
CALL OMEGA_BOUNDARY_CONDITIONS (4,U,V,OMEGA)
CALL DEVIATORIC_HYPERELASTIC_TENSOR (PHIS,BXX,BXY,BYX,BYY,SIGMA_SH_XX,SIGMA_SH_XY,SIGMA_SH_YX,SIGMA_SH_YY)
CALL CURL_DIVERGENCE_HYPERELASTIC_STRESS_TENSOR (SIGMA_SH_XX,SIGMA_SH_XY,SIGMA_SH_YX,SIGMA_SH_YY,CDHST)
CALL VOR_TRANSPORT_RHS (PHIS,U,V,OMEGA,CDHST,OMEGA_RHS)
CALL LEFT_CAUCHY_GREEN_DEFORMATION (U,V,PHIS,BXX,BXY,BYX,BYY)
CALL VELOCITY_DIVERGENCE (U,V,VEL_DIV)
IF (WRITE_FIELD_STATUS) THEN
IF (((K/NPRINT1)*NPRINT1)==K) CALL WRITE_TOTAL_SOL_VV (X,Y,U,V,OMEGA,PHIS,BXX,BXY,BYX,BYY,VEL_DIV,K)
END IF
CALL LEVEL_SET_FUNCTION_MCCORMACK (U,V,PHIS)
!CALL LEVEL_SET_FUNCTION_LAXWENDROFF (U,V,PHIS)
CALL CONVERGENCE_TEST (K,PHIS,U,V,UNORM)
OMEGA=OMEGA+DT*OMEGA_RHS
!------------------------------------------------------------------------------
END SUBROUTINE EULER_VV
!==============================================================================
SUBROUTINE CONVERGENCE_TEST (K,PHIS,U,V,UNORM)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: K,I,J,N
REAL(8)                      :: UNORM,DIFF,UNORMP
REAL(8),DIMENSION(IM,JM)     :: U,V,PHIS
!------------------------------------------------------------------------------
UNORM=0.D0
N=0
DO I=1,IM
DO J=1,JM
IF (PHIS(I,J).GE.0.5D0) THEN
N=N+1
UNORM=UNORM+(U(I,J)**2.D0+V(I,J)**2.D0)
!UNORM=MAX(UNORM,DSQRT(U(I,J)**2.D0+V(I,J)**2.D0))
END IF
END DO
END DO
!UNORM=DSQRT(UNORM)/(IMJM)
UNORM=UNORM/N
!DIFF=ABS(UNORM-UNORMP)
IF (WRITE_POINT_STATUS) THEN
IF (((K/NPRINT2)*NPRINT2)==K) CALL WRITE_VEL_NORM (UNORM,K)
END IF
!IF (DIFF.LE.1.D-10) STOP
!UNORMP=UNORM
!------------------------------------------------------------------------------
END SUBROUTINE CONVERGENCE_TEST
!==============================================================================
SUBROUTINE WRITE_VEL_NORM (UNORM,K)
USE MODULE1
IMPLICIT NONE
INTEGER                             :: K,COUNT
REAL(8)                             :: UNORM
CHARACTER*7                         :: EXT
CHARACTER*3                         :: FN1
CHARACTER*18                        :: FNAME
!------------------------------------------------------------------------------
FN1='URESULT'
!COUNT=0
!------------------------------------------------------------------------------
COUNT=COUNT+1
WRITE(EXT,'(I7)') K
FNAME=FN1//EXT//'.DAT'
!------------------------------------------------------------------------------
OPEN(COUNT,FILE=FNAME,POSITION='REWIND')
WRITE(COUNT,*) K*DT,UNORM
CLOSE(COUNT)
!------------------------------------------------------------------------------
WRITE(*,*) '========================'
WRITE(*,*) 'PRINTING ON ',FNAME
WRITE(*,*) '========================'
!------------------------------------------------------------------------------
END SUBROUTINE WRITE_VEL_NORM
!==============================================================================
SUBROUTINE SOLVER (PHI,F)
USE MODULE1
IMPLICIT NONE
REAL(8),DIMENSION (IM,JM)    :: PHI,F
!------------------------------------------------------------------------------
!CALL PSOR (PHI,F)
!CALL GS_FOURTH (PHI,F)
!CALL GS_SECOND (PHI,F)
CALL MGLIN (PHI,F)
!------------------------------------------------------------------------------
END SUBROUTINE SOLVER
!==============================================================================
SUBROUTINE PSOR (PHI,F)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J,ITER
REAL(8),DIMENSION(IM,JM)     :: PHI,PHIP,F
REAL(8)                      :: E,A,ORF
!------------------------------------------------------------------------------
A=((DCOS(PI/(IM-1))+BETA2*DCOS(PI/(JM-1)))/(1+BETA2))**2
ORF=(2.D0-2.D0*DSQRT(1-A))/A
!------------------------------------------------------------------------------
ITER=1
E=1.D0
DO WHILE (ITER<ITERM.AND.E>EPS)
!PRINT*,ITER
ITER=ITER+1
E=0.D0
PHIP=PHI
DO I=2,IM-1
DO J=2,JM-1
PHI(I,J)=(1.D0-ORF)*PHI(I,J)+ORF*(BETA2*(PHI(I,J+1)+PHI(I,J-1)+PHI(I+1,J)+PHI(I-1,J)-DX2*F(I,J))/(2.D0*(1.D0+BETA2)))
E=MAX(E,ABS(PHIP(I,J)-PHI(I,J)))
END DO
END DO
END DO
IF (E.GT.EPS) THEN
PRINT*, "!!!Convergence can not be met with", ITERM-1, "iterations..."
END IF
!------------------------------------------------------------------------------
END SUBROUTINE PSOR
!==============================================================================
SUBROUTINE FORCING_FUNCTION_VV (ORDER_FLAG,OMEGA,F1,F2)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J,ORDER_FLAG
REAL(8),DIMENSION(IM,JM)     :: F1,F2,OMEGA
!------------------------------------------------------------------------------
IF (ORDER_FLAG==2) THEN
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
F1(I,J)=-(OMEGA(I,J+1)-OMEGA(I,J-1))/(2.D0*DY)
F2(I,J)=(OMEGA(I+1,J)-OMEGA(I-1,J))/(2.D0*DX)
END DO
END DO
!------------------------------------------------------------------------------
ELSE IF (ORDER_FLAG==4) THEN
!------------------------------------------------------------------------------
F1(2,2)=-(-3.D0*OMEGA(2,1)-10.D0*OMEGA(2,2)+18.D0*OMEGA(2,3)-6.D0*OMEGA(2,4)+OMEGA(2,5))/(12.D0*DY)
F1(IM-1,2)=-(-3.D0*OMEGA(IM-1,1)-10.D0*OMEGA(IM-1,2)+18.D0*OMEGA(IM-1,3)-6.D0*OMEGA(IM-1,4)+OMEGA(IM-1,5))/(12.D0*DY)
F1(2,JM-1)=-(-OMEGA(2,JM-4)+6.D0*OMEGA(2,JM-3)-18.D0*OMEGA(2,JM-2)+10.D0*OMEGA(2,JM-1)+3.D0*OMEGA(2,JM))/(12.D0*DY)
F1(IM-1,JM-1)=-(-OMEGA(IM-1,JM-4)+6.D0*OMEGA(IM-1,JM-3)-18.D0*OMEGA(IM-1,JM-2)+10.D0*OMEGA(IM-1,JM-1)&
              +3.D0*OMEGA(IM-1,JM))/(12.D0*DY)
DO J=3,JM-2
F1(IM-1,J)=-(OMEGA(IM-1,J-2)-8.D0*OMEGA(IM-1,J-1)+8.D0*OMEGA(IM-1,J+1)-OMEGA(IM-1,J+2))/(12.D0*DY)
F1(2,J)=-(OMEGA(2,J-2)-8.D0*OMEGA(2,J-1)+8.D0*OMEGA(2,J+1)-OMEGA(2,J+2))/(12.D0*DY)
END DO
DO I=3,IM-2
F1(I,2)=-(-3.D0*OMEGA(I,1)-10.D0*OMEGA(I,2)+18.D0*OMEGA(I,3)-6.D0*OMEGA(I,4)+OMEGA(I,5))/(12.D0*DY)
F1(I,JM-1)=-(-OMEGA(I,JM-4)+6.D0*OMEGA(I,JM-3)-18.D0*OMEGA(I,JM-2)+10.D0*OMEGA(I,JM-1)+3.D0*OMEGA(I,JM))/(12.D0*DY)
END DO
DO I=3,IM-2
DO J=3,JM-2
F1(I,J)=-(OMEGA(I,J-2)-8.D0*OMEGA(I,J-1)+8.D0*OMEGA(I,J+1)-OMEGA(I,J+2))/(12.D0*DY)
END DO
END DO
!------------------------------------------------------------------------------
F2(2,2)=(-3.D0*OMEGA(1,2)-10.D0*OMEGA(2,2)+18.D0*OMEGA(3,2)-6.D0*OMEGA(4,2)+OMEGA(5,2))/(12.D0*DX)
F2(IM-1,2)=(-OMEGA(IM-4,2)+6.D0*OMEGA(IM-3,2)-18.D0*OMEGA(IM-2,2)+10.D0*OMEGA(IM-1,2)+3.D0*OMEGA(IM,2))/(12.D0*DX)
F2(2,JM-1)=(-3.D0*OMEGA(1,JM-1)-10.D0*OMEGA(2,JM-1)+18.D0*OMEGA(3,JM-1)-6.D0*OMEGA(4,JM-1)+OMEGA(5,JM-1))/(12.D0*DX)
F2(IM-1,JM-1)=(-OMEGA(IM-4,JM-1)+6.D0*OMEGA(IM-3,JM-1)-18.D0*OMEGA(IM-2,JM-1)+10.D0*OMEGA(IM-1,JM-1)+3.D0*OMEGA(IM,JM-1))/(12.D0*DX)
DO J=3,JM-2
F2(IM-1,J)=(-OMEGA(IM-4,J)+6.D0*OMEGA(IM-3,J)-18.D0*OMEGA(IM-2,J)+10.D0*OMEGA(IM-1,J)+3.D0*OMEGA(IM,J))/(12.D0*DX)
F2(2,J)=(-3.D0*OMEGA(1,J)-10.D0*OMEGA(2,J)+18.D0*OMEGA(3,J)-6.D0*OMEGA(4,J)+OMEGA(5,J))/(12.D0*DX)
END DO
DO I=3,IM-2
F2(I,2)=(OMEGA(I-2,2)-8.D0*OMEGA(I-1,2)+8.D0*OMEGA(I+1,2)-OMEGA(I+2,2))/(12.D0*DX)
F2(I,JM-1)=(OMEGA(I-2,JM-1)-8.D0*OMEGA(I-1,JM-1)+8.D0*OMEGA(I+1,JM-1)-OMEGA(I+2,JM-1))/(12.D0*DX)
END DO
DO I=3,IM-2
DO J=3,JM-2
F2(I,J)=(OMEGA(I-2,J)-8.D0*OMEGA(I-1,J)+8.D0*OMEGA(I+1,J)-OMEGA(I+2,J))/(12.D0*DX)
END DO
END DO
!------------------------------------------------------------------------------
ELSE
!------------------------------------------------------------------------------
PRINT*,'!!!Invalid input for order of accuracy please revise the ORDER_FLAG'
!------------------------------------------------------------------------------
END IF
!------------------------------------------------------------------------------
END SUBROUTINE FORCING_FUNCTION_VV
!==============================================================================
SUBROUTINE FORCING_FUNCTION_VS (OMEGA,F)
USE MODULE1
IMPLICIT NONE
REAL(8),DIMENSION(IM,JM)     :: F,OMEGA
!------------------------------------------------------------------------------
F=-OMEGA
!------------------------------------------------------------------------------
END SUBROUTINE FORCING_FUNCTION_VS
!==============================================================================
SUBROUTINE U_V (SI,U,V)
USE MODULE1
IMPLICIT NONE
INTEGER                  :: I,J
REAL(8),DIMENSION(IM,JM) :: SI,U,V
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
U(I,J)=+(SI(I,J+1)-SI(I,J-1))/(2*DY)
V(I,J)=-(SI(I+1,J)-SI(I-1,J))/(2*DX)
END DO
END DO
!V(IM,:)=-(-4.D0*SI(IM-1,:)+SI(IM-2,:))/(2.D0*DX)
!V(1,:)=-(4.D0*SI(2,:)-SI(3,:))/(2.D0*DX)
!U(:,1)=(4.D0*SI(:,2)-SI(:,3))/(2.D0*DY)
!U(:,JM)=(-4.D0*SI(:,JM-1)+SI(:,JM-2))/(2.D0*DY)
!------------------------------------------------------------------------------
!U(2,2)=(-3.D0*SI(2,1)-10.D0*SI(2,2)+18.D0*SI(2,3)-6.D0*SI(2,4)+SI(2,5))/(12.D0*DY)
!------------------------------------------------------------------------------
!U(IM-1,2)=(-3.D0*SI(IM-1,1)-10.D0*SI(IM-1,2)+18.D0*SI(IM-1,3)-6.D0*SI(IM-1,4)+SI(IM-1,5))/(12.D0*DY)
!------------------------------------------------------------------------------
!U(2,JM-1)=(-SI(2,JM-4)+6.D0*SI(2,JM-3)-18.D0*SI(2,JM-2)+10.D0*SI(2,JM-1)+3.D0*SI(2,JM))/(12.D0*DY)
!------------------------------------------------------------------------------
!U(IM-1,JM-1)=(-SI(IM-1,JM-4)+6.D0*SI(IM-1,JM-3)-18.D0*SI(IM-1,JM-2)+10.D0*SI(IM-1,JM-1)+3.D0*SI(IM-1,JM))/(12.D0*DY)
!------------------------------------------------------------------------------
!DO J=3,JM-2
!U(2,J)=(SI(2,J-2)-8.D0*SI(2,J-1)+8.D0*SI(2,J+1)-SI(2,J+2))/(12.D0*DY)
!END DO
!------------------------------------------------------------------------------
!DO J=3,JM-2
!U(IM-1,J)=(SI(IM-1,J-2)-8.D0*SI(IM-1,J-1)+8.D0*SI(IM-1,J+1)-SI(IM-1,J+2))/(12.D0*DY)
!END DO
!------------------------------------------------------------------------------
!DO I=3,IM-2
!U(I,JM-1)=(-SI(I,JM-4)+6.D0*SI(I,JM-3)-18.D0*SI(I,JM-2)+10.D0*SI(I,JM-1)+3.D0*SI(I,JM))/(12.D0*DY)
!END DO
!------------------------------------------------------------------------------
!DO I=3,IM-2
!U(I,2)=(-3.D0*SI(I,1)-10.D0*SI(I,2)+18.D0*SI(I,3)-6.D0*SI(I,4)+SI(I,5))/(12.D0*DY)
!END DO
!------------------------------------------------------------------------------
!DO I=3,IM-2
!DO J=3,JM-2
!U(I,J)=(SI(I,J-2)-8.D0*SI(I,J-1)+8.D0*SI(I,J+1)-SI(I,J+2))/(12.D0*DY)
!END DO
!END DO
!------------------------------------------------------------------------------
!V(2,2)=-(-3.D0*SI(1,2)-10.D0*SI(2,2)+18.D0*SI(3,2)-6.D0*SI(4,2)+SI(5,2))/(12.D0*DX)
!------------------------------------------------------------------------------
!V(IM-1,2)=-(-SI(IM-4,2)+6.D0*SI(IM-3,2)-18.D0*SI(IM-2,2)+10.D0*SI(IM-1,2)+3.D0*SI(IM,2))/(12.D0*DX)
!------------------------------------------------------------------------------
!V(2,JM-1)=-(-3.D0*SI(1,JM-1)-10.D0*SI(2,JM-1)+18.D0*SI(3,JM-1)-6.D0*SI(4,JM-1)+SI(5,JM-1))/(12.D0*DX)
!------------------------------------------------------------------------------
!V(IM-1,JM-1)=-(-SI(IM-4,JM-1)+6.D0*SI(IM-3,JM-1)-18.D0*SI(IM-2,JM-1)+10.D0*SI(IM-1,JM-1)+3.D0*SI(IM,JM-1))/(12.D0*DX)
!------------------------------------------------------------------------------
!DO J=3,JM-2
!V(2,J)=-(-3.D0*SI(1,J)-10.D0*SI(2,J)+18.D0*SI(3,J)-6.D0*SI(4,J)+SI(5,J))/(12.D0*DX)
!END DO
!------------------------------------------------------------------------------
!DO J=3,JM-2
!V(IM-1,J)=-(-SI(IM-4,J)+6.D0*SI(IM-3,J)-18.D0*SI(IM-2,J)+10.D0*SI(IM-1,J)+3.D0*SI(IM,J))/(12.D0*DX)
!END DO
!------------------------------------------------------------------------------
!DO I=3,IM-2
!V(I,JM-1)=-(SI(I-2,JM-1)-8.D0*SI(I-1,JM-1)+8.D0*SI(I+1,JM-1)-SI(I+2,JM-1))/(12.D0*DX)
!END DO
!------------------------------------------------------------------------------
!DO I=3,IM-2
!V(I,2)=-(SI(I-2,2)-8.D0*SI(I-1,2)+8.D0*SI(I+1,2)-SI(I+2,2))/(12.D0*DX)
!END DO
!------------------------------------------------------------------------------
!DO I=3,IM-2
!DO J=3,JM-2
!V(I,J)=-(SI(I-2,J)-8.D0*SI(I-1,J)+8.D0*SI(I+1,J)-SI(I+2,J))/(12.D0*DX)
!END DO
!END DO
!------------------------------------------------------------------------------
END SUBROUTINE U_V
!==============================================================================
SUBROUTINE ZERO_FIELD_VELOCITY (X,Y,PHIS,U,V)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: U,V,PHIS,X,Y
!------------------------------------------------------------------------------
!U(IM1:IM2,JM1:JM2)=0.D0
!V(IM1:IM2,JM1:JM2)=0.D0
!------------------------------------------------------------------------------
DO I=1,IM
DO J=1,JM
IF ((((X(I,J))**2.D0+(Y(I,J)+1.D0)**2.D0).LE.RC**2.D0).OR.(((X(I,J))**2.D0+(Y(I,J)-1.D0)**2.D0).LE.RC**2.D0)) THEN
U(I,J)=0.D0
V(I,J)=0.D0
END IF
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE ZERO_FIELD_VELOCITY
!==============================================================================
SUBROUTINE OUTLET_BOUNDARY_CONDITIONS (U,V,SI)
USE MODULE1
IMPLICIT NONE
REAL(8),DIMENSION(IM,JM)     :: U,V,SI
!------------------------------------------------------------------------------
IF (BCFLAG) THEN
U(IM,2:JM-1)=U(IM,2:JM-1)-(UMAX*DT/DX)*(U(IM,2:JM-1)-U(IM-1,2:JM-1))
V(IM,2:JM-1)=V(IM,2:JM-1)-(UMAX*DT/DX)*(V(IM,2:JM-1)-V(IM-1,2:JM-1))
SI(IM,2:JM-1)=SI(IM-1,2:JM-1)
ELSE
U(IM,2:JM-1)=U(IM-1,2:JM-1)
V(IM,2:JM-1)=V(IM-1,2:JM-1)
SI(IM,2:JM-1)=SI(IM-1,2:JM-1)
END IF
!------------------------------------------------------------------------------
END SUBROUTINE OUTLET_BOUNDARY_CONDITIONS
!==============================================================================
SUBROUTINE LEFT_CAUCHY_GREEN_DEFORMATION (U,V,PHIS,BXX,BXY,BYX,BYY)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: U,V,BXX,BYY,BXY,BYX,BXX_RHS,BYY_RHS,BXY_RHS,BYX_RHS,PHIS
!------------------------------------------------------------------------------
DO I=1,IM
DO J=1,JM
IF (PHIS(I,J).LE.PHIS_MIN) THEN
BXX(I,J)=0.D0
BXY(I,J)=0.D0
BYX(I,J)=0.D0
BYY(I,J)=0.D0
!ELSE
!BXX(I,J)=BXX(I,J)*PHIS(I,J)**0.5D0
!BXY(I,J)=BXY(I,J)*PHIS(I,J)**0.5D0
!BYX(I,J)=BYX(I,J)*PHIS(I,J)**0.5D0
!BYY(I,J)=BYY(I,J)*PHIS(I,J)**0.5D0
END IF
BXX_RHS(I,J)=BXX(I,J)*(U(I+1,J)-U(I-1,J))/DX+(BYX(I,J)+BXY(I,J))*(U(I,J+1)-U(I,J-1))/(2.D0*DY)&
            -U(I,J)*(BXX(I+1,J)-BXX(I-1,J))/(2.D0*DX)-V(I,J)*(BXX(I,J+1)-BXX(I,J-1))/(2.D0*DY)
BYY_RHS(I,J)=BYY(I,J)*(V(I,J+1)-V(I,J-1))/DY+(BXY(I,J)+BYX(I,J))*(V(I+1,J)-V(I-1,J))/(2.D0*DX)&
            -U(I,J)*(BYY(I+1,J)-BYY(I-1,J))/(2.D0*DX)-V(I,J)*(BYY(I,J+1)-BYY(I,J-1))/(2.D0*DY)
BXY_RHS(I,J)=BXX(I,J)*(V(I+1,J)-V(I-1,J))/(2.D0*DX)+BYY(I,J)*(U(I,J+1)-U(I,J-1))/(2.D0*DY)&
            +BXY(I,J)*((U(I+1,J)-U(I-1,J))/(2.D0*DX)+(V(I,J+1)-V(I,J-1))/(2.D0*DY))&
            -U(I,J)*(BXY(I+1,J)-BXY(I-1,J))/(2.D0*DX)-V(I,J)*(BXY(I,J+1)-BXY(I,J-1))/(2.D0*DY)
BYX_RHS(I,J)=BXX(I,J)*(V(I+1,J)-V(I-1,J))/(2.D0*DX)+BYY(I,J)*(U(I,J+1)-U(I,J-1))/(2.D0*DY)&
            +BYX(I,J)*((U(I+1,J)-U(I-1,J))/(2.D0*DX)+(V(I,J+1)-V(I,J-1))/(2.D0*DY))&
            -U(I,J)*(BYX(I+1,J)-BYX(I-1,J))/(2.D0*DX)-V(I,J)*(BYX(I,J+1)-BYX(I,J-1))/(2.D0*DY)
END DO
END DO
BXX=BXX+DT*BXX_RHS
BXY=BXY+DT*BXY_RHS
BYX=BYX+DT*BYX_RHS
BYY=BYY+DT*BYY_RHS
!------------------------------------------------------------------------------
END SUBROUTINE LEFT_CAUCHY_GREEN_DEFORMATION
!==============================================================================
SUBROUTINE VOR_TRANSPORT_RHS (PHIS,U,V,OMEGA,CDHST,OMEGA_RHS)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: CDHST,U,V,PHIS,OMEGA,OMEGA_RHS
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
OMEGA_RHS(I,J)=NUF*(1.D0-PHIS(I,J))*((OMEGA(I+1,J)-2.D0*OMEGA(I,J)+OMEGA(I-1,J))/DX2&!
              +(OMEGA(I,J+1)-2.D0*OMEGA(I,J)+OMEGA(I,J-1))/DY2)&
			  +INV_RHO*CDHST(I,J)&
			  -0.5D0*U(I,J)*(OMEGA(I+1,J)-OMEGA(I-1,J))/DX&
			  -0.5D0*V(I,J)*(OMEGA(I,J+1)-OMEGA(I,J-1))/DY
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE VOR_TRANSPORT_RHS
!==============================================================================
SUBROUTINE LEVEL_SET_FUNCTION_LAXWENDROFF (U,V,PHIS)
USE MODULE1
IMPLICIT NONE
REAL(8),DIMENSION(IM,JM)     :: PHIS,U,V
INTEGER                      :: I,J
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
PHIS(I,J)=PHIS(I,J)-0.5D0*CX*U(I,J)*(PHIS(I+1,J)-PHIS(I-1,J))&
                  -0.5D0*CY*V(I,J)*(PHIS(I,J+1)-PHIS(I,J-1))&
                  +0.5D0*CX2*U(I,J)*U(I,J)*(PHIS(I+1,J)-2.D0*PHIS(I,J)+PHIS(I-1,J))&
                  +0.5D0*CY2*V(I,J)*V(I,J)*(PHIS(I,J+1)-2.D0*PHIS(I,J)+PHIS(I,J-1))
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE LEVEL_SET_FUNCTION_LAXWENDROFF
!==============================================================================
SUBROUTINE LEVEL_SET_FUNCTION_MCCORMACK (U,V,PHIS)
USE MODULE1
IMPLICIT NONE
REAL(8),DIMENSION(IM,JM)     :: PHIS,U,V,PHISS
INTEGER                      :: I,J
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
PHISS(I,J)=PHIS(I,J)-CX*U(I,J)*(PHIS(I+1,J)-PHIS(I,J))&
                  -CY*V(I,J)*(PHIS(I,J+1)-PHIS(I,J))
END DO
END DO
DO I=2,IM-1
DO J=2,JM-1
PHIS(I,J)=0.5D0*(PHIS(I,J)+PHISS(I,J)-CX*U(I,J)*(PHISS(I,J)-PHISS(I-1,J))&
                                  -CY*V(I,J)*(PHISS(I,J)-PHISS(I,J-1)))
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE LEVEL_SET_FUNCTION_MCCORMACK
!==============================================================================
SUBROUTINE DEVIATORIC_HYPERELASTIC_TENSOR (PHIS,BXX,BXY,BYX,BYY,SIGMA_SH_XX,SIGMA_SH_XY,&
                                           SIGMA_SH_YX,SIGMA_SH_YY)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: BXX,BXY,BYX,BYY,SIGMA_SH_XX,SIGMA_SH_XY,SIGMA_SH_YX,SIGMA_SH_YY
REAL(8),DIMENSION(IM,JM)     :: PHIS_AUX,PHIS
!------------------------------------------------------------------------------
DO I=1,IM
DO J=1,JM
IF (PHIS(I,J).LE.PHIS_MIN) THEN
PHIS_AUX(I,J)=0.D0
ELSE
PHIS_AUX(I,J)=PHIS(I,J)**0.5D0
END IF
SIGMA_SH_XX(I,J)=(2.D0/3.D0)*C1*PHIS_AUX(I,J)*(2.D0*BXX(I,J)-BYY(I,J))&
                -(2.D0/3.D0)*C2*(BXY(I,J)*BYX(I,J)-BXX(I,J)*BYY(I,J))&
                +(4.D0/3.D0)*C3*(2.D0*BXX(I,J)**2.D0-BYY(I,J)**2.D0)&
                +(4.D0/3.D0)*C3*BXX(I,J)*BYY(I,J)&
                -4.D0*C3*PHIS_AUX(I,J)*(2.D0*BXX(I,J)-BYY(I,J))
SIGMA_SH_YY(I,J)=(2.D0/3.D0)*C1*PHIS_AUX(I,J)*(2.D0*BYY(I,J)-BXX(I,J))&
                -(2.D0/3.D0)*C2*(BXY(I,J)*BYX(I,J)-BXX(I,J)*BYY(I,J))&
                +(4.D0/3.D0)*C3*(2.D0*BYY(I,J)**2.D0-BXX(I,J)**2.D0)&
                +(4.D0/3.D0)*C3*BXX(I,J)*BYY(I,J)&
                -4.D0*C3*PHIS_AUX(I,J)*(2.D0*BYY(I,J)-BXX(I,J))
SIGMA_SH_XY(I,J)=2.D0*C1*PHIS_AUX(I,J)*BXY(I,J)+4.D0*C3*BXY(I,J)*(BXX(I,J)+BYY(I,J)-3.D0*PHIS_AUX(I,J))
SIGMA_SH_YX(I,J)=2.D0*C1*PHIS_AUX(I,J)*BYX(I,J)+4.D0*C3*BYX(I,J)*(BXX(I,J)+BYY(I,J)-3.D0*PHIS_AUX(I,J))
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE DEVIATORIC_HYPERELASTIC_TENSOR
!==============================================================================
SUBROUTINE CURL_DIVERGENCE_HYPERELASTIC_STRESS_TENSOR (SIGMA_SH_XX,SIGMA_SH_XY,&
                                                       SIGMA_SH_YX,SIGMA_SH_YY,&
                                                       CDHST)
USE MODULE1
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: SIGMA_SH_XX,SIGMA_SH_XY,SIGMA_SH_YX,SIGMA_SH_YY,CDHST
!------------------------------------------------------------------------------
DO I=2,IM-1
DO J=2,JM-1
CDHST(I,J)=(SIGMA_SH_XY(I+1,J)-2.D0*SIGMA_SH_XY(I,J)+SIGMA_SH_XY(I-1,J))/DX2&
          -(SIGMA_SH_YX(I,J+1)-2.D0*SIGMA_SH_YX(I,J)+SIGMA_SH_YX(I,J-1))/DY2&
          +(SIGMA_SH_YY(I+1,J+1)-SIGMA_SH_YY(I-1,J+1)-SIGMA_SH_YY(I+1,J-1)+SIGMA_SH_YY(I-1,J-1))/FDXDY&
          -(SIGMA_SH_XX(I+1,J+1)-SIGMA_SH_XX(I-1,J+1)-SIGMA_SH_XX(I+1,J-1)+SIGMA_SH_XX(I-1,J-1))/FDXDY
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE CURL_DIVERGENCE_HYPERELASTIC_STRESS_TENSOR
!==============================================================================
SUBROUTINE MGLIN (U,F)
USE MODULE1,ONLY:IM,JM,DX,DY
USE MULTIGRID,ONLY:NCYCLE,NPRE,NPOST,NFINAL,NG_X,NG_Y,NG_D,NG_L,MEMLEN
IMPLICIT NONE
INTEGER J,JCYCLE,JJ,JPOST,JPRE,MEM,NF_X,NF_Y,NGRID,NN_X,NN_Y,IRES(NG_L),IRHO(NG_L),IRHS(NG_L),IU(NG_L),MALOC,I_1,I_2,I_3
REAL(8) Z,H,F(IM,JM),U(IM,JM),E(IM,JM),SUM_E
COMMON /MEMORY/ Z(MEMLEN),MEM
!------------------------------------------------------------------------------
MEM=0
DO I_1=NFINAL,NG_L
   I_2=2**I_1+1
   I_3=2**(I_1+NG_D)+1
   IF (I_1<NG_L) IRHO(I_1)=MALOC(I_3*I_2)
   IF (I_1>NFINAL) IRES(I_1)=MALOC(I_3*I_2)
   IRHS(I_1)=MALOC(I_3*I_2)
   IU(I_1)=MALOC(I_3*I_2)
ENDDO
!------------------------------------------------------------------------------
NN_X=IM/2+1
NN_Y=JM/2+1
NGRID=NG_L-1
!------------------------------------------------------------------------------
CALL RSTRCT(Z(IRHO(NGRID)),F,NN_X,NN_Y)
!------------------------------------------------------------------------------
DO I_1=NG_L-2,NFINAL,-1
   NN_X=NN_X/2+1
   NN_Y=NN_Y/2+1
   NGRID=NGRID-1
   CALL RSTRCT(Z(IRHO(NGRID)),Z(IRHO(NGRID+1)),NN_X,NN_Y)
ENDDO
!------------------------------------------------------------------------------
NN_X=2**(NG_X-NG_L+NFINAL)+1
NN_Y=2**(NG_Y-NG_L+NFINAL)+1
!------------------------------------------------------------------------------
CALL SLVSML(Z(IU(NFINAL)),Z(IRHO(NFINAL)),U,1)
!------------------------------------------------------------------------------
NGRID=NG_L
DO J=NFINAL+1,NGRID
   NN_X=2*NN_X-1
   NN_Y=2*NN_Y-1
!------------------------------------------------------------------------------
   CALL SETBND(U,Z(IU(J)),NN_X-1,NN_Y-1)
!------------------------------------------------------------------------------
   CALL INTERP(Z(IU(J)),Z(IU(J-1)),NN_X,NN_Y)
!------------------------------------------------------------------------------
   IF (J /= NGRID) THEN
      CALL COPY(Z(IRHS(J)),Z(IRHO(J)),NN_X,NN_Y)
   ELSE
      CALL COPY(Z(IRHS(J)),F,NN_X,NN_Y)
   ENDIF
!------------------------------------------------------------------------------
   DO JCYCLE=1,NCYCLE
      NF_X=NN_X
      NF_Y=NN_Y
!------------------------------------------------------------------------------
	  DO JJ=J,NFINAL+1,-1
!------------------------------------------------------------------------------
         DO JPRE=1,NPRE
            CALL RELAX(Z(IU(JJ)),Z(IRHS(JJ)),NF_X,NF_Y)
         ENDDO
!------------------------------------------------------------------------------
         CALL RESID(Z(IRES(JJ)),Z(IU(JJ)),Z(IRHS(JJ)),NF_X,NF_Y)
!------------------------------------------------------------------------------
         NF_X=NF_X/2+1
         NF_Y=NF_Y/2+1
         CALL RSTRCT(Z(IRHS(JJ-1)),Z(IRES(JJ)),NF_X,NF_Y)
!------------------------------------------------------------------------------
         CALL FILL0(Z(IU(JJ-1)),NF_X,NF_Y)
!------------------------------------------------------------------------------
      ENDDO
!------------------------------------------------------------------------------
	  CALL SLVSML(Z(IU(NFINAL)),Z(IRHS(NFINAL)),U,0)
!------------------------------------------------------------------------------
      NF_X=2**(NG_X-NG_L+NFINAL)+1
      NF_Y=2**(NG_Y-NG_L+NFINAL)+1
!------------------------------------------------------------------------------
      DO JJ=NFINAL+1,J
         NF_X=2*NF_X-1
         NF_Y=2*NF_Y-1
!------------------------------------------------------------------------------
         CALL ADDINT(U,Z(IU(JJ)),Z(IU(JJ-1)),Z(IRES(JJ)),NF_X,NF_Y)
!------------------------------------------------------------------------------
         DO JPOST=1,NPOST
            CALL RELAX(Z(IU(JJ)),Z(IRHS(JJ)),NF_X,NF_Y)
         ENDDO
!------------------------------------------------------------------------------
      ENDDO
!------------------------------------------------------------------------------
   ENDDO
!------------------------------------------------------------------------------
ENDDO
!------------------------------------------------------------------------------
CALL COPY(U,Z(IU(NGRID)),IM,JM)
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
FUNCTION MALOC(LEN)
!USE MEM_MAN,ONLY :MEMLEN
USE MULTIGRID,ONLY :MEMLEN
IMPLICIT NONE
INTEGER MALOC,LEN
INTEGER MEM
REAL(8) Z
COMMON /memory/ Z(MEMLEN),MEM
!------------------------------------------------------------------------------
IF(MEM+LEN+1 > MEMLEN) PAUSE 'INSUFFICIENT MEMORY IN MALOC'
!------------------------------------------------------------------------------
Z(MEM+1)=LEN
MALOC=MEM+2
MEM=MEM+LEN+1
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE RSTRCT(UC,UF,NC_X,NC_Y)
IMPLICIT NONE
INTEGER NC_X,NC_Y
REAL(8) UC(NC_X,NC_Y),UF(2*NC_X-1,2*NC_Y-1)
INTEGER IC,IF,JC,JF
!------------------------------------------------------------------------------
DO JC=2,NC_Y-1
   JF=2*JC-1
   DO IC=2,NC_X-1
      IF=2*IC-1
      UC(IC,JC)=0.5*UF(IF,JF)+0.125*(UF(IF+1,JF)+UF(IF-1,JF)+UF(IF,JF+1)+UF(IF,JF-1))
   ENDDO
ENDDO
!------------------------------------------------------------------------------
RETURN
DO IC=1,NC_X
  IF=2*IC-1
  UC(IC,1)=UF(IF,1)
  UC(IC,NC_Y)=UF(IF,2*NC_Y-1)
ENDDO
DO JC=1,NC_Y
  JF=2*JC-1
  UC(1,JC)=UF(1,JF)
  UC(NC_X,JC)=UF(2*NC_X-1,JF)
ENDDO
!------------------------------------------------------------------------------
DO JC=2,NC_Y-1
   JF=2*JC-1
   DO IC=2,NC_X-1
      IF=2*IC-1
      UC(IC,JC)=0.5*UF(IF,JF)+0.125*(UF(IF+1,JF)+UF(IF-1,JF)+UF(IF,JF+1)+UF(IF,JF-1))
   ENDDO
ENDDO
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE SLVSML(U,RHS,UBND,TEST)
USE MULTIGRID,ONLY:NFINAL,NG_L,NG_X,NG_Y
USE MODULE1,ONLY:IM,JM
IMPLICIT NONE
REAL(8) RHS(2**(NG_X-NG_L+NFINAL)+1,2**(NG_Y-NG_L+NFINAL)+1),U(2**(NG_X-NG_L+NFINAL)+1,2**(NG_Y-NG_L+NFINAL)+1),UBND(IM,JM),H
INTEGER TEST,N_X,N_Y
!------------------------------------------------------------------------------
N_X=2**(NG_X-NG_L+NFINAL)
N_Y=2**(NG_Y-NG_L+NFINAL)
!------------------------------------------------------------------------------
CALL FILL0(U,N_X+1,N_Y+1)
!------------------------------------------------------------------------------
IF (TEST==1) CALL SETPRE(UBND,U,N_X,N_Y)
!------------------------------------------------------------------------------
CALL RELAX(U,RHS,N_X+1,N_Y+1)
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE SETPRE(UF,UC,NPOINT_X,NPOINT_Y)
USE MODULE1,ONLY:IM,JM
IMPLICIT NONE
INTEGER I,J,IF,JF,NPOINT_X,NPOINT_Y
REAL(8) UC(NPOINT_X+1,NPOINT_Y+1),UF(IM,JM)
!------------------------------------------------------------------------------
DO J=1,NPOINT_Y+1
   JF=((J-1)*(JM-1)/(NPOINT_Y))+1
   DO I=1,NPOINT_X+1
      IF=((I-1)*(IM-1)/(NPOINT_X))+1
      UC(I,J)=UF(IF,JF)
   ENDDO
ENDDO
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE SETBND(UF,UC,NPOINT_X,NPOINT_Y)
USE MODULE1,ONLY:IM,JM
IMPLICIT NONE
INTEGER I,J,IF,JF,NPOINT_X,NPOINT_Y
REAL(8) UC(NPOINT_X+1,NPOINT_Y+1),UF(IM,JM)
!------------------------------------------------------------------------------ I ROWS SET BOUNDARY
JF=1
J=1
DO I=1,NPOINT_X+1
   IF=((I-1)*(IM-1)/(NPOINT_X))+1
   UC(I,J)=UF(IF,JF)
ENDDO
!------------------------------------------------------------------------------
JF=JM
J=NPOINT_Y+1
DO I=1,NPOINT_X+1
   IF=((I-1)*(IM-1)/(NPOINT_X))+1
   UC(I,J)=UF(IF,JF)
ENDDO
!------------------------------------------------------------------------------ J ROWS SET BOUNDARY
IF=1
I=1
DO J=2,NPOINT_Y
   JF=((J-1)*(JM-1)/(NPOINT_Y))+1
   UC(I,J)=UF(IF,JF)
ENDDO
!------------------------------------------------------------------------------
IF=IM
I=NPOINT_X+1
DO J=2,NPOINT_Y
   JF=((J-1)*(JM-1)/(NPOINT_Y))+1
   UC(I,J)=UF(IF,JF)
ENDDO
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE INTERP(UF,UC,NF_X,NF_Y)
USE MODULE1,ONLY:IM,JM
IMPLICIT NONE
INTEGER NF_X,NF_Y
REAL(8) UC(NF_X/2+1,NF_Y/2+1),UF(NF_X,NF_Y)
INTEGER IC,IF,JC,JF,NC_X,NC_Y
NC_X=NF_X/2+1
NC_Y=NF_Y/2+1
!------------------------------------------------------------------------------ SET SIMILAR POINT FROM COARSE GRID TO FINE GRID
DO JC=2,NC_Y-1
   JF=2*JC-1
   DO IC=2,NC_X-1
      UF(2*IC-1,JF)=UC(IC,JC)
   ENDDO
ENDDO
!------------------------------------------------------------------------------ SET NEW VALUE FOR J ROW
DO IF=3,NF_X-1,2
   DO JF=2,NF_Y,2
      UF(IF,JF)=(UF(IF,JF-1)+UF(IF,JF+1))/2
   ENDDO
ENDDO
!------------------------------------------------------------------------------ SET NEW VALUE FOR I ROW
DO JF=3,NF_Y-1,2
   DO IF=2,NF_X,2
      UF(IF,JF)=(UF(IF+1,JF)+UF(IF-1,JF))/2
   ENDDO
ENDDO
!------------------------------------------------------------------------------ SET MID POINTS
DO JF=2,NF_Y,2
   DO IF=2,NF_X,2
      UF(IF,JF)=(UF(IF+1,JF+1)+UF(IF+1,JF-1)+UF(IF-1,JF+1)+UF(IF-1,JF-1))/4
   ENDDO
ENDDO
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE RELAX(U,RHS,N_X,N_Y)
USE MODULE1,ONLY:LX,LY
IMPLICIT NONE
INTEGER N_X,N_Y
REAL(8) RHS(N_X,N_Y),U(N_X,N_Y)
INTEGER I,ISW,J
REAL(8) H2
H2=LX*LX/(N_X-1)/(N_X-1)
!------------------------------------------------------------------------------
ISW=1
!------------------------------------------------------------------------------
DO J=2,N_Y-1
   DO I=ISW+1,N_X-1,2
      U(I,J)=0.25D0*(U(I+1,J)+U(I-1,J)+U(I,J+1)+U(I,J-1)-H2*RHS(I,J))
   ENDDO
   ISW=3-ISW
ENDDO
!------------------------------------------------------------------------------
ISW=2
!------------------------------------------------------------------------------
DO J=2,N_Y-1
   DO I=ISW+1,N_X-1,2
      U(I,J)=0.25D0*(U(I+1,J)+U(I-1,J)+U(I,J+1)+U(I,J-1)-H2*RHS(I,J))
   ENDDO
   ISW=3-ISW
ENDDO
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE RESID(RES,U,RHS,N_X,N_Y)
USE MODULE1,ONLY:LX,LY
IMPLICIT NONE
INTEGER N_X,N_Y
REAL(8) RES(N_X,N_Y),RHS(N_X,N_Y),U(N_X,N_Y)
INTEGER I,J
REAL(8) H2I
H2I=(N_X-1)*(N_X-1)/LX/LX
!------------------------------------------------------------------------------
DO J=2,N_Y-1
   DO I=2,N_X-1
      RES(I,J)=-H2I*(U(I+1,J)+U(I-1,J)+U(I,J+1)+U(I,J-1)-4.D0*U(I,J))+RHS(I,J)
   ENDDO
ENDDO
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE ADDINT(U,UF,UC,RES,NF_X,NF_Y)
USE MODULE1,ONLY:IM,JM
IMPLICIT NONE
INTEGER NF_X,NF_Y,I,J
REAL(8) RES(NF_X,NF_Y),UC(NF_X/2+1,NF_Y/2+1),UF(NF_X,NF_Y),U(IM,JM)
!------------------------------------------------------------------------------
CALL INTERP(RES,UC,NF_X,NF_Y)
!------------------------------------------------------------------------------
!UF=UF+RES
DO J=2,NF_Y-1
   DO I=2,NF_X-1
      UF(I,J)=UF(I,J)+RES(I,J)
   ENDDO
ENDDO
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE COPY(AOUT,AIN,N_X,N_Y)
IMPLICIT NONE
INTEGER N_X,N_Y
REAL(8) AIN(N_X,N_Y),AOUT(N_X,N_Y)
!------------------------------------------------------------------------------
AOUT=AIN
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================
SUBROUTINE FILL0(U,N_X,N_Y)
IMPLICIT NONE
INTEGER N_X,N_Y
REAL(8) U(N_X,N_Y)
!------------------------------------------------------------------------------
U=0.D0
!------------------------------------------------------------------------------
RETURN
END
!==============================================================================

