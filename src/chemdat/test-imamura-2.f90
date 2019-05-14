program main
  ! 
  ! Imamura and Hashimoto (1998) に基づいた金星雲微物理計算モジュールのテスト
  ! 
  
  use Imamura1998
    
  !暗黙の型宣言禁止
  implicit none
  
  !変数の宣言
  integer            :: i, j
  integer, parameter :: Num = 6850
  real(8), parameter :: dx = 1.0d-1
  real(8), parameter :: Boltz = 1.38064852d-23  !!Boltzmann constant  
  real(8)            :: alt, temp, press, n1, n2, p1, p2, x1
  real(8)            :: molfr
  real(8)            :: n1_gas, n2_gas, n1_liq, n2_liq
  integer, parameter :: sw = 2  !ニュートン法
  logical            :: flag
  
  write(*,*) "===================================="
  write(*,*) "  Bisection vs. Newton     "
  write(*,*) "===================================="

  open (17, file='fukuhara_data.csv', status='old')
  read (17, '()')       ! ヘッダ行の読み飛ばし

  do i = 1, Num     

     !データ読み込み
     read (17, *) alt, press, temp, n1, n2, p1, p2, x1

!     !データを間引く
!     if ( mod(i, 10) .ne. 0) cycle

     !平衡状態の計算
     call Imamura1998_EquivState( Temp, n1, n2, n1_gas, n2_gas, n1_liq, n2_liq, molfr, sw, flag )

     !出力
     write(*,*) "*",i,",",real(temp),",",real(n1_gas),",",real(n1_liq),",",real(n2_gas),",",real(n2_liq),",",real(molfr),",",real(x1),","
     
  end do

  close (17)
  
end program main
