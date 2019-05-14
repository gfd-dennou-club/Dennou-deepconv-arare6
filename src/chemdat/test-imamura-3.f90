program main
  ! 
  
  use Imamura1998
    
  !暗黙の型宣言禁止
  implicit none
  
  !変数の宣言
  integer            :: i, j
  integer, parameter :: Num = 6850
  real(8), parameter :: dx = 1.0d-1
  real(8), parameter :: Boltz = 1.38064852d-23  !!Boltzmann constant  
  real(8)            :: alt(Num), temp(Num), press(Num), tmp1(Num), tmp2(Num), p1(Num), p2(Num), tmp3(Num)
  real(8), parameter :: molfr1_sfc = 10.0d-6   !大気下端でのモル比. H2SO4 : 10 ppm, H2O : 30 ppm.
  real(8), parameter :: molfr2_sfc = 30.0d-6   !大気下端でのモル比. H2SO4 : 10 ppm, H2O : 30 ppm.
  real(8)            :: molfr1(Num), molfr2(Num), x1(Num)
  real(8)            :: n1_gas(Num), n2_gas(Num), n1_liq(Num), n2_liq(Num), n1(Num), n2(Num)
  integer, parameter :: sw = 2  !ニュートン法
  logical            :: flag

  !csv からデータ読み込み
  open (17, file='test-imamura-data_from_fukuhara.csv', status='old')
  read (17, '()')       ! ヘッダ行の読み飛ばし
  
  !初期化
  molfr1(1) = molfr1_sfc
  molfr2(1) = molfr2_sfc

  !データ読み込み. 余計なものもあるが気にしない.   
  !温度・圧力は外部から与えるためにデータ読み込みが必要. 
  do i = 1, Num
     read (17, *) alt(i), press(i), temp(i), tmp1(i), tmp2(i), p1(i), p2(i), tmp3(i)
  end do
  
!!!
!!! 各高度での計算
!!!
  do i = 2, Num

     !数密度の計算.
     !モル比が保存するとして, 1 つ前の高度でのモル比を用いる. 
     n1(i) = press(i) * molfr1(i-1) / ( Boltz * temp(i) )
     n2(i) = press(i) * molfr2(i-1) / ( Boltz * temp(i) )
     
     !平衡状態の計算
     call Imamura1998_EquivState( Temp(i), n1(i), n2(i), n1_gas(i), n2_gas(i), n1_liq(i), n2_liq(i), x1(i), sw, flag )
     
     if ( x1(i) == -1.0d0 ) then 
        !凝結が生じる (flag == .true.) の場合は出力されたモル比を使う
        !
        molfr1(i) = n1_gas(i) * Boltz * temp(i) / press(i)
        molfr2(i) = n2_gas(i) * Boltz * temp(i) / press(i)
        
     else
        !凝結が生じない (flag == .false.) の場合はモル比は保存.
        !
        molfr1(i) = molfr1(i-1)
        molfr2(i) = molfr2(i-1)
     end if
     
     !出力
     write(*,*) "*",i,",",real(temp(i)),",",real(n1(i)),",",real(n1_gas(i)),",",real(n1_liq(i)),",",real(n2(i)),",",real(n2_gas(i)),",",real(n2_liq(i)),",",real(x1(i)),","
     
  end do

  close (17)
  
end program main

     
