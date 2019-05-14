program main
  ! 
  
  use Imamura1998
    
  !暗黙の型宣言禁止
  implicit none
  
  !変数の宣言
  integer            :: i
  integer, parameter :: Num = 6850
  real(8)            :: alt(Num), temp(Num), press(Num), tmp1(Num), tmp2(Num), p1(Num), p2(Num), tmp3(Num)
  real(8)            :: VelZ(Num)
  real(8)            :: Height

  !csv からデータ読み込み
  open (17, file='test-imamura-data_from_fukuhara.csv', status='old')
  read (17, '()')       ! ヘッダ行の読み飛ばし
  
  !データ読み込み. 余計なものもあるが気にしない.   
  !温度・圧力は外部から与えるためにデータ読み込みが必要. 
  do i = 1, Num
     read (17, *) alt(i), press(i), temp(i), tmp1(i), tmp2(i), p1(i), p2(i), tmp3(i)
  end do

  alt = alt * 1.0d3
  
!!!
!!! 各高度での計算
!!!
  do i = 1, Num

     if ( mod(i, 100) /= 0.0d0 ) cycle
     
     !落下速度
     call Imamura1998_Sediment( alt(i), VelZ(i) )
     
     !出力
     write(*,*) "*",i,",",real(alt(i)),",",real(velz(i)),","
     
  end do
  
  close (17)
  
end program main

     
