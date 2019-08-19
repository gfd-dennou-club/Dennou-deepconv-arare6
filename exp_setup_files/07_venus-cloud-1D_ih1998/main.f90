program main
  !
  ! Imamura and Hashimoto (1998) に準じた鉛直 1 次元計算
  ! 途中, write 文で deepconv に読み込むためのファイルを標準出力に出力している. 
  ! 

  !モジュールの呼び出し
  use gtool_history
  use Cloudphys_IH1998
  
  !暗黙の型宣言禁止
  implicit none
  
  !変数の宣言
  integer            :: i, j, k, n
  integer, parameter :: Num = 365
!  integer, parameter :: Num = 8000
  real(8), parameter :: Boltz = 1.38064852d-23  !!Boltzmann constant  
  real(8)            :: alt(Num), temp(Num), press(Num)
  real(8), parameter :: molfr1_sfc = 10.0d-6   !大気下端でのモル比. H2SO4 : 10 ppm, H2O : 30 ppm.
  real(8), parameter :: molfr2_sfc = 30.0d-6   !大気下端でのモル比. H2SO4 : 10 ppm, H2O : 30 ppm.
  real(8)            :: molfr1(Num), molfr2(Num), x1(Num)
  real(8)            :: n1_sat(Num), n2_sat(Num)
  real(8)            :: n1_gas(Num), n2_gas(Num), n1_liq(Num), n2_liq(Num), n1(Num), n2(Num)
  real(8)            :: n1_a(Num), n2_a(Num)
  real(8)            :: Wsed(Num)
  real(8)            :: Dn1Dt_prdt(Num), Dn1Dt_loss(Num), Dn1Dt_fall(Num), Dn1Dt(Num)
  real(8)            :: Dn2Dt_prdt(Num), Dn2Dt_loss(Num), Dn2Dt_fall(Num), Dn2Dt(Num)
  integer, parameter :: sw = 2  !ニュートン法
  logical            :: flag
  real(8)            :: tinit = 0
  integer            :: ndisp = 100
  real(8)            :: dt = 2000.0d0
  real(8)            :: DelZ
  integer            :: CalNum = 100000
  real(8)            :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
  real(8)            :: p1, p2
  
  Dn1Dt_prdt = 0.0d0
  Dn1Dt_loss = 0.0d0
  Dn1Dt_fall = 0.0d0
  Dn1Dt      = 0.0d0
  Dn2Dt_prdt = 0.0d0
  Dn2Dt_loss = 0.0d0
  Dn2Dt_fall = 0.0d0
  Dn2Dt      = 0.0d0
 
  !csv からデータ読み込み.
  ! venus_andoh_initial.dat は
!  open (17, file='test-imamura-data_from_fukuhara2.csv', status='old')
  open (17, file='venus_andoh_initial.dat', status='old')
!  read (17, '()')       ! ヘッダ行の読み飛ばし
  
  !初期化
  molfr1(1) = molfr1_sfc
  molfr2(1) = molfr2_sfc
  
  !データ読み込み. 余計なものもあるが気にしない.   
  !高度・温度・圧力は外部から与えるためにデータ読み込みが必要. 
  do i = 1, Num
    read (17, *) alt(i), press(i), temp(i), tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
!    read (17, *) alt(i), temp(i), press(i)
!    write(*,*) alt(i), temp(i), press(i)
  end do

  !高度刻み幅
!  alt  = alt * 1.0d3
  DelZ = alt(2) - alt(1)
  
  n1 = 0.0d0
  n2 = 0.0d0
  n1_a = 0.0d0
  n2_a = 0.0d0

  call HistoryCreate( &                        ! ヒストリー作成
    & file='imamura1998-20190819.nc', title='Imamura and Hashimoto 1998', &
    & source='deepconv/arare6',                          &
    & institution='GFD_Dennou Club deepconv project',    &
    & dims=(/'z','t'/), dimsizes=(/Num,0/),               &
    & longnames=(/'Z-coordinate','time        '/),       &
    & units=(/'m','s'/),                                 &
    & origin=real(tinit), interval=real(ndisp*dt) )

  call HistoryPut('z',alt)                     ! 次元変数出力

  call HistoryAddVariable( &                   ! 変数定義
    & varname='Temp', dims=(/'z'/), &
    & longname='Temperature', units='K', xtype='float')
  
  call HistoryAddVariable( &                   ! 変数定義
    & varname='Press', dims=(/'z'/), &
    & longname='Pressure', units='Pa', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='n1', dims=(/'z','t'/), &
    & longname='number density (H2SO4)', units='1', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='n2', dims=(/'z','t'/), &
    & longname='number density (H2O)', units='1', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='n1_sat', dims=(/'z','t'/), &
    & longname='saturated number density (H2SO4)', units='1', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='n2_sat', dims=(/'z','t'/), &
    & longname='saturated number density (H2O)', units='1', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='n1_gas', dims=(/'z','t'/), &
    & longname='number density (H2SO4 gas)', units='1', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='n2_gas', dims=(/'z','t'/), &
    & longname='number density (H2O gas)', units='1', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='n1_liq', dims=(/'z','t'/), &
    & longname='number density (H2SO4 liq)', units='1', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='n2_liq', dims=(/'z','t'/), &
    & longname='number density (H2O liq)', units='1', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='x1', dims=(/'z','t'/), &
    & longname='concentration of H2SO4 solution', units='1', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='molfr1', dims=(/'z','t'/), &
    & longname='mole fraction (H2SO4)', units='1', xtype='float')  

  call HistoryAddVariable( &                   ! 変数定義
    & varname='molfr2', dims=(/'z','t'/), &
    & longname='mole fraction (H2O)', units='1', xtype='float')  

  call HistoryAddVariable( &                   ! 変数定義
    & varname='H2SO4_prdt', dims=(/'z','t'/), &
    & longname='H2SO4 production rate', units='1/s', xtype='float')  
  
  call HistoryAddVariable( &                   ! 変数定義
    & varname='H2SO4_loss', dims=(/'z','t'/), &
    & longname='H2SO4 loss rate', units='1/s', xtype='float')  

  call HistoryAddVariable( &                   ! 変数定義
    & varname='H2SO4_fall', dims=(/'z','t'/), &
    & longname='H2SO4 sedimentation rate', units='1/s', xtype='float')    

  call HistoryAddVariable( &                   ! 変数定義
    & varname='H2SO4_tend', dims=(/'z','t'/), &
    & longname='H2SO4 tendency', units='1/s', xtype='float')  
  
  call HistoryAddVariable( &                   ! 変数定義
    & varname='H2O_prdt', dims=(/'z','t'/), &
    & longname='H2O production rate', units='1/s', xtype='float')  
  
  call HistoryAddVariable( &                   ! 変数定義
    & varname='H2O_loss', dims=(/'z','t'/), &
    & longname='H2O loss rate', units='1/s', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='H2O_fall', dims=(/'z','t'/), &
    & longname='H2O sedimentation rate', units='1/s', xtype='float')

  call HistoryAddVariable( &                   ! 変数定義
    & varname='H2O_tend', dims=(/'z','t'/), &
    & longname='H2O tendency', units='1/s', xtype='float')    

  call HistoryAddVariable( &                   ! 変数定義
    & varname='Wsed', dims=(/'z','t'/), &
    & longname='Wsed', units='m/s', xtype='float')    

  call HistoryPut('Temp',  temp)               ! 変数出力
  call HistoryPut('Press', press)              ! 変数出力  
  
!!!1
!!! 各高度での計算
!!!
  do i = 2, Num
    
    !数密度の計算.
    !モル比が保存するとして, 1 つ前の高度でのモル比を用いる. 
    n1(i) = press(i) * molfr1(i-1) / ( Boltz * temp(i) )
    n2(i) = press(i) * molfr2(i-1) / ( Boltz * temp(i) )
    
    !平衡状態の計算
    call IH1998_EquivState( Temp(i), n1(i), n2(i), n1_gas(i), n2_gas(i), n1_liq(i), n2_liq(i), x1(i), 1, flag )
    
    !気体のモル比
    molfr1(i) = n1_gas(i) * Boltz * temp(i) / press(i)
    molfr2(i) = n2_gas(i) * Boltz * temp(i) / press(i)

    !飽和蒸気圧
    call IH1998_SatPressRef( Temp(i), p1, p2 )
    n1_sat(i) = p1 / ( Boltz * temp(i) )
    n2_sat(i) = p2 / ( Boltz * temp(i) )
   
    !出力
!    write(*,*) "*",i,",",real(temp(i)),",",real(n1(i)),",",real(n1_gas(i)),",",real(n1_liq(i)),",",real(n2(i)),",",real(n2_gas(i)),",",real(n2_liq(i)),",",real(x1(i)),",",real(molfr1(i)),",",real(molfr2(i)),","
  end do
  n1_sat(1) = n1_sat(2)
  n2_sat(1) = n2_sat(2)

  call HistoryPut('n1',  max( n1, 1e10 ))
  call HistoryPut('n2',  max( n2, 1e10 ))
  call HistoryPut('n1_gas', max( n1_gas , 1e10 ))
  call HistoryPut('n2_gas', max( n2_gas , 1e10 ))
  call HistoryPut('n1_liq', max( n1_liq , 1e10 ))
  call HistoryPut('n2_liq', max( n2_liq , 1e10 ))
  call HistoryPut('n1_sat', max( n1_sat , 1e10 ))
  call HistoryPut('n2_sat', max( n2_sat , 1e10 ))
  call HistoryPut('x1',     x1 + 1e-20)
  call HistoryPut('molfr1', molfr1 + 1e2)
  call HistoryPut('molfr2', molfr2 + 1e2)
  call HistoryPut('H2SO4_prdt', Dn1Dt_prdt)
  call HistoryPut('H2SO4_loss', Dn1Dt_loss)
  call HistoryPut('H2SO4_fall', Dn1Dt_fall)
  call HistoryPut('H2SO4_tend', Dn1Dt)
  call HistoryPut('H2O_prdt',   Dn2Dt_prdt)
  call HistoryPut('H2O_loss',   Dn2Dt_loss)        
  call HistoryPut('H2O_fall',   Dn2Dt_fall)  
  call HistoryPut('H2O_tend',   Dn2Dt)
  call HistoryPut('Wsed',       Wsed)


  do i = 2, Num
    
    !H2SO4 モル比の調整
    !H2SO4 の分解は高度 39 km 以下で生じる (Jenkins et al., 1994). tanh で H2SO4 のモル比を落とす. 
    if ( molfr1(i) > molfr1_sfc * ( 0.5d0 * tanh( ( alt(i) - 37.5d3 ) * 1.0d-3 ) + 0.5d0 ) ) then
      molfr1(i) = molfr1_sfc * ( 0.5d0 * tanh( ( alt(i) - 40.0d3 ) * 1.0d-3 ) + 0.5d0 )
      n1_gas(i) = press(i) * molfr1(i) / ( Boltz * temp(i) )
      n1(i)     = n1_gas(i)
    end if
!   if ( mod(i, 100) == 0 ) then
!      write(*,*) i, alt(i)
!      write(*,*) molfr1(i), molfr1_sfc * ( 0.5d0 * tanh( ( alt(i) - 40.0d3 ) * 1.0d-3 ) + 0.5d0 ) 
!    end if

    !出力
    !    write(*,*) "*",i,",",real(temp(i)),",",real(n1(i)),",",real(n1_gas(i)),",",real(n1_liq(i)),",",real(n2(i)),,",real(n2_gas(i)),",",real(n2_liq(i)),",",real(x1(i)),","
    write(*,*) "*",real(alt(i)),",",real(n1(i)),",",real(n1_gas(i)),",",real(n2(i)),",",real(n2_gas(i)),","
!    write(*,*) Temp(i), molfr1(i), molfr2(i),  molfr1_sat(i), molfr2_sat(i)
  end do

  !test
  n1 = n1_gas
  n2 = n2_gas
  n1_liq = 0.0d0
  n2_liq = 0.0d0
  
  do i = 2, Num
    !落下項の計算. 高度のみの関数なので, 最初に一度行えば良い. 
    call IH1998_Sediment( alt(i), Wsed(i) )
    
    !H2SO4 の生成項・H2O の消滅項. 高度のみの関数なので, 最初に一度行えば良い. 
    call IH1998_H2SO4Prdt( alt(i), DelZ, Dn1Dt_prdt(i), Dn2Dt_loss(i) )
  end do

  n1_sat(1) = n1_sat(2)
  n2_sat(1) = n2_sat(2)
  n1_gas(1) = n1_gas(2)
  n2_gas(1) = n2_gas(2)
    
  call HistoryPut('n1',  max( n1, 1e10 ))
  call HistoryPut('n2',  max( n2, 1e10 ))
  call HistoryPut('n1_gas', max( n1_gas , 1e10 ))
  call HistoryPut('n2_gas', max( n2_gas , 1e10 ))
  call HistoryPut('n1_liq', max( n1_liq , 1e10 ))
  call HistoryPut('n2_liq', max( n2_liq , 1e10 ))
  call HistoryPut('n1_sat', max( n1_sat , 1e10 ))
  call HistoryPut('n2_sat', max( n2_sat , 1e10 ))
  call HistoryPut('x1',     x1 + 1e-20)
  call HistoryPut('molfr1', molfr1 + 1e2)
  call HistoryPut('molfr2', molfr2 + 1e2)
  call HistoryPut('H2SO4_prdt', Dn1Dt_prdt)
  call HistoryPut('H2SO4_loss', Dn1Dt_loss)
  call HistoryPut('H2SO4_fall', Dn1Dt_fall)
  call HistoryPut('H2SO4_tend', Dn1Dt)
  call HistoryPut('H2O_prdt',   Dn2Dt_prdt)
  call HistoryPut('H2O_loss',   Dn2Dt_loss)        
  call HistoryPut('H2O_fall',   Dn2Dt_fall)  
  call HistoryPut('H2O_tend',   Dn2Dt)
  call HistoryPut('Wsed',       Wsed)

  
!!!
!!! 時間発展
!!!  
  
  !時間発展のループ
  TimeLoop: do n = 1, CalNum
    
    
    HeightLoop: do i = 2, Num-1
      
      !H2SO4 の消滅項・H2O の生成項. 現在の凝結成分の存在度が必要. 
      call IH1998_H2SO4Loss( Temp(i), press(i), n1_gas(i), n2_gas(i), Dn1Dt_loss(i), Dn2Dt_prdt(i) )
      
      !落下 
      Dn1Dt_fall(i) = - ( Wsed(i+1) * n1_liq(i+1) - Wsed(i) * n1_liq(i) ) / (alt(i) - alt(i-1))
      Dn2Dt_fall(i) = - ( Wsed(i+1) * n2_liq(i+1) - Wsed(i) * n2_liq(i) ) / (alt(i) - alt(i-1))
      
    end do HeightLoop
    
    !tendency
    Dn1Dt = Dn1Dt_fall + Dn1Dt_prdt + Dn1Dt_loss
    Dn2Dt = Dn2Dt_fall + Dn2Dt_prdt + Dn2Dt_loss

!    Dn1Dt = Dn1Dt_prdt + Dn1Dt_loss
!    Dn2Dt = Dn2Dt_prdt + Dn2Dt_loss    
    
!    Dn1Dt = Dn1Dt_fall
!    Dn2Dt = Dn2Dt_fall
    
    n1_a = max( n1 + Dn1Dt * dt, 0.0d0)
    n2_a = max( n2 + Dn2Dt * dt, 0.0d0)


    HeightLoop2: do i = 2, Num

      !平衡状態の計算
      call IH1998_EquivState( Temp(i), n1_a(i), n2_a(i), n1_gas(i), n2_gas(i), n1_liq(i), n2_liq(i), x1(i), sw, flag )

      molfr1(i) = n1_gas(i) * Boltz * temp(i) / press(i)
      molfr2(i) = n2_gas(i) * Boltz * temp(i) / press(i)

      !飽和蒸気圧
      call IH1998_SatPressRef( Temp(i), p1, p2 )
      n1_sat(i) = p1 / ( Boltz * temp(i) )
      n2_sat(i) = p2 / ( Boltz * temp(i) )
    
    end do HeightLoop2

    n1_sat(1) = n1_sat(2)
    n2_sat(1) = n2_sat(2)
    n1_gas(1) = n1_gas(2)
    n2_gas(1) = n2_gas(2)

    
    write(*,*) "----  t = ", n * dt, " ----"
    if ( mod(n, ndisp) == 0 ) then
      write(*,*) "++++  t = ", n * dt, " ++++"

      call HistoryPut('n1',  max( n1, 1e10 ))
      call HistoryPut('n2',  max( n2, 1e10 ))
      call HistoryPut('n1_gas', max( n1_gas , 1e10 ))
      call HistoryPut('n2_gas', max( n2_gas , 1e10 ))
      call HistoryPut('n1_liq', max( n1_liq , 1e10 ))
      call HistoryPut('n2_liq', max( n2_liq , 1e10 ))
      call HistoryPut('n1_sat', max( n1_sat , 1e10 ))
      call HistoryPut('n2_sat', max( n2_sat , 1e10 ))
      call HistoryPut('x1',     x1 + 1e-20)
      call HistoryPut('molfr1', molfr1 + 1e2)
      call HistoryPut('molfr2', molfr2 + 1e2)
      call HistoryPut('H2SO4_prdt', Dn1Dt_prdt)
      call HistoryPut('H2SO4_loss', Dn1Dt_loss)
      call HistoryPut('H2SO4_fall', Dn1Dt_fall)
      call HistoryPut('H2SO4_tend', Dn1Dt)
      call HistoryPut('H2O_prdt',   Dn2Dt_prdt)
      call HistoryPut('H2O_loss',   Dn2Dt_loss)        
      call HistoryPut('H2O_fall',   Dn2Dt_fall)  
      call HistoryPut('H2O_tend',   Dn2Dt)
      call HistoryPut('Wsed',       Wsed)

    end if

    n1 = n1_a
    n2 = n2_a
    
  end do TimeLoop
  
  close (17)
  call HistoryClose
  
end program main

     
