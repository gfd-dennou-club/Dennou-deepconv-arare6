
module Cloudphys_IH1998

  != 金星雲微物計算モジュール
  !
  ! 元論文：Imamura and Hashimoto (1998), JGR103, E13, 31349. 
  !
  ! Authors::   SUGIYAMA Ko-ichiro, FUKUHARA Nozomu 
  !

  !暗黙の型宣言禁止
  implicit none

  !属性の指定
  private

  !H2SO4 水溶液のギブス自由エネルギーを求めるための係数.
  !Zeleznik (1991), J. Phys. Chem. Data の Table 6 に基づく.  
  real(8), private, save :: mu111
  real(8), private, save :: mu211
  real(8), private, save :: mu121
  real(8), private, save :: mu221
  real(8), private, save :: mu122
  real(8), private, save :: mu212

  real(8), private, save :: eps111
  real(8), private, save :: eps211
  real(8), private, save :: eps121
  real(8), private, save :: eps221
  real(8), private, save :: eps122
  real(8), private, save :: eps212


  !共通変数
  real(8), private, parameter :: Boltz = 1.38064852d-23  !!Boltzmann constant  
  real(8), private, parameter :: TempRef = 385.0d0
  real(8), private, parameter :: TempCr= 647.26d0

  !公開
  public IH1998_GibbsRDivRT
  public IH1998_DelChemPotRDivRT
  public IH1998_SatPress
  public IH1998_SatPressRef
  public IH1998_Newton
  public IH1998_Bisection
  public IH1998_EquivState
  public IH1998_Sediment
  public IH1998_H2SO4Prdt
  public IH1998_H2SO4Loss

  !関数を public に
  public Cloudphys_IH1998_Init
  public Cloudphys_IH1998_Forcing
  public Cloudphys_IH1998_Sounding
  
contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  subroutine IH1998_coefficient( temp )
    !
    !計算に利用する係数 \mu_{ijk}, \varepsilon_{ijk} を計算.
    !Zeleznik (1991), J. Phys. Chem. Data の Table 6 に基づく.  
    !係数は温度の関数. 
    !
    
    ! 暗黙の型宣言禁止
    implicit none
    
    ! 変数宣言
    real(8), intent(in) :: temp    !! 温度 (temperature)
    
    !H2SO4 水溶液のギブス自由エネルギーを求めるための係数.
    real(8), parameter :: mu111_a = -0.235245033870d+2
    real(8), parameter :: mu111_b =  0.406889449841d-1
    real(8), parameter :: mu111_c = -0.151369362907d-4
    real(8), parameter :: mu111_d =  0.296144445015d+4
    real(8), parameter :: mu111_e =  0.492476973663d+0
    real(8), parameter :: mu121_a =  0.111458541077d+4
    real(8), parameter :: mu121_b = -0.118330789360d+1
    real(8), parameter :: mu121_c = -0.209946114412d-2
    real(8), parameter :: mu121_d = -0.246749842271d+6
    real(8), parameter :: mu121_e =  0.341234558134d+2
    real(8), parameter :: mu221_a = -0.801488100747d+2
    real(8), parameter :: mu221_b = -0.116246143257d-1
    real(8), parameter :: mu221_c =  0.606767928954d-5
    real(8), parameter :: mu221_d =  0.309272150882d+4
    real(8), parameter :: mu221_e =  0.127601667471d+2
    real(8), parameter :: mu122_a =  0.888711613784d+3
    real(8), parameter :: mu122_b = -0.250531359687d+1
    real(8), parameter :: mu122_c =  0.605638824061d-3
    real(8), parameter :: mu122_d = -0.196985296431d+6
    real(8), parameter :: mu122_e =  0.745500643380d+2
    real(8), parameter :: eps111_a =  0.288731663295d+4
    real(8), parameter :: eps111_b = -0.332602457749d+1
    real(8), parameter :: eps111_c = -0.282047283300d-2
    real(8), parameter :: eps111_d = -0.528216112353d+6
    real(8), parameter :: eps111_e =  0.686997433564d+0
    real(8), parameter :: eps211_a =  0.383025318809d+2
    real(8), parameter :: eps211_b = -0.295997878789d-1
    real(8), parameter :: eps211_c =  0.120999746782d-4
    real(8), parameter :: eps211_d = -0.324697498999d+4
    real(8), parameter :: eps211_e = -0.383566039532d+1
    real(8), parameter :: eps121_a = -0.370944593249d+3
    real(8), parameter :: eps121_b = -0.690310834523d+0
    real(8), parameter :: eps121_c =  0.563455068422d-3
    real(8), parameter :: eps121_d = -0.382252997064d+4
    real(8), parameter :: eps121_e =  0.942682037574d+2
    real(8), parameter :: eps221_a =  0.232476399402d+4
    real(8), parameter :: eps221_b = -0.141626921317d+0
    real(8), parameter :: eps221_c = -0.626760562881d-2
    real(8), parameter :: eps221_d = -0.450590687961d+6
    real(8), parameter :: eps221_e = -0.612339472744d+2
    real(8), parameter :: eps122_a = -0.163385547832d+4
    real(8), parameter :: eps122_b = -0.335344369968d+1
    real(8), parameter :: eps122_c =  0.710978119903d-2
    real(8), parameter :: eps122_d =  0.198200003569d+6
    real(8), parameter :: eps122_e =  0.246693619189d+3 
    real(8), parameter :: eps212_a =  0.127375159848d+4
    real(8), parameter :: eps212_b =  0.103333898148d+1
    real(8), parameter :: eps212_c =  0.341400487633d-2
    real(8), parameter :: eps212_d =  0.195290667051d+6
    real(8), parameter :: eps212_e = -0.431737442782d+3
  
    !係数 \mu_ijk の計算
    mu111 =   mu111_a               &
      &     + mu111_b * temp        &
      &     + mu111_c * temp* temp  &
      &     + mu111_d / temp        &
      &     + mu111_e * dlog(temp)
    mu121 =   mu121_a               &
      &     + mu121_b * temp        &
      &     + mu121_c * temp* temp  &
      &     + mu121_d / temp        &
      &     + mu121_e * dlog(temp)
    mu221 =   mu221_a               &
      &     + mu221_b * temp        &
      &     + mu221_c * temp* temp  &
      &     + mu221_d / temp        &
      &     + mu221_e * dlog(temp)
    mu122 =   mu122_a               &
      &     + mu122_b * temp        &
      &     + mu122_c * temp* temp  &
      &     + mu122_d / temp        &
      &     + mu122_e * dlog(temp)
    mu211 = mu121
    mu212 = mu122
    
    !係数 \varepsilon_ijk の計算
    eps111 =   eps111_a              &
      &      + eps111_b * temp       &
      &      + eps111_c * temp* temp &
      &      + eps111_d / temp       &
      &      + eps111_e * dlog(temp)
    eps121 =   eps121_a              &
      &      + eps121_b * temp       &
      &      + eps121_c * temp* temp &
      &      + eps121_d / temp       &
      &      + eps121_e * dlog(temp)
    eps211 =   eps211_a              &
      &      + eps211_b * temp       &
      &      + eps211_c * temp* temp &
      &      + eps211_d / temp       &
      &      + eps211_e * dlog(temp)
    eps221 =   eps221_a              &
      &      + eps221_b * temp       &
      &      + eps221_c * temp* temp &
      &      + eps221_d / temp       &
      &      + eps221_e * dlog(temp)
    eps122 =   eps122_a              &
      &      + eps122_b * temp       &
      &      + eps122_c * temp* temp &
      &      + eps122_d / temp       &
      &      + eps122_e * dlog(temp)
    eps212 =   eps212_a              &
      &      + eps212_b * temp       &
      &      + eps212_c * temp* temp &
      &      + eps212_d / temp       &
      &      + eps212_e * dlog(temp)
    
  end subroutine IH1998_coefficient
  
  
  subroutine IH1998_SatPressRef( TempIN, SatPress1, SatPress2 )
    !
    !飽和蒸気圧(単成分)の計算.     
    !
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数宣言
    real(8), intent(in)    :: TempIN  !! 温度
                                      !! temperature
    real(8), intent(out)   :: SatPress1
    real(8), intent(out)   :: SatPress2
    real(8), parameter     :: alpha = 24.021415d0
    real(8), parameter     :: beta  = 4616.9134d0
    real(8), parameter     :: gamm  = 3.1934553d-4
    real(8), parameter     :: delta = 2.7550431d-11
    real(8), parameter     :: eps   = 1.0131374d-2
    real(8), parameter     :: ita   = 1.3158813d-2
    real(8)                :: xi
    real(8)                :: zeta
    real(8)                :: Temp
    real(8)                :: TempRefDivTemp
    real(8)                :: Log_SatPressRef1
    real(8)                :: Log_SatPressRef2    
    
    !温度は TempCr 以上は取れない.
    Temp = min(TempIN, TempCr)

    !!
    !! H2SO4 の飽和蒸気圧 
    !!

    ! 定数の設定
    xi   = ( ( Temp + 0.01d0 ) ** 2.0d0 ) - 293700.0d0
    zeta = (TempCr - Temp) ** 1.25d0  
    TempRefDivTemp = TempRef / Temp 
    
    ! H2SO4 (純物質) の飽和蒸気圧 (単位 : bar )
    ! なぜか log 2底なので....
    !    Log_SatPressRef1 = - ( 10156.0d0 / Temp ) + 16.259d0
    Log_SatPressRef1 = &         
      & +  16.259d0 &
      & -  10156.0d0 / Temp  &
      & +  7.42d0 * ( 1.0d0 + log( TempRefDivTemp ) / log(2.0d0) - TempRefDivTemp )
    
    ! H2SO4 (溶液) の飽和蒸気圧 (単位 : Pa)
    SatPress1 = exp( Log_SatPressRef1 ) * 1.0d5
    
    !!
    !! H20 の飽和蒸気圧 
    !!
    
    ! H2O (純物質) の飽和蒸気圧 (単位 : Pa)
    Log_SatPressRef2 = &
      &   alpha &
      & - beta / Temp &
      & + gamm * xi * ( dexp( delta * ( xi ** 2.0d0 ) ) - 1.0d0 ) / (Temp + 0.01d0) &
      & - eps * dexp( -1.0d0 * ita * zeta )
    
    ! H2O (溶液) の飽和蒸気圧
    SatPress2 = exp( Log_SatPressRef2 )
    
  end subroutine IH1998_SatPressRef
  

  subroutine IH1998_SatPress( Temp, con, SatPress1, SatPress2 )
    !
    ! H2SO4 水溶液の飽和蒸気圧
    !

    !暗黙の型宣言禁止
    implicit none

    !変数宣言
    real(8), intent(in)    :: con   !! H2SO4 水溶液の濃度
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: Temp  !! 温度
                                    !! temperature
    real(8), intent(out)   :: SatPress1
    real(8), intent(out)   :: SatPress2
    real(8)                :: SatPressRef1
    real(8)                :: SatPressRef2    
    real(8)                :: DelChemPotRDivRT1
    real(8)                :: DelChemPotRDivRT2 

    !純物質の飽和蒸気圧
    call IH1998_SatPressRef( Temp, SatPressRef1, SatPressRef2 )
    
    !混合の化学ポテンシャル変化を計算
    call IH1998_DelChemPotRDivRT( &
      & Temp, con, DelChemPotRDivRT1, DelChemPotRDivRT2 )

    ! H2SO4 (溶液) の飽和蒸気圧
    SatPress1 = SatPressRef1 * exp( DelChemPotRDivRT1 )
    
    ! H2O (溶液) の飽和蒸気圧
    SatPress2 = SatPressRef2 * exp( DelChemPotRDivRT2 )
    
  end subroutine IH1998_SatPress
  

  subroutine IH1998_GibbsRDivRT( temp, con, GibbsRDivRT ) 
    !
    ! - G^{(r)} / RT の計算.
    ! G^{(r)} = G - G^{\circ}(T,p,x)
    ! (G^{(r)} : relative molar Gibbs energy of aquenous sulfuric asid)
    !
    
    !暗黙の型宣言禁止
    implicit none

    !引数の宣言
    real(8), intent(in)    :: con   !! H2SO4 水溶液の濃度
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: temp  !! 温度
                                    !! temperature
    real(8), intent(out)   :: GibbsRDivRT

    
    !変数の宣言
    real(8)                :: x1, x2

    !モル比がマシンイプシロンより小さい場合の処理.
    !\ln(0.0) は発散するので. 
    x1 = max( epsilon(con), con )
    x1 = min( 1.0d0 - epsilon(x1), x1 )

    ! H2SO4 と H2O のモル比
    x2 = 1.0d0 - x1

    !係数を決める
    call IH1998_coefficient(temp)   

    !-\frac{G^{(r)}}{RT}
    !  = \sum^2_{i=1} \Phi_{i} \sum^2_{j=1}\sum^2_{k=1}
    !     (\mu_{jki} + \varepsilon_{jki} \ln x_{j}) x_j x_k
    !
    ! mu112 = mu222 = eps112 = eps122 = 0.0
    ! Ph1_1 = 1, Phi2 = x1 * x2
    !
    GibbsRDivRT =                                                 &
      & -1.0d0 * (                                                &
      &    (                                                      & !i=1
      &       + ( mu111 + eps111 * dlog(x1) ) * x1 * x1           &   !j=1,k=1
      &       + ( mu121 + eps121 * dlog(x1) ) * x1 * x2           &   !j=1,k=2
      &       + ( mu211 + eps211 * dlog(x2) ) * x2 * x1           &   !j=2,k=1
      &       + ( mu221 + eps221 * dlog(x2) ) * x2 * x2           &   !j=2,k=2
      &    )                                                      &
      &  +                                                        &
      &    (                                                      & !i=2
      &       + ( mu122 + eps122 * dlog(x1) ) * x1 * x1 * x2 * x2 &   !j=1,k=2
      &       + ( mu212 + eps212 * dlog(x2) ) * x1 * x1 * x2 * x2 &   !j=2,k=1
      &    )                                                      &
      & )
    
  end subroutine IH1998_GibbsRDivRT
  

  subroutine IH1998_DelChemPotRDivRT( &
    & temp, con, DelChemPotRDivRT1, DelChemPotRDivRT2 ) 
    !
    ! \frac{\Del\mu_{i}}{RT} = \frac{\mu_i(T,x_i) - \mu^{\circ}(T)}{RT}
    ! = (G^{(r)} - x_1 \DP{G^{(r)}}{x_1} - x_2 \DP{G^{(r)}}{x_2} + \DP{G^{(r)}}{x_i})/RT
    !    
    
    !暗黙の型宣言禁止
    implicit none

    !引数の宣言
    real(8), intent(in)    :: con   !! H2SO4 水溶液の濃度
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: temp  !! 温度
                                    !! temperature
    real(8), intent(out)   :: DelChemPotRDivRT1
    real(8), intent(out)   :: DelChemPotRDivRT2

    !変数の宣言 
    real(8)                :: GibbsRDivRT
    real(8)                :: DGibbsRDivRT_Dx1
    real(8)                :: DGibbsRDivRT_Dx2
    real(8)                :: x1, x2
    
    !モル比がマシンイプシロンより小さい場合の処理.
    !\ln(0.0) は発散するので. 
    x1 = max( epsilon(con), con )
    x1 = min( 1.0d0 - epsilon(x1), x1 )
    
    ! H2SO4 と H2O のモル比
    x2 = 1.0d0 - x1

    !係数を決める
    call IH1998_coefficient(temp)

    ! 関数の呼び出し
    call IH1998_GibbsRDivRT(temp, con, GibbsRDivRT)
   
    ! \DP{G^{(r)}/RT}{x_1} 
    DGibbsRDivRT_Dx1 =                                               &
         & -1.0d0 * (                                                &
         &  (                                                        & !i=1
         &    + ( mu111 + eps111 * dlog(x1) ) * x1 * 2.0d0           &   !j=1,k=1 
         &              + eps111 * x1                                &
         &    + ( mu121 + eps121 * dlog(x1) ) * x2                   &   !j=1,k=2
         &              + eps121 * x2                                &
         &    + ( mu211 + eps211 * dlog(x2) ) * x2                   &   !j=2,k=1
         &  )                                                        &
         &  +                                                        &
         &  (                                                        & !i=2
         &    + ( mu122 + eps122 * dlog(x1) ) * x1 * x2 * x2 * 2.0d0 &   !j=1,k=2
         &              + eps122 * x1 * x2 * x2                      &
         &    + ( mu212 + eps212 * dlog(x2) ) * x1 * x2 * x2 * 2.0d0 &   !j=2,k=1
         &  )                                                        &
         & )

    ! \DP{G^{(r)}/RT}{x_2} 
    DGibbsRDivRT_Dx2 =                                               &
         & -1.0d0 * (                                                &
         &  (                                                        & !i=1
         &    + ( mu121 + eps121 * dlog(x1) ) * x1                   &   !j=1,k=2
         &    + ( mu211 + eps211 * dlog(x2) ) * x1                   &   !j=2,k=1
         &              + eps211 * x1                                &
         &    + ( mu221 + eps221 * dlog(x2) ) * x2 * 2.0d0           &   !j=2,k=2
         &              + eps221 * x2                                &
         &  )                                                        &
         &  +                                                        &
         &  (                                                        & !i=2
         &   + ( mu122 + eps122 * dlog(x1) ) * x1 * x1 * x2 * 2.0d0  &   !j=1,k=2
         &   + ( mu212 + eps212 * dlog(x2) ) * x1 * x1 * x2 * 2.0d0  &   !j=2,k=1
         &             + eps212 * x1 * x1 * x2                       &
         &  )                                                        &
         & )
    
    DelChemPotRDivRT1 =                   &
         &  min(   0.0d0,                   &
         &       + GibbsRDivRT            &
         &       - x1 * DGibbsRDivRT_Dx1  &
         &       - x2 * DGibbsRDivRT_Dx2  &
         &       + DGibbsRDivRT_Dx1 )

    DelChemPotRDivRT2 =   &
         &  min(   0.0d0, &
         &       + GibbsRDivRT            &
         &       - x1 * DGibbsRDivRT_Dx1  &
         &       - x2 * DGibbsRDivRT_Dx2  &
         &       + DGibbsRDivRT_Dx2 )
    
  end subroutine IH1998_DelChemPotRDivRT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

  subroutine IH1998_DDelChemPotRDivRTDx( &
       & temp, con, DDelChemPotRDivRT1_Dx1, DDelChemPotRDivRT2_Dx1 )
    !
    ! \frac{\Del\mu_{i}}{RT} = \frac{\mu_i(T,x_i) - \mu^{\circ}(T)}{RT}
    ! = (G^{(r)} - x_1 \DP{G^{(r)}}{x_1} - x_2 \DP{G^{(r)}}{x_2} + \DP{G^{(r)}}{x_i})/RT
    ! x1 微分
    !
    
    !暗黙の型宣言禁止
    implicit none

    !引数の宣言
    real(8), intent(in)    :: con   !! H2SO4 水溶液の濃度
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: temp  !! 温度
                                    !! temperature
    real(8), intent(out)   :: DDelChemPotRDivRT1_Dx1
    real(8), intent(out)   :: DDelChemPotRDivRT2_Dx1

    !変数の宣言
    real(8)                :: DGibbsRDivRT_Dx1Dx1
    real(8)                :: DGibbsRDivRT_Dx1Dx2
!    real(8)                :: DGibbsRDivRT_Dx2Dx2
    real(8)                :: x1, x2
    
    !モル比がマシンイプシロンより小さい場合の処理.
    !\ln(0.0) は発散するので. 
    x1 = max( epsilon(con), con )
    x1 = min( 1.0d0 - epsilon(x1), x1 )
    
    ! H2SO4 と H2O のモル比
    x2 = 1.0d0 - x1

    !係数を決める
    call IH1998_coefficient(temp)

    ! DP{\DP{G^{(r)}/RT}{x_1}}{x_1}
    DGibbsRDivRT_Dx1Dx1 =                                                         &
         &    - ( mu111 + eps111 * dlog(x1) ) * 2.0d0                             &   !j=1,k=1 
         &    -           eps111              * 3.0d0                             &
         &    + ( mu121 + eps121 * dlog(x1) )                                     &   !j=1,k=2
         &    -           eps121              * ( x2 / x1  - 1.0d0 )              &
         &    + ( mu211 + eps211 * dlog(x2) )                                     &   !j=2,k=1
         &    +           eps211                                                  &
         &    - ( mu122 + eps122 * dlog(x1) ) * x2 * (1.0d0 - x1 * 3.0d0) * 2.0d0 &   !j=1,k=2
         &    -           eps122              * x2 * (3.0d0 - 5.0d0 * x1)         &
         &    - ( mu212 + eps212 * dlog(x2) ) * x2 * (1.0d0 - x1 * 3.0d0) * 2.0d0 &   !j=2,k=1
         &    +           eps212              * x1 * x2 * 2.0d0                      !j=2,k=1

    ! \DP{\DP{G^{(r)}/RT}{x_1}}{x_2}
    DGibbsRDivRT_Dx1Dx2 =                                                         &
         &    - ( mu121 + eps121 * dlog(x1) )                                     &   !j=1,k=2
         &    -           eps121                                                  &
         &    - ( mu211 + eps211 * dlog(x2) )                                     &   !j=2,k=1
         &    +           eps211              * (x1 / x2 - 1.0d0)                 &
         &    + ( mu221 + eps221 * dlog(x2) ) * 2.0d0                             &   !j=2,k=2
         &    +           eps221              * 3.0d0                             &
         &    - ( mu122 + eps122 * dlog(x1) ) * x1 * (2.0d0 - 3.0d0 * x1) * 2.0d0 &   !j=1,k=2
         &    -           eps122              * x1 * x2 * 2.0d0                   &   !j=1,k=2
         &    - ( mu212 + eps212 * dlog(x2) ) * x1 * (2.0d0 - 3.0d0 * x1) * 2.0d0 &   !j=2,k=1
         &    +           eps212              * x1 * ( 5.0d0 * x1 - 2.0d0 )

    DDelChemPotRDivRT1_Dx1 = + x2 * DGibbsRDivRT_Dx1Dx1  &
         &                   - x2 * DGibbsRDivRT_Dx1Dx2  
    
    DDelChemPotRDivRT2_Dx1 = - x1 * DGibbsRDivRT_Dx1Dx1  &
         &                   + x1 * DGibbsRDivRT_Dx1Dx2
    
  end subroutine IH1998_DDelChemPotRDivRTDx
  
  

  subroutine IH1998_newton( NumDens1, NumDens2, Temp, con, SatPress1, SatPress2, Flag )

    !暗黙の型宣言
    implicit none

    !変数宣言
    real(8), intent(inout) :: con   !! H2SO4 水溶液の濃度
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: Temp  !! 温度
                                    !! temperature
    real(8), intent(in)    :: NumDens1, NumDens2
    logical, intent(out)   :: Flag
    real(8)                :: Func, DFuncDx1
    real(8)                :: SatPress1, SatPress2
    real(8)                :: DSatPress1Dx1, DSatPress2Dx1
    real(8)                :: DDelChemPotRDivRT1_Dx1, DDelChemPotRDivRT2_Dx1
    real(8)                :: BoltzTemp
    real(8)                :: x1, x1A, x1P, x1L, x1R
    integer                :: CalNum
    integer, parameter     :: MaxNum = 10
    real(8)                :: x1save
    
    !初期化
    CalNum = 0
    x1save = con
    
    !定数の設定
    BoltzTemp = Boltz * Temp
    
    !変数名の置き換え
    x1 = con
    
    Loop: do
      !ループ回数
      CalNum = CalNum + 1
      
      !水溶液の飽和蒸気圧
      call IH1998_SatPress( Temp, x1, SatPress1, SatPress2 )
      
      !混合の化学ポテンシャルの x1 微分を計算
      call IH1998_DDelChemPotRDivRTDx( &
        & Temp, x1, DDelChemPotRDivRT1_Dx1, DDelChemPotRDivRT2_Dx1 )
      
      !H2SO4 (溶液) の飽和蒸気圧の x1 (H2SO4 モル比) 微分
      DSatPress1Dx1 = DDelChemPotRDivRT1_Dx1 * SatPress1
      
      !H2O (溶液) の飽和蒸気圧の x1 (H2SO4 モル比) 微分
      DSatPress2Dx1 = DDelChemPotRDivRT2_Dx1 * SatPress2       
      
      !関数
      Func = &
        &   ( 1.0d0 - x1 ) * (NumDens1 - SatPress1 / BoltzTemp ) &
        & - x1             * (NumDens2 - SatPress2 / BoltzTemp )
      
      !関数の x1 に関する 1 階微分
      DFuncDx1 = &
        & - ( NumDens1 - SatPress1 / BoltzTemp ) &
        & - ( NumDens2 - SatPress2 / BoltzTemp ) &
        & - ( 1.0d0 - x1 ) * DSatPress1Dx1 / BoltzTemp  &
        & +   x1           * DSatPress2Dx1 / BoltzTemp         
      
      !真に近い解. 値は 0 ~ 1 の範囲. 
      x1A =  x1 - Func / DFuncDx1
      
      !次の値が x1 の満たすべき範囲 0 <= x1 <= 1 を超えた場合の処置.
      if (x1A >= 1.0d0 ) then 
        !          write(*,*) "CALL Bisection (1)"
        x1L = x1
        x1R = 1.0d0
        call IH1998_bisection  &
          & (NumDens1, NumDens2, Temp, x1L, x1R, con, SatPress1, SatPress2, Flag)
        exit Loop
      end if
!       if (x1A <= 0.0d0 ) then 
!          write(*,*) "CALL Bisection (2)"
!          x1L = 0.0d0
!          x1R = x1
!          call IH1998_bisection(NumDens1, NumDens2, Temp, x1L, x1R, con, Flag)
!          exit Loop
!       end if

      !ループの終了条件
      if ( abs((x1A - x1) / x1) < 1.0d-12 .OR. CalNum > MaxNum ) then 
        
        if( NumDens1 - SatPress1 / BoltzTemp > 0 .AND. NumDens2 - SatPress2 / BoltzTemp > 0) then
          !
          ! 凝結の条件が OK の時の処理
          !
          con = x1A
          Flag = .true.
          
        elseif(      abs(NumDens1 - SatPress1 / BoltzTemp) / NumDens1 < 1.0d-7 &
          &     .OR. abs(NumDens2 - SatPress2 / BoltzTemp)/NumDens2 < 1.0d-7   &
          &     .OR. CalNum > MaxNum                                           &
          &     .OR. x1A > 1.0d0                                               &
          &     .OR. x1A < 0.0d0                                               ) then
          !
          ! 精度を確保するために, 二分法でもダメか念押しする.
          !
          !             write(*,*) "CALL Bisection (3)"
          x1L = max(0.0d0, x1save - 1.0d-1)
          x1R = min(1.0d0, x1save + 1.0d-1)
          call IH1998_bisection &
            & (NumDens1, NumDens2, Temp, x1L, x1R, con, SatPress1, SatPress2, Flag)
          
        else
          !
          ! 凝結の条件を満たせない場合の処理
          !
          con = -1.0d0
          Flag = .true.
        end if
        exit
      end if
      
      !ループを回すための処理
      x1P = x1
      x1  = x1A
      
    end do Loop
    
  end subroutine IH1998_Newton
  
  
  subroutine IH1998_bisection &
    & ( NumDens1, NumDens2, Temp, x1L, x1R, x1C, SatPress1, SatPress2, Flag )

    !暗黙の型宣言禁止
    implicit none

    !変数宣言
    real(8), intent(inout) :: x1L, x1R, x1C    !! H2SO4 のモル比 
    real(8), intent(in)    :: Temp  !! 温度
                                    !! temperature
    real(8), intent(in)    :: NumDens1, NumDens2

    real(8), intent(out)   :: SatPress1, SatPress2
    logical, intent(out)   :: Flag
    real(8)                :: FuncL, FuncC, FuncR
    real(8)                :: BoltzTemp
    integer                :: CalNum

    !初期化
    CalNum = 0
    Flag = .false.
    
    !定数の設定
    BoltzTemp = Boltz * Temp

    !以下の条件が同時に満たされないといけない.
    ! NumDens1 - SatPress1 / BoltzTemp > 0
    ! NumDens2 - SatPress2 / BoltzTemp > 0
    
    Loop: do
      !ループ回数
      CalNum = CalNum + 1
      
      !中点. 真に近い解
      x1C = ( x1R + x1L ) * 5.0d-1
!       write(*,*) CalNum, real(x1C), real(abs((x1R - x1L) / x1R))

      !関数
      call IH1998_SatPress( Temp, x1L, SatPress1, SatPress2 )
      FuncL = &
        &  ( 1.0d0 - x1L ) * max( (NumDens1 - SatPress1 / BoltzTemp ), 0.0d0) &
        & - x1L            * max( (NumDens2 - SatPress2 / BoltzTemp ), 0.0d0)
      
      call IH1998_SatPress( Temp, x1R, SatPress1, SatPress2 )
      FuncR = &
        &  ( 1.0d0 - x1R ) * max( (NumDens1 - SatPress1 / BoltzTemp ), 0.0d0 ) &
        & - x1R            * max( (NumDens2 - SatPress2 / BoltzTemp ), 0.0d0 )
      
      call IH1998_SatPress( Temp, x1C, SatPress1, SatPress2 )
      FuncC = &
        &  ( 1.0d0 - x1C ) * max( (NumDens1 - SatPress1 / BoltzTemp ), 0.0d0 ) &
        & - x1C            * max( (NumDens2 - SatPress2 / BoltzTemp ), 0.0d0 )
      
      !ループを回すための処理
      if ( FuncL * FuncC < 0.0d0 ) then 
        x1R = x1C          
      elseif ( FuncC * FuncR < 0.0d0 ) then
        x1L = x1C
      end if
      
      !ループの終了条件
      if ( abs((x1R - x1L) / x1R) < 1.0d-12 .OR. CalNum >= 60 .OR. FuncC == 0.0d0 ) then 
        
        if( NumDens1 - SatPress1 / BoltzTemp > 0 .AND. NumDens2 - SatPress2 / BoltzTemp > 0) then
          ! 凝結の条件を満たしている場合は x1C を戻り値とする. 
          flag = .true.
          
        else
          ! 凝結の条件を満たしていない場合は flag を .false. にする.
          x1C = -1.0d0
          flag = .false.
        end if
        exit
      end if
      
       
    end do Loop
    
  end subroutine IH1998_Bisection
  

  subroutine IH1998_EquivState &
    & ( Temp, n1, n2, n1_gas, n2_gas, n1_liq, n2_liq, con, sw, flag )
    !
    ! 平衡状態を得るためのプログラム.
    !

    !暗黙の型宣言禁止
    implicit none

    !変数宣言
    real(8), intent(in)  :: Temp
    real(8), intent(in)  :: n1, n2
    integer, intent(in)  :: sw
    real(8), intent(out) :: n1_gas, n2_gas, n1_liq, n2_liq, con
    logical, intent(out) :: flag    
    
    integer              :: j 
    integer, parameter   :: NumTest = 6
    integer              :: idx(1), cnt(1)
    real(8)              :: FuncSign(NumTest)
    real(8)              :: FuncTest1, FuncTest2
    real(8)              :: del 
    real(8), parameter   :: TempCr= 647.26d0
    real(8)              :: x1N   
    real(8)              :: SatPressRef1, SatPressRef2
    real(8)              :: SatPress1, SatPress2
    real(8)              :: x1R, x1L

    !初期化
    flag = .true.
    del = 1.0d0 / real(NumTest - 1)
    
    !データのあたりを付ける.      
    do j = 1, NumTest
      con = real(j-1) * del
      call IH1998_SatPress( Temp, con, SatPress1, SatPress2)
      
      FuncTest2 = &
        &   ( 1.0d0 - con ) * max( (n1 - SatPress1 / ( boltz * Temp ) ), 0.0d0 ) &
        & -   con           * max( (n2 - SatPress2 / ( boltz * Temp ) ), 0.0d0 )
      
      if ( j==1 ) FuncTest1 = FuncTest2
      
      FuncSign(j)= FuncTest1 * FuncTest2
      FuncTest1  = FuncTest2        
    end do
    
    ! 0 <= x1 <= 1 の範囲に解があるのなら,
    ! dx 間隔で隣同士の掛け算 (ex. Func(0) * Func(0.2)) で負になるのは一箇所.
    ! 負になる箇所の配列添え字は 1 以下にならない. 
    !
    cnt = count( FuncSign < 0.0d0)
    idx = minloc(FuncSign, FuncSign < 0.0d0)
    
    ! もしも 0 <= x1 <= 1 の範囲に解があるという条件を満たさないならば,
    ! 水溶液は生成されない. 純物質の飽和条件のみ考慮すれば良い. 
    if ( (cnt(1) /= 1 .OR. idx(1) <= 1) .OR. n1 == 0.0d0 .OR. n2 == 0.0d0 ) then
      !
      !水溶液が凝結しない場合.
      !
      
      call IH1998_SatPressRef( Temp, SatPressRef1, SatPressRef2)
      
      con = -1.0d0
      
      if ( n1 > (SatPressRef1 / ( boltz * Temp ))) then 
        n1_gas =      SatPressRef1 / ( boltz * Temp )
        n1_liq = n1 - SatPressRef1 / ( boltz * Temp ) 
        n2_gas = n2
        n2_liq = 0.0d0
      elseif ( n2 > (SatPressRef2 / ( boltz * Temp ))) then 
        n1_gas = n1
        n1_liq = 0.0d0
        n2_gas =      SatPressRef2 / ( boltz * Temp )
        n2_liq = n2 - SatPressRef2 / ( boltz * Temp )
      else
        n1_gas = n1
        n1_liq = 0.0d0
        n2_gas = n2
        n2_liq = 0.0d0
      end if
      
    else     
      
      !
      !水溶液が凝結する場合
      !
      
      !二分法の初期値
      x1L = real(idx(1)-2) * del
      x1R = real(idx(1)-1) * del
      
      if ( sw == 1 ) then 
        !二分法で計算
        call IH1998_bisection(n1, n2, Temp, x1L, x1R, con, SatPress1, SatPress2, flag)
      else
        ! ニュートン法の初期値
        con = (real(idx(1)-2) + real(idx(1)-1))* del * 5.0d-1
        call IH1998_newton(n1, n2, Temp, con, SatPress1, SatPress2, flag)
      end if
      
      ! 凝結量の計算
      if ( n1 > (SatPress1 / ( boltz * Temp )) .AND. n2 > (SatPress2 / ( boltz * Temp ))) then 
        n1_gas =       SatPress1 / ( boltz * Temp )
        n1_liq =  n1 - SatPress1 / ( boltz * Temp )
        n2_gas =       SatPress2 / ( boltz * Temp )
        n2_liq =  n2 - SatPress2 / ( boltz * Temp )                   
      else
        n1_gas = n1
        n1_liq = 0.0d0
        n2_gas = n2
        n2_liq = 0.0d0
      end if
    end if
    
  end subroutine IH1998_EquivState
  

  subroutine IH1998_Sediment(Height, VelZ)
    !    
    ! 落下項の計算. 高度のみの関数. 
    !
    
    !暗黙の型宣言禁止
    implicit none

    !変数宣言
    real(8), intent(in)  :: Height
    real(8), intent(out) :: VelZ
    real(8), parameter   :: Grav = 8.7d0
    real(8), parameter   :: Dens = 1.8d3
    real(8), parameter   :: R_2  = 1.15d-6  !モード2 の半径
    real(8), parameter   :: R_3  = 3.65d-6
    real(8), parameter   :: Vis  = 1.5d-5
    real(8)              :: Molfr_2, Molfr_3
    real(8)              :: VelZ_2,  VelZ_3
    
    !mode3 の分布を tanh で近似
    Molfr_3 = 0.455d0 * tanh( ( 57.5d3 - Height ) * 1.0d-3 ) + 0.455d0
    
    !mode2 の分布
    Molfr_2 = 1.d0 - Molfr_3
    
    !mode 2, 3 の速度
    VelZ_3  = - 2.0d0 * Grav * Dens * R_3 * R_3 / Vis / 9.0d0
    VelZ_2  = - 2.0d0 * Grav * Dens * R_2 * R_2 / Vis / 9.0d0
    
    !平均速度
    VelZ = Molfr_2 * VelZ_2 + Molfr_3 * VelZ_3
    
  end subroutine IH1998_Sediment
  
  
  subroutine IH1998_H2SO4Prdt( Height, DelZ, H2SO4_Prdt, H2O_Loss )
    !
    ! H2SO4 の生成項 & H2O の消滅項
    !
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数宣言
    real(8), intent(in) :: Height
    real(8), intent(in) :: DelZ
    real(8), intent(out):: H2SO4_Prdt 
    real(8), intent(out):: H2O_Loss
    
    real(8), parameter  :: PI = 3.1415926535897932d0
    real(8), parameter  :: phi = 30.0d0 * PI / 180.0d0
    real(8), parameter  :: sigma = 1.5d3
    real(8), parameter  :: Height0 = 62.0d3
    real(8)             :: DH
    
    DH = Height - Height0
    
    !単位を揃えるために DelZ で割り算しておく必要があるはず.
    !
    H2SO4_Prdt = & 
      & 3.5d0  * 1.0d16  / sigma / sqrt( 2.0d0 * PI ) &
      & * exp( -1.0d0 * DH * DH  / ( 2.0d0 * sigma * sigma ) ) &
      & * cos( Phi ) / DelZ
    
    H2O_Loss = -1.0d0 * H2SO4_Prdt
    
  end subroutine IH1998_H2SO4Prdt
  
  
  subroutine IH1998_H2SO4Loss( Temp, Press, n1, n2, H2SO4_Loss, H2O_Prdt )
    !
    ! H2SO4 の消滅項 & H2O の生成項
    !
    
    !暗黙の型宣言禁止
    implicit none
    
    !変数宣言
    real(8), intent(in) :: Temp
    real(8), intent(in) :: Press
    real(8), intent(in) :: n1, n2      !数密度
    real(8), intent(out):: H2SO4_Loss
    real(8), intent(out):: H2O_Prdt
    
    real(8), parameter  :: k13 = 5.5d-29
    real(8), parameter  :: f_co = 25.0d-6  !緯度 > 35
!     real(8), parameter  :: f_co = 30.0d-6  !35 < 緯度 < 60
!     real(8), parameter  :: f_co = 35.0d-6  !緯度 > 60
    real(8)              :: mean
    real(8)              :: kp
    
    !作業変数
    mean = press / (Boltz * Temp)
    kp = 10.0d0 ** ( 7.584d0 - ( 5060.0d0 / Temp ) )
    
    !H2SO4 の消滅項 (圧力は bar 単位に換算)
    H2SO4_Loss = - f_co * k13 * kp * mean * mean  * n1 / n2 / ( press * 1.0d-5 )
    
    !H2O の生成項
    H2O_Prdt = -1.0d0 * H2SO4_Loss
    
  end subroutine IH1998_H2SO4Loss
  
!!!-----------------------------------------------------------------------

  subroutine Cloudphys_IH1998_init

    ! モジュール呼び出し
    use gtool_historyauto, only: HistoryAutoAddVariable

    !暗黙の型宣言禁止
    implicit none

    !ヒストリの準備
    call HistoryAutoAddVariable(                     &
      & varname='H2O', dims=(/'x','y','z','t'/),  &
      & longname='Number Density of H2O (all)',      &
      & units='m-3',                                 &
      & xtype='float')

    call HistoryAutoAddVariable(                     &
      & varname='H2SO4', dims=(/'x','y','z','t'/),   &
      & longname='Number Density of H2SO4 (all)',    &
      & units='m-3',                                 &
      & xtype='float')

    call HistoryAutoAddVariable(                     &
      & varname='H2Ogas', dims=(/'x','y','z','t'/),  &
      & longname='Number Density of H2O (gas)',      &
      & units='m-3',                                 &
      & xtype='float')

    call HistoryAutoAddVariable(                     &
      & varname='H2SO4gas', dims=(/'x','y','z','t'/),&
      & longname='Number Density of H2SO4 (gas)',    &
      & units='m-3',                                 &
      & xtype='float')

    call HistoryAutoAddVariable(                     &
      & varname='H2Oliq', dims=(/'x','y','z','t'/),  &
      & longname='Number Density of H2O (liq)',      &
      & units='m-3',                                 &
      & xtype='float')

    call HistoryAutoAddVariable(                     &
      & varname='H2SO4liq', dims=(/'x','y','z','t'/),&
      & longname='Number Density of H2SO4 (liq)',    &
      & units='m-3',                                 &
      & xtype='float')
    
  end subroutine Cloudphys_IH1998_init


  
  subroutine Cloudphys_IH1998_Sounding( &
    & xyzf_NDens1, xyzf_NDens2          & !(out)
    )
    !
    !金星雲微物理 (Imamura and Hashimoto, 1998) 用の初期値生成
    !初期値ファイルから初期値を生成する.
    !

    ! モジュール呼び出し
    !
    use dc_types,          only : STRING, DP
    use setmargin,         only : SetMargin_xyz
    use gtool_historyauto, only : HistoryAutoAddVariable
    use gridset,           only : imin, imax, jmin, jmax, kmin, kmax, &
      &                           nx, ny, nz
    use axesset,           only : r_Z, z_Z
  
    !暗黙の型宣言禁止
    implicit none

    !引数
    real(DP), intent(out) :: xyzf_NDens1(imin:imax, jmin:jmax, kmin:kmax, 1:3)
    real(DP), intent(out) :: xyzf_NDens2(imin:imax, jmin:jmax, kmin:kmax, 1:3)
    
    !変数宣言
    real(DP) :: r_alt(0:nz)
    real(DP) :: r_H2SO4(0:nz)
    real(DP) :: r_H2SO4gas(0:nz)
    real(DP) :: r_H2O(0:nz)
    real(DP) :: r_H2Ogas(0:nz)
    integer  :: i, j, k

    !初期化
    xyzf_NDens1 = 0.0d0
    xyzf_NDens2 = 0.0d0
    
    !初期値ファイル読み込み. 
    open (17, file='H2SO4-H2O_initial.dat', status='old')
!    read (17, '()')       ! ヘッダ行の読み飛ばし

    ! データ読み込み
    do k = 1, nz
      read (17, *) r_alt(k), r_H2SO4(k), r_H2SO4gas(k), r_H2O(k), r_H2Ogas(k)
!      write(*,*) "**", r_alt(k), r_H2SO4(k), r_H2O(k)
    end do

    ! 添字 0 での値を決める
    r_alt(0)      = r_alt(1) - 125.0d0
    r_H2SO4(0)    = r_H2SO4(1)
    r_H2SO4gas(0) = r_H2SO4gas(1)
    r_H2O(0)      = r_H2O(1)
    r_H2Ogas(0)   = r_H2Ogas(1)

    !ガス成分だけ利用する. 位置合わせのために平均する (z -> r). 
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          xyzf_NDens1(i,j,k,1) = ( r_H2SO4gas(k) + r_H2SO4gas(k-1) ) * 0.5d0
          xyzf_NDens2(i,j,k,1) = ( r_H2Ogas(k)   + r_H2Ogas(k-1)   ) * 0.5d0
        end do
      end do
    end do

    !現在は格子間隔を決め打ち
    if ( r_alt(2) - r_alt(1) /= 125.0d0 .OR. r_alt(1) /= r_Z(1) .OR. r_alt(2) /= r_Z(2) ) then
      write(*,*) "XXXXXXXXXX  NG, dz = ", r_alt(2) - r_alt(1)
      write(*,*)  r_alt(1), r_alt(2), r_alt(3), r_alt(4), r_alt(5)

      write(*,*) "READ:        ", r_alt(1), r_alt(2)
      write(*,*) "READ:        ", r_alt(1) - 62.5d0, r_alt(1) + 62.5d0
      write(*,*) "MODEL (r_Z): ", r_z(1), r_z(2)
      write(*,*) "MODEL (z_Z): ", z_z(1), z_z(2)      
    end if
    
    !境界条件
    call SetMargin_xyz( xyzf_NDens1(:,:,:,1) )
    call SetMargin_xyz( xyzf_NDens2(:,:,:,1) )
    
    !全量と気体の値を同じにする. 液相はゼロとする. 
    xyzf_NDens1(:,:,:,2) = xyzf_NDens1(:,:,:,1)
    xyzf_NDens2(:,:,:,2) = xyzf_NDens2(:,:,:,1)

    !液相はゼロとする. 
    xyzf_NDens1(:,:,:,3) = 0.0d0
    xyzf_NDens2(:,:,:,3) = 0.0d0
    
  end subroutine Cloudphys_IH1998_Sounding
  

  subroutine Cloudphys_IH1998_forcing(         &
    &   xyz_PtempAl,  xyz_PtempBZ,             & !(in)
    &   xyz_ExnerNl,  xyz_ExnerBZ,             & !(in)
    &   pyz_VelXNl,   xqz_VelYNl,  xyr_VelZNl, & !(in)
    &   xyz_Khbl,                              & !(in)
    &   xyzf_NDens1B, xyzf_NDens2B,            & !(in)
    &   xyzf_NDens1N, xyzf_NDens2N,            & !(in)
    &   xyzf_NDens1A, xyzf_NDens2A             & !(out)
    & )
    !
    ! 数密度の時間発展方程式を解く
    !

    ! モジュール呼び出し
    !
    use dc_types,     only : STRING, DP
    use timeset,      only : DelTimeLong
    use setmargin,    only : SetMargin_xyz
    use constants,    only : PressBasis,        &!温位の基準圧力
      &                      CpDry,             &!乾燥成分の比熱
      &                      GasRDry             !乾燥成分の気体定数
    use axesset,      only :  z_Z,              &!Z 座標
      &                       dx, dy, dz         !格子間隔
    use gridset,      only : FlagCalc3D,        &
      &                      imin, imax, jmin, jmax, kmin, kmax, &
      &                      nx, ny, nz
    use differentiate_center2,  &
      &           only: xyz_dx_pyz, xyz_dy_xqz, xyz_dz_xyr, &
      &                 pyz_dx_xyz, xqz_dy_xyz, xyr_dz_xyz, &
      &                 pqz_dx_xqz, pqz_dy_pyz, pyz_dy_pqz, &
      &                 pyr_dx_xyr, pyz_dz_pyr, pyr_dz_pyz, &
      &                 pyr_dz_pyz, pyr_dx_xyr, xyr_dx_pyr, &
      &                 xqr_dz_xqz, xqr_dy_xyr, xyr_dy_xqr, &
      &                 pqz_dy_pyz, pqz_dx_xqz, xqz_dx_pqz, &
      &                 xqr_dy_xyr, xqr_dz_xqz, xqz_dz_xqr      
    use average,  only: xyz_pyz, xyr_pyr, xqz_pqz, &
      &                 pyz_xyz, pyr_xyr, pqz_xqz, &
      &                 xyz_xqz, pyz_pqz, xyr_xqr, &
      &                 xqz_xyz, pqz_pyz, xqr_xyr, &
      &                 xyz_xyr, pyz_pyr, xqz_xqr, &
      &                 xyr_xyz, pyr_pyz, xqr_xqz, &
      &                 pqz_xyz, pyr_xyz, xqr_xyz, &
      &                 xyz_pqz, xyz_pyr, xyz_xqr
    
    !暗黙の型宣言禁止
    implicit none

    !変数の定義
    real(DP), intent(in)  :: xyz_PtempAl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyz_PtempBZ(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyz_ExnerNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyz_ExnerBZ(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xqz_VelyNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: pyz_VelXNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyr_VelZNl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyz_KhBl(imin:imax,jmin:jmax,kmin:kmax)
    real(DP), intent(in)  :: xyzf_NDens1B(imin:imax,jmin:jmax,kmin:kmax,1:3)
    real(DP), intent(in)  :: xyzf_NDens2B(imin:imax,jmin:jmax,kmin:kmax,1:3)
    real(DP), intent(in)  :: xyzf_NDens1N(imin:imax,jmin:jmax,kmin:kmax,1:3)
    real(DP), intent(in)  :: xyzf_NDens2N(imin:imax,jmin:jmax,kmin:kmax,1:3)
    real(DP), intent(out) :: xyzf_NDens1A(imin:imax,jmin:jmax,kmin:kmax,1:3)
    real(DP), intent(out) :: xyzf_NDens2A(imin:imax,jmin:jmax,kmin:kmax,1:3)

    real(DP) :: xyz_Dn1Dt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn1Dt_prdt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn1Dt_loss(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn1Dt_fall(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_Dn1DtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn1DtDif(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn1DtTrb(imin:imax,jmin:jmax,kmin:kmax)
    
    real(DP) :: xyz_Dn2Dt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn2Dt_prdt(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn2Dt_loss(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn2Dt_fall(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_Dn2DtAdv(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn2DtDif(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Dn2DtTrb(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_Temp(imin:imax,jmin:jmax,kmin:kmax)
    real(DP) :: xyz_Press(imin:imax,jmin:jmax,kmin:kmax)

    real(DP) :: xyz_Wsed(imin:imax,jmin:jmax,kmin:kmax)
    
    integer  :: i, j, k, s
    integer, parameter :: sw = 2  !ニュートン法
    logical            :: flag
    real(DP) :: x1
    real(DP) :: NuHh, NuVh
    real(DP), parameter :: AlphaNDiff = 1.0d-3  !数値拡散の係数
    
    !初期化
    xyz_Dn1Dt_prdt = 0.0d0
    xyz_Dn1Dt_loss = 0.0d0
    xyz_Dn1Dt_fall = 0.0d0
    xyz_Dn1Dt      = 0.0d0
    xyz_Dn2Dt_prdt = 0.0d0
    xyz_Dn2Dt_loss = 0.0d0
    xyz_Dn2Dt_fall = 0.0d0
    xyz_Dn2Dt      = 0.0d0

    xyz_Dn1DtAdv = 0.0d0
    xyz_Dn1DtDif = 0.0d0
    xyz_Dn1DtTrb = 0.0d0

    xyz_Dn2DtAdv = 0.0d0
    xyz_Dn2DtDif = 0.0d0
    xyz_Dn2DtTrb = 0.0d0

    xyzf_NDens1A = 0.0d0
    xyzf_NDens2A = 0.0d0

    
    if ( FlagCalc3D ) then 
      NuHh = AlphaNDiff * ( SQRT( dx * dy ) ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    else
      NuHh = AlphaNDiff * ( dx ** 4.0d0 ) / (2.0d0 * DelTimeLong)
    end if
    NuVh = AlphaNDiff * ( dz ** 4.0d0 ) / (2.0d0 * DelTimeLong)
        
    !温度・圧力に換算
    xyz_Temp  = ( xyz_PTempAl + xyz_PTempBZ ) * ( xyz_ExnerNl + xyz_ExnerBZ )
    xyz_Press = PressBasis * (( xyz_ExnerNl + xyz_ExnerBZ ) ** (CpDry / GasRDry))
    
    !移流
    xyz_Dn1DtAdv =                                                     &
      & - xyz_pyz( pyz_VelXNl * pyz_dx_xyz( xyzf_NDens1N(:,:,:,1) ) )  &
      & - xyz_xqz( xqz_VelYNl * xqz_dy_xyz( xyzf_NDens1N(:,:,:,1) ) )  &
      & - xyz_xyr( xyr_VelZNl * xyr_dz_xyz( xyzf_NDens1N(:,:,:,1) ) )    
    xyz_Dn2DtAdv =                                                     &
      & - xyz_pyz( pyz_VelXNl * pyz_dx_xyz( xyzf_NDens2N(:,:,:,1) ) )  &
      & - xyz_xqz( xqz_VelYNl * xqz_dy_xyz( xyzf_NDens2N(:,:,:,1) ) )  &
      & - xyz_xyr( xyr_VelZNl * xyr_dz_xyz( xyzf_NDens2N(:,:,:,1) ) )
    
    !数値拡散
    xyz_Dn1DtDif = &
      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyzf_NDens1B(:,:,:,1) ))))) &
      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyzf_NDens1B(:,:,:,1) ))))) &
      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyzf_NDens1B(:,:,:,1) )))))
    xyz_Dn2DtDif = &
      &  - NuHh * (xyz_dx_pyz(pyz_dx_xyz(xyz_dx_pyz(pyz_dx_xyz( xyzf_NDens2B(:,:,:,1) ))))) &
      &  - NuHh * (xyz_dy_xqz(xqz_dy_xyz(xyz_dy_xqz(xqz_dy_xyz( xyzf_NDens2B(:,:,:,1) ))))) &
      &  - NuVh * (xyz_dz_xyr(xyr_dz_xyz(xyz_dz_xyr(xyr_dz_xyz( xyzf_NDens2B(:,:,:,1) )))))
    
    !乱流拡散
    xyz_Dn1DtTrb =                                                                 &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyzf_NDens1B(:,:,:,1) ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyzf_NDens1B(:,:,:,1) ) )  &
      & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyzf_NDens1B(:,:,:,1) ) )
    xyz_Dn2DtTrb =                                                                 &
      &   xyz_dx_pyz( pyz_xyz( xyz_KhBl ) * pyz_dx_xyz( xyzf_NDens2B(:,:,:,1) ) )  &
      & + xyz_dy_xqz( xqz_xyz( xyz_KhBl ) * xqz_dy_xyz( xyzf_NDens2B(:,:,:,1) ) )  &
      & + xyz_dz_xyr( xyr_xyz( xyz_KhBl ) * xyr_dz_xyz( xyzf_NDens2B(:,:,:,1) ) )    

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          !落下速度
          call IH1998_Sediment( z_Z(k), xyz_Wsed(i,j,k) )
        end do
      end do
    end do
    
    
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          
          !落下項. 風上差分
          xyz_Dn1Dt_fall(i,j,k) = - (   xyz_Wsed(i,j,k+1) * xyzf_NDens1N(i,j,k+1,3) &
            &                         - xyz_Wsed(i,j,k)   * xyzf_NDens1N(i,j,k  ,3) ) / dz
          xyz_Dn2Dt_fall(i,j,k) = - (   xyz_Wsed(i,j,k+1) * xyzf_NDens2N(i,j,k+1,3) &
            &                         - xyz_Wsed(i,j,k)   * xyzf_NDens2N(i,j,k  ,3) ) / dz
          
          !H2SO4 生成・H2O 消滅
          call IH1998_H2SO4Prdt( z_Z(k), dz, xyz_Dn1Dt_prdt(i,j,k), xyz_Dn2Dt_loss(i,j,k) )
          
          !H2SO4 消滅・H2O 生成
          call IH1998_H2SO4Loss( xyz_Temp(i,j,k),       xyz_Press(i,j,k),       &
            &                         xyzf_NDens1N(i,j,k,2), xyzf_NDens2N(i,j,k,2),  &
            &                         xyz_Dn1Dt_loss(i,j,k), xyz_Dn2Dt_prdt(i,j,k) )
        end do
      end do
    end do

    xyz_Dn1Dt =   xyz_Dn1DtAdv   + xyz_Dn1DtDif   + xyz_Dn1DtTrb    &
      &         + xyz_Dn1Dt_fall + xyz_Dn1Dt_prdt + xyz_Dn1Dt_loss
    xyz_Dn2Dt =   xyz_Dn2DtAdv   + xyz_Dn2DtDif   + xyz_Dn2DtTrb    &
      &         + xyz_Dn2Dt_fall + xyz_Dn2Dt_prdt + xyz_Dn2Dt_loss    
              
    !時間積分
    xyzf_NDens1A(:,:,:,1) = &
      & max( xyzf_NDens1B(:,:,:,1) + 2.0d0 * DelTimeLong * xyz_Dn1Dt, 0.0d0 )
    xyzf_NDens2A(:,:,:,1) = &
      & max( xyzf_NDens2B(:,:,:,1) + 2.0d0 * DelTimeLong * xyz_Dn2Dt, 0.0d0 )

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          
          !凝結判定
          call IH1998_EquivState( xyz_Temp(i,j,k),            & !(IN)
            &   xyzf_NDens1A(i,j,k,1), xyzf_NDens2A(i,j,k,1), & !(IN)
            &   xyzf_NDens1A(i,j,k,2), xyzf_NDens2A(i,j,k,2), & !(OUT)
            &   xyzf_NDens1A(i,j,k,3), xyzf_NDens2A(i,j,k,3), & !(OUT)
            &   x1,                                           & !(OUT)
            &   sw,                                           & !(IN)
            &   flag                                          & !(OUT)
            & )
          
!          molfr1(i) = n1_gas(i) * Boltz * temp(i) / press(i)
!          molfr2(i) = n2_gas(i) * Boltz * temp(i) / press(i)

        end do
      end do
   end do
   
   !境界条件
   do s = 1, 3
     call SetMargin_xyz( xyzf_NDens1A(:,:,:,s) )
     call SetMargin_xyz( xyzf_NDens2A(:,:,:,s) )
   end do
   
 end subroutine Cloudphys_IH1998_forcing

end module cloudphys_IH1998
