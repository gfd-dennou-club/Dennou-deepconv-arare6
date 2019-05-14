module Imamura1998
  != ��������ʪ�׻��⥸�塼��
  !
  ! ����ʸ��Imamura and Hashimoto (1998), JGR103, E13, 31349. 
  !
  ! Authors::   SUGIYAMA Ko-ichiro, FUKUHARA Nozomu 
  !

  !���ۤη�����ػ�
  implicit none

  !°���λ���
  private

  !H2SO4 ���ϱդΥ��֥���ͳ���ͥ륮������뤿��η���.
  !Zeleznik (1991), J. Phys. Chem. Data �� Table 6 �˴�Ť�.  
  real(8)            :: mu111
  real(8)            :: mu211
  real(8)            :: mu121
  real(8)            :: mu221
  real(8)            :: mu122
  real(8)            :: mu212

  real(8)            :: eps111
  real(8)            :: eps211
  real(8)            :: eps121
  real(8)            :: eps221
  real(8)            :: eps122
  real(8)            :: eps212


  !�����ѿ�
  real(8), parameter :: Boltz = 1.38064852d-23  !!Boltzmann constant  
  real(8), parameter :: TempRef = 385.0d0
  real(8), parameter :: TempCr= 647.26d0

  !����
  public imamura1998_GibbsRDivRT
  public imamura1998_DelChemPotRDivRT
  public imamura1998_SatPress
  public imamura1998_SatPressRef
  public imamura1998_newton
  public imamura1998_bisection
  public imamura1998_Imamura1998
  public Imamura1998_Sediment
  public Imamura1998_H2SO4Prdt
  public Imamura1998_H2SO4Loss
  
contains

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  subroutine imamura1998_coefficient( temp )
    !
    !�׻������Ѥ��뷸�� \mu_{ijk}, \varepsilon_{ijk} ��׻�.
    !Zeleznik (1991), J. Phys. Chem. Data �� Table 6 �˴�Ť�.  
    !�����ϲ��٤δؿ�. 
    !
    
    ! ���ۤη�����ػ�
    implicit none
    
    ! �ѿ����
    real(8), intent(in) :: temp    !! ���� (temperature)
    
    !H2SO4 ���ϱդΥ��֥���ͳ���ͥ륮������뤿��η���.
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
  
    !���� \mu_ijk �η׻�
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
    
    !���� \varepsilon_ijk �η׻�
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
    
  end subroutine imamura1998_coefficient
  
  
  subroutine imamura1998_SatPressRef( TempIN, SatPress1, SatPress2 )
    !
    !˰�¾�����(ñ��ʬ)�η׻�.     
    !
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in)    :: TempIN  !! ����
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
    
    !���٤� TempCr �ʾ�ϼ��ʤ�.
    Temp = min(TempIN, TempCr)

    !!
    !! H2SO4 ��˰�¾����� 
    !!

    ! ���������
    xi   = ( ( Temp + 0.01d0 ) ** 2.0d0 ) - 293700.0d0
    zeta = (TempCr - Temp) ** 1.25d0  
    TempRefDivTemp = TempRef / Temp 
    
    ! H2SO4 (��ʪ��) ��˰�¾����� (ñ�� : bar )
    ! �ʤ��� log 2��ʤΤ�....
    !    Log_SatPressRef1 = - ( 10156.0d0 / Temp ) + 16.259d0
    Log_SatPressRef1 = &         
      & +  16.259d0 &
      & -  10156.0d0 / Temp  &
      & +  7.42d0 * ( 1.0d0 + log( TempRefDivTemp ) / log(2.0d0) - TempRefDivTemp )
    
    ! H2SO4 (�ϱ�) ��˰�¾����� (ñ�� : Pa)
    SatPress1 = exp( Log_SatPressRef1 ) * 1.0d5
    
    !!
    !! H20 ��˰�¾����� 
    !!
    
    ! H2O (��ʪ��) ��˰�¾����� (ñ�� : Pa)
    Log_SatPressRef2 = &
      &   alpha &
      & - beta / Temp &
      & + gamm * xi * ( dexp( delta * ( xi ** 2.0d0 ) ) - 1.0d0 ) / (Temp + 0.01d0) &
      & - eps * dexp( -1.0d0 * ita * zeta )
    
    ! H2O (�ϱ�) ��˰�¾�����
    SatPress2 = exp( Log_SatPressRef2 )
    
  end subroutine Imamura1998_SatPressRef
  

  subroutine imamura1998_SatPress( Temp, con, SatPress1, SatPress2 )
    !
    ! H2SO4 ���ϱդ�˰�¾�����
    !

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(8), intent(in)    :: con   !! H2SO4 ���ϱդ�ǻ��
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: Temp  !! ����
                                    !! temperature
    real(8), intent(out)   :: SatPress1
    real(8), intent(out)   :: SatPress2
    real(8)                :: SatPressRef1
    real(8)                :: SatPressRef2    
    real(8)                :: DelChemPotRDivRT1
    real(8)                :: DelChemPotRDivRT2 

    !��ʪ����˰�¾�����
    call imamura1998_SatPressRef( Temp, SatPressRef1, SatPressRef2 )
    
    !����β��إݥƥ󥷥���Ѳ���׻�
    call imamura1998_DelChemPotRDivRT( &
      & Temp, con, DelChemPotRDivRT1, DelChemPotRDivRT2 )

    ! H2SO4 (�ϱ�) ��˰�¾�����
    SatPress1 = SatPressRef1 * exp( DelChemPotRDivRT1 )
    
    ! H2O (�ϱ�) ��˰�¾�����
    SatPress2 = SatPressRef2 * exp( DelChemPotRDivRT2 )
    
  end subroutine Imamura1998_SatPress
  

  subroutine imamura1998_GibbsRDivRT( temp, con, GibbsRDivRT ) 
    !
    ! - G^{(r)} / RT �η׻�.
    ! G^{(r)} = G - G^{\circ}(T,p,x)
    ! (G^{(r)} : relative molar Gibbs energy of aquenous sulfuric asid)
    !
    
    !���ۤη�����ػ�
    implicit none

    !���������
    real(8), intent(in)    :: con   !! H2SO4 ���ϱդ�ǻ��
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: temp  !! ����
                                    !! temperature
    real(8), intent(out)   :: GibbsRDivRT

    
    !�ѿ������
    real(8)                :: x1, x2

    !����椬�ޥ��󥤥ץ�����꾮�������ν���.
    !\ln(0.0) ��ȯ������Τ�. 
    x1 = max( epsilon(con), con )
    x1 = min( 1.0d0 - epsilon(x1), x1 )

    ! H2SO4 �� H2O �Υ����
    x2 = 1.0d0 - x1

    !���������
    call imamura1998_coefficient(temp)   

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
    
  end subroutine imamura1998_GibbsRDivRT
  

  subroutine imamura1998_DelChemPotRDivRT( &
    & temp, con, DelChemPotRDivRT1, DelChemPotRDivRT2 ) 
    !
    ! \frac{\Del\mu_{i}}{RT} = \frac{\mu_i(T,x_i) - \mu^{\circ}(T)}{RT}
    ! = (G^{(r)} - x_1 \DP{G^{(r)}}{x_1} - x_2 \DP{G^{(r)}}{x_2} + \DP{G^{(r)}}{x_i})/RT
    !    
    
    !���ۤη�����ػ�
    implicit none

    !���������
    real(8), intent(in)    :: con   !! H2SO4 ���ϱդ�ǻ��
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: temp  !! ����
                                    !! temperature
    real(8), intent(out)   :: DelChemPotRDivRT1
    real(8), intent(out)   :: DelChemPotRDivRT2

    !�ѿ������ 
    real(8)                :: GibbsRDivRT
    real(8)                :: DGibbsRDivRT_Dx1
    real(8)                :: DGibbsRDivRT_Dx2
    real(8)                :: x1, x2
    
    !����椬�ޥ��󥤥ץ�����꾮�������ν���.
    !\ln(0.0) ��ȯ������Τ�. 
    x1 = max( epsilon(con), con )
    x1 = min( 1.0d0 - epsilon(x1), x1 )
    
    ! H2SO4 �� H2O �Υ����
    x2 = 1.0d0 - x1

    !���������
    call imamura1998_coefficient(temp)

    ! �ؿ��θƤӽФ�
    call imamura1998_GibbsRDivRT(temp, con, GibbsRDivRT)
   
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
    
  end subroutine imamura1998_DelChemPotRDivRT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

  subroutine imamura1998_DDelChemPotRDivRTDx( &
       & temp, con, DDelChemPotRDivRT1_Dx1, DDelChemPotRDivRT2_Dx1 )
    !
    ! \frac{\Del\mu_{i}}{RT} = \frac{\mu_i(T,x_i) - \mu^{\circ}(T)}{RT}
    ! = (G^{(r)} - x_1 \DP{G^{(r)}}{x_1} - x_2 \DP{G^{(r)}}{x_2} + \DP{G^{(r)}}{x_i})/RT
    ! x1 ��ʬ
    !
    
    !���ۤη�����ػ�
    implicit none

    !���������
    real(8), intent(in)    :: con   !! H2SO4 ���ϱդ�ǻ��
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: temp  !! ����
                                    !! temperature
    real(8), intent(out)   :: DDelChemPotRDivRT1_Dx1
    real(8), intent(out)   :: DDelChemPotRDivRT2_Dx1

    !�ѿ������
    real(8)                :: GibbsRDivRT
    real(8)                :: DGibbsRDivRT_Dx1Dx1
    real(8)                :: DGibbsRDivRT_Dx1Dx2
!    real(8)                :: DGibbsRDivRT_Dx2Dx2
    real(8)                :: x1, x2
    
    !����椬�ޥ��󥤥ץ�����꾮�������ν���.
    !\ln(0.0) ��ȯ������Τ�. 
    x1 = max( epsilon(con), con )
    x1 = min( 1.0d0 - epsilon(x1), x1 )
    
    ! H2SO4 �� H2O �Υ����
    x2 = 1.0d0 - x1

!    write(*,*) con
    
    !���������
    call imamura1998_coefficient(temp)

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
    
  end subroutine imamura1998_DDelChemPotRDivRTDx
  
  

  subroutine imamura1998_newton( NumDens1, NumDens2, Temp, con, SatPress1, SatPress2, Flag )

    !���ۤη����
    implicit none

    !�ѿ����
    real(8), intent(inout) :: con   !! H2SO4 ���ϱդ�ǻ��
                                    !! concentration of H2SO4 solution
    real(8), intent(in)    :: Temp  !! ����
                                    !! temperature
    real(8), intent(in)    :: NumDens1, NumDens2
    logical, intent(out)   :: Flag
    real(8)                :: Func, DFuncDx1
    real(8)                :: SatPress1, SatPress2
    real(8)                :: DSatPress1Dx1, DSatPress2Dx1
    real(8)                :: DDelChemPotRDivRT1_Dx1, DDelChemPotRDivRT2_Dx1
    real(8)                :: DelChemPotRDivRT1, DelChemPotRDivRT2
    real(8)                :: BoltzTemp
    real(8)                :: x1, x1A, x1P, x1L, x1R
    integer                :: CalNum
    integer, parameter     :: MaxNum = 10
    real(8)                :: x1save
    real(8)                :: gosa(10)
    
    !�����
    CalNum = 0
    x1save = con
    
    !���������
    BoltzTemp = Boltz * Temp
    
    !�ѿ�̾���֤�����
    x1 = con
    
    Loop: do
      !�롼�ײ��
      CalNum = CalNum + 1
      
      !���ϱդ�˰�¾�����
      call imamura1998_SatPress( Temp, x1, SatPress1, SatPress2 )
      
      !����β��إݥƥ󥷥��� x1 ��ʬ��׻�
      call imamura1998_DDelChemPotRDivRTDx( &
        & Temp, x1, DDelChemPotRDivRT1_Dx1, DDelChemPotRDivRT2_Dx1 )
      
      !H2SO4 (�ϱ�) ��˰�¾������� x1 (H2SO4 �����) ��ʬ
      DSatPress1Dx1 = DDelChemPotRDivRT1_Dx1 * SatPress1
      
      !H2O (�ϱ�) ��˰�¾������� x1 (H2SO4 �����) ��ʬ
      DSatPress2Dx1 = DDelChemPotRDivRT2_Dx1 * SatPress2       
      
      !�ؿ�
      Func = &
        &   ( 1.0d0 - x1 ) * (NumDens1 - SatPress1 / BoltzTemp ) &
        & - x1             * (NumDens2 - SatPress2 / BoltzTemp )
      
      !�ؿ��� x1 �˴ؤ��� 1 ����ʬ
      DFuncDx1 = &
        & - ( NumDens1 - SatPress1 / BoltzTemp ) &
        & - ( NumDens2 - SatPress2 / BoltzTemp ) &
        & - ( 1.0d0 - x1 ) * DSatPress1Dx1 / BoltzTemp  &
        & +   x1           * DSatPress2Dx1 / BoltzTemp         
      
      !���˶ᤤ��. �ͤ� 0 ~ 1 ���ϰ�. 
      x1A =  x1 - Func / DFuncDx1
      
      !�����ͤ� x1 ���������٤��ϰ� 0 <= x1 <= 1 ��Ķ�������ν���.
      if (x1A >= 1.0d0 ) then 
        !          write(*,*) "CALL Bisection (1)"
        x1L = x1
        x1R = 1.0d0
        call imamura1998_bisection  &
          & (NumDens1, NumDens2, Temp, x1L, x1R, con, SatPress1, SatPress2, Flag)
        exit Loop
      end if
!       if (x1A <= 0.0d0 ) then 
!          write(*,*) "CALL Bisection (2)"
!          x1L = 0.0d0
!          x1R = x1
!          call imamura1998_bisection(NumDens1, NumDens2, Temp, x1L, x1R, con, Flag)
!          exit Loop
!       end if

      !�롼�פν�λ���
      if ( abs((x1A - x1) / x1) < 1.0d-12 .OR. CalNum > MaxNum ) then 
        
        if( NumDens1 - SatPress1 / BoltzTemp > 0 .AND. NumDens2 - SatPress2 / BoltzTemp > 0) then
          !
          ! �ŷ�ξ�郎 OK �λ��ν���
          !
          con = x1A
          Flag = .true.
          
        elseif( abs(NumDens1 - SatPress1 / BoltzTemp) / NumDens1 < 1.0d-7 .OR. abs(NumDens2 - SatPress2 / BoltzTemp)/NumDens2 < 1.0d-7 .OR. CalNum > MaxNum .OR. x1A > 1.0d0 .OR. x1A < 0.0d0) then
          !
          ! ���٤���ݤ��뤿���, ��ʬˡ�Ǥ���ᤫǰ��������.
          !
          !             write(*,*) "CALL Bisection (3)"
          x1L = max(0.0d0, x1save - 1.0d-1)
          x1R = min(1.0d0, x1save + 1.0d-1)
          call imamura1998_bisection &
            & (NumDens1, NumDens2, Temp, x1L, x1R, con, SatPress1, SatPress2, Flag)
          
        else
          !
          ! �ŷ�ξ����������ʤ����ν���
          !
          con = -1.0d0
          Flag = .true.
        end if
        exit
      end if
      
      !�롼�פ�󤹤���ν���
      x1P = x1
      x1  = x1A
      
    end do Loop
    
  end subroutine Imamura1998_Newton
  
  
  subroutine imamura1998_bisection &
    & ( NumDens1, NumDens2, Temp, x1L, x1R, x1C, SatPress1, SatPress2, Flag )

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(8), intent(inout) :: x1L, x1R, x1C    !! H2SO4 �Υ���� 
    real(8), intent(in)    :: Temp  !! ����
                                    !! temperature
    real(8), intent(in)    :: NumDens1, NumDens2

    real(8), intent(out)   :: SatPress1, SatPress2
    logical, intent(out)   :: Flag
    real(8)                :: FuncL, FuncC, FuncR
    real(8)                :: BoltzTemp
    real(8)                :: Func
    integer                :: CalNum
    integer                :: i

    !�����
    CalNum = 0
    Flag = .false.
    
    !���������
    BoltzTemp = Boltz * Temp

    !�ʲ��ξ�郎Ʊ������������ʤ��Ȥ����ʤ�.
    ! NumDens1 - SatPress1 / BoltzTemp > 0
    ! NumDens2 - SatPress2 / BoltzTemp > 0
    
    Loop: do
      !�롼�ײ��
      CalNum = CalNum + 1
      
      !����. ���˶ᤤ��
      x1C = ( x1R + x1L ) * 5.0d-1
!       write(*,*) CalNum, real(x1C), real(abs((x1R - x1L) / x1R))

      !�ؿ�
      call imamura1998_SatPress( Temp, x1L, SatPress1, SatPress2 )
      FuncL = &
        &  ( 1.0d0 - x1L ) * max( (NumDens1 - SatPress1 / BoltzTemp ), 0.0d0) &
        & - x1L            * max( (NumDens2 - SatPress2 / BoltzTemp ), 0.0d0)
      
      call imamura1998_SatPress( Temp, x1R, SatPress1, SatPress2 )
      FuncR = &
        &  ( 1.0d0 - x1R ) * max( (NumDens1 - SatPress1 / BoltzTemp ), 0.0d0 ) &
        & - x1R            * max( (NumDens2 - SatPress2 / BoltzTemp ), 0.0d0 )
      
      call imamura1998_SatPress( Temp, x1C, SatPress1, SatPress2 )
      FuncC = &
        &  ( 1.0d0 - x1C ) * max( (NumDens1 - SatPress1 / BoltzTemp ), 0.0d0 ) &
        & - x1C            * max( (NumDens2 - SatPress2 / BoltzTemp ), 0.0d0 )
      
      !�롼�פ�󤹤���ν���
      if ( FuncL * FuncC < 0.0d0 ) then 
        x1R = x1C          
      elseif ( FuncC * FuncR < 0.0d0 ) then
        x1L = x1C
      end if
      
      !�롼�פν�λ���
      if ( abs((x1R - x1L) / x1R) < 1.0d-12 .OR. CalNum >= 60 .OR. FuncC == 0.0d0 ) then 
        
        if( NumDens1 - SatPress1 / BoltzTemp > 0 .AND. NumDens2 - SatPress2 / BoltzTemp > 0) then
          ! �ŷ�ξ����������Ƥ������ x1C ������ͤȤ���. 
          flag = .true.
          
        else
          ! �ŷ�ξ����������Ƥ��ʤ����� flag �� .false. �ˤ���.
          x1C = -1.0d0
          flag = .false.
        end if
        exit
      end if
      
       
    end do Loop
    
  end subroutine Imamura1998_Bisection
  

!!!
!!! �ޤȤ�Ʋ򤯤���Υץ����
!!!
  subroutine Imamura1998_Imamura1998( Temp, n1, n2, n1_gas, n2_gas, n1_liq, n2_liq, con, sw, flag )

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
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
    real(8)              :: del = 1.0d0 / real(NumTest - 1)
    real(8), parameter   :: TempCr= 647.26d0
    real(8)              :: x1N   
    real(8)              :: SatPressRef1, SatPressRef2
    real(8)              :: SatPress1, SatPress2
    real(8)              :: x1R, x1L

    !�����
    flag = .true. 
    
    !�ǡ����Τ�������դ���.      
    do j = 1, NumTest
      con = real(j-1) * del
      call imamura1998_SatPress( Temp, con, SatPress1, SatPress2)
      
      FuncTest2 = &
        &   ( 1.0d0 - con ) * max( (n1 - SatPress1 / ( boltz * Temp ) ), 0.0d0 ) &
        & -   con           * max( (n2 - SatPress2 / ( boltz * Temp ) ), 0.0d0 )
      
      if ( j==1 ) FuncTest1 = FuncTest2
      
      FuncSign(j)= FuncTest1 * FuncTest2
      FuncTest1  = FuncTest2        
    end do
    
    ! 0 <= x1 <= 1 ���ϰϤ˲򤬤���Τʤ�,
    ! dx �ֳ֤���Ʊ�Τγݤ��� (ex. Func(0) * Func(0.2)) ����ˤʤ�Τϰ�ս�.
    ! ��ˤʤ�ս������ź������ 1 �ʲ��ˤʤ�ʤ�. 
    !
    cnt = count( FuncSign < 0.0d0)
    idx = minloc(FuncSign, FuncSign < 0.0d0)
    
    ! �⤷�� 0 <= x1 <= 1 ���ϰϤ˲򤬤���Ȥ��������������ʤ��ʤ��,
    ! ���ϱդ���������ʤ�. ��ʪ����˰�¾��Τ߹�θ������ɤ�. 
    if ( (cnt(1) /= 1 .OR. idx(1) <= 1) .OR. n1 == 0.0d0 .OR. n2 == 0.0d0 ) then
      !
      !���ϱդ��ŷ뤷�ʤ����.
      !
      
      call imamura1998_SatPressRef( Temp, SatPressRef1, SatPressRef2)
      
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
      !���ϱդ��ŷ뤹����
      !
      
      !��ʬˡ�ν����
      x1L = real(idx(1)-2) * del
      x1R = real(idx(1)-1) * del
      
      if ( sw == 1 ) then 
        !��ʬˡ�Ƿ׻�
        call imamura1998_bisection(n1, n2, Temp, x1L, x1R, con, SatPress1, SatPress2, flag)
      else
        ! �˥塼�ȥ�ˡ�ν����
        con = (real(idx(1)-2) + real(idx(1)-1))* del * 5.0d-1
        call imamura1998_newton(n1, n2, Temp, con, SatPress1, SatPress2, flag)
      end if
      
      ! �ŷ��̤η׻�
!      call imamura1998_SatPress( Temp, con, SatPress1, SatPress2)
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
    
  end subroutine Imamura1998_Imamura1998
  

  subroutine Imamura1998_Sediment(Height, VelZ)

    !���ۤη�����ػ�
    implicit none

    !�ѿ����
    real(8), intent(in)  :: Height
    real(8), intent(out) :: VelZ
    real(8), parameter   :: Grav = 8.7d0
    real(8), parameter   :: Dens = 1.8d3
    real(8), parameter   :: R_2  = 1.15d-6  !�⡼��2 ��Ⱦ��
    real(8), parameter   :: R_3  = 3.65d-6
    real(8), parameter   :: Vis  = 1.5d-5
    real(8)              :: Molfr_2, Molfr_3
    real(8)              :: VelZ_2,  VelZ_3
    
    !mode3 ��ʬ�ۤ� tanh �Ƕ��
    Molfr_3 = 0.455d0 * tanh( 57.5d3 - Height ) + 0.455d0
    
    !mode2 ��ʬ��
    Molfr_2 = 1.d0 - Molfr_3
    
    !mode 2, 3 ��®��
    VelZ_3  = - 2.0d0 * Grav * Dens * R_3 * R_3 / Vis / 9.0d0
    VelZ_2  = - 2.0d0 * Grav * Dens * R_2 * R_2 / Vis / 9.0d0
    
    !ʿ��®��
    VelZ = Molfr_2 * VelZ_2 + Molfr_3 * VelZ_3
    
  end subroutine Imamura1998_Sediment
  
  
  subroutine Imamura1998_H2SO4Prdt( Height, DelZ, H2SO4_Prdt, H2O_Loss )
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
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
    
    !ñ�̤�·���뤿��� DelZ �ǳ�껻���Ƥ���ɬ�פ�����Ϥ�.
    !
    H2SO4_Prdt = & 
      & 3.5d0  * 1.0d16  / sigma / sqrt( 2.0d0 * PI ) &
      & * exp( -1.0d0 * DH * DH  / ( 2.0d0 * sigma * sigma ) ) &
      & * cos( Phi ) / DelZ
    
    H2O_Loss = -1.0d0 * H2SO4_Prdt
    
  end subroutine Imamura1998_H2SO4Prdt
  
  
  subroutine Imamura1998_H2SO4Loss( Temp, Press, n1, n2, H2SO4_Loss, H2O_Prdt )
    
    !���ۤη�����ػ�
    implicit none
    
    !�ѿ����
    real(8), intent(in) :: Temp
    real(8), intent(in) :: Press
    real(8), intent(in) :: n1, n2      !��̩��
    real(8), intent(out):: H2SO4_Loss
    real(8), intent(out):: H2O_Prdt
    
    real(8), parameter  :: k13 = 5.5d-29
    real(8), parameter  :: f_co = 25.0d-6  !���� > 35
!     real(8), parameter  :: f_co = 30.0d-6  !35 < ���� < 60
!     real(8), parameter  :: f_co = 35.0d-6  !���� > 60
    real(8)              :: mean
    real(8)              :: kp
    
    !����ѿ�
    mean = press / (Boltz * Temp)
    kp = 10.0d0 ** ( 7.584d0 - ( 5060.0d0 / Temp ) )
    
    !H2SO4 �ξ��ǹ�
    H2SO4_Loss = f_co * k13 * kp * mean * mean  * n1 / n2 / press
    
    !H2O ��������
    H2O_Prdt = -1.0d0 * H2SO4_Loss
    
  end subroutine Imamura1998_H2SO4Loss
  
end module imamura1998
