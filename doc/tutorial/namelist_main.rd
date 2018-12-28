#= Totorial of deepconv/arare: rakuraku deepconv
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: namelist_main.rd,v 1.1 2014/03/04 04:44:05 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
########################################################################
#
=begin TOPLINK
[((<GFD Dennou Club|URL:http://www.gfd-dennou.org>)) |
((<deepconv Project|URL:http://www.gfd-dennou.org/library/deepconv>))]
[((<deepconv Reference Manual|URL:../index.htm>))]
=end TOPLINK

=begin

= ����ե�����γ���


== ʪ������������

deepconv_main_nml �ˤ�����, ���ͼ¸������Ѥ���ʪ�����������򤹤�. 

  &deepconv_main_nml
    !
    ! ��ή����
    !FlagTurbMethod  = "KW1978",
    !  KW1978 (Klemp and Wilhelmson 1978)
    !  ConstKm (��ή�Ȼ���������)
    !
    ! ����ʪ������
    !FlagCloudMethod = "K1969",
    !   K1969 (�Ȥ������Υѥ�᥿�ꥼ������� (Kessler, 1969))
    !   MarsCond (�Ȼ���Ĺ)
    !
    ! ���Ͳ���
    !FlagRadMethod = "HeatConst"
    !  HeatConst (���ꤷ�����٤�������)
    !  HeatVary  (���ꤷ�����٤�������, ���Ψ�Ϲ��٤δؿ�)
    !  HeatBalance (�Ϥ�Ϳ�����Ǯ��Ѥ����礦�褦��Ǯ������Ϳ����)
    !  Baker1998 (������, Baker et al., 1998)
    !  Sounding  (�ե����뤫���ǮΨ�����Ψ��Ϳ����)
    !
    ! ��ɽ�̲���
    !FlagSurfaceMethod = "Diff"
    !  Diff (��������Ȼ�Ū��Ǯ��ʪ����Ϳ����, ���������β��١�ʪ���̸���)
    !  Bulk (�Х륯ˡ)
    !  Const (����, ������������ΰ����Ǯ����ư�̡�ʪ���ե�å�����Ϳ����)
    !  Baker1998 (������, Baker et al., 1998)
    !
    ! �ǥХå���
    !FlagDebugMethod = "Const"
    !  WindConst (®�ٰ���, Ϳ����줿®�٤Ƿ׻���̤��񤭤���.)
    !  NoTendencyLong (Ĺ�����֥��ƥåפǷ׻����줿 tendency �򥼥�Ȥ���.)
    !    
   /    


== ��ή����

=== FlagTurbMethod  = "KW1978",

Klemp and Wilhelmson (1978) �� 1.5 ���Υ������㡼�˴�Ť�����ή������׻�������ˤ�
�ʲ����ͤ����ꤹ�뤳�Ȥ������. �ǥե�����ͤΤޤޤ��ɤ�����, �ͤ����ꤷ�ʤ����ɤ�. 

  &turbulence_kw1978_nml
     Cm     = 2.0d-1          !��ή���ͥ륮�����Ǽ��η���
     KmMax  = 800.0d0         !��ή�Ȼ������κ�����
     FlagDExnerDtTurb =.true. !��������������ή�Ȼ�����θ���뤫�Υ����å�
  /

=== FlagTurbMethod  = "ConstKm"

��ή�Ȼ��������ͤ�����ͤˤ�����ˤ�, �ʲ����ͤ����ꤹ�뤳�Ȥ������. 
�ǥե�����ͤΤޤޤ��ɤ�����, �ͤ����ꤷ�ʤ����ɤ�. 

  &turbulence_constKm_nml
     Cm     = 2.0d-1          !��ή���ͥ륮�����Ǽ��η���
     MixLen = 0.0d0           !ʿ�Ѻ����Υ
     ConstKm = 0.0d0          !��ư�̤��Ф�����ή�Ȼ�����
     ConstKh = 0.0d0          !Ǯ���Ф�����ή�Ȼ�����
     FlagDispHeat = .false.   !�����Ǯ���θ���뤫�Υ����å�
     FlagDExnerDtTurb =.true. !��������������ή�Ȼ�����θ���뤫�Υ����å�
  /


== ��ʪ���ѥ�᥿�ꥼ������������

=== FlagCloudMethod = "K1969"

�Ȥ������Υѥ�᥿�ꥼ������� (Kessler, 1969) ���Ѥ����������. 

  &cloudphys_k1969_nml
    Planet            = ""        ! "Earth" of "Jupiter"
                                  !   Earth:   FactorJ = 1.0 �����ꤵ���
                                  !   Jupiter: FactorJ = 3.0 �����ꤵ���
    FactorJ           = 1.0d0     ! ��ʪ�������Υѥ�᡼��
                                  !   �����Ǥ� 3.0d0
                                  !   �ϵ�Ǥ� 1.0d0 �Ȥ���
    AutoConvTime      = 1000.0d0  ! ʻ����Ĺ�λ���� [sec]
    QMixCr            = 1.0d-3    ! ʻ����Ĺ���������׳������� [kg/kg]
    FlagDExnerDtCloud = .true.    !���Ϥμ��˶ŷ�θ��̤��θ���뤫�ݤ�
                                  !��θ���ʤ������ͤ� .false. �ˤ���.
    FlagDExnerDtFall  = .true.    !���Ϥμ�����θ��̤��θ���뤫�ݤ�
                                  !��θ���ʤ������ͤ� .false. �ˤ���.
    FactorFallRain    = 1.0d0     !�������̵ͭ
                                  !��θ���ʤ������ͤ򥼥�ˤ���.
    FactorCloud2Rain  = 1.0d0     !�����鱫�ؤ��Ѵ���̵ͭ
                                  !��θ���ʤ������ͤ򥼥�ˤ���.
    FactorRain2Gas    = 1.0d0     !����������ؤ��Ѵ���̵ͭ
                                  !��θ���ʤ������ͤ򥼥�ˤ���.
    FactorCloud2Gas   = 1.0d0     !����������ؤ��Ѵ���̵ͭ 
                                  !��θ���ʤ������ͤ򥼥�ˤ���.
  /

=== FlagCloudMethod = "MarsCond"

Yamasita et al (��ƽ�����) ���Ѥ��Ƥ���Ȼ���Ĺ

  &cloudphys_marscond_nml
    DensIce     = 1.565d3  ! �����̩�� [kg/m^3]
    NumAerosol  = 0.0d0    ! ��������ο�̩�� [1/kg]
    RadiAerosol = 0.0d0    ! ��������ο�̩�� [1/kg]
    Kd          = 0.0d0    ! �絤��Ǯ��Ƴ���� [W/K m]
    SatRatioCr  = 0.0d0    ! �׳�˰���� 
    SatRtWetAdia = 0.0d0   ! ������Ǯ����˰���� 
    CO2LatHeat  = 0.0d0    ! ñ�̼��̤�����ζŷ�Ǯ [J/kg]
    CDensCr     = 5.0d-5   ! ����
  /

== ���Ͳ���

=== FlagRadMethod = "HeatConst"

��ñ����: �Ȥ�����٤���Ͳ�Ǯ�����.

  &radiation_simple_nml
    RadHeatRate= 0.0d0,     !�������Ͷ������礭�� [K/day]
    HeightUp   = 10.0d3,    !���Ͷ�����Ϳ�����ľ�ΰ�ξ��
    HeightDown = 0.0d3,     !���Ͷ�����Ϳ�����ľ�ΰ�β��� 
    FlagDExnerDtRad = .true. !���������������ͤˤ���Ǯ��Ѥδ�Ϳ���θ���뤫�Υ����å�
  /

=== FlagRadMethod = "HeatVary"

��ñ����: ��ɽ�̤��� HeightDown �ޤǤ� RadHeatRate �����. 
HeightUp ������ϲ�ǮΨ����ˤʤ�褦�˲�ǮΨ�򸺾�������. 

  &radiation_simple_nml
    RadHeatRate= 0.0d0,      !���Ͷ������礭�� [K/day]
    HeightUp   = 10.0d3,     !���Ͷ�����Ϳ�����ľ�ΰ�ξ��
    HeightDown = 0.0d3,      !���Ͷ�����Ϳ�����ľ�ΰ�β��� 
    FlagDExnerDtRad = .true. !���������������ͤˤ���Ǯ��Ѥδ�Ϳ���θ���뤫�Υ����å�
  /

=== FlagRadMethod = "HeatBalance"

��ñ���ͥ⥸�塼��: �Ȥ�������ΰ�������ѡ���Ǯ����. 
���Ψ������ե������Ϳ����. ��ǮΨ����Ѥ����礦�褦�˥�ǥ������Ƿ���. 

  &radiation_heatbalance_nml
    RadCoolRate = 0.0d0      !�������Ͳ�ǮΨ [K/day]                 
    HeightHeatUp   = 0.0d0   !��Ǯ�ΰ�ξ�ü�ι���
    HeightHeatDown = 0.0d0   !��Ǯ�ΰ�β�ü�ι���
    HeightCoolUp   = 0.0d0   !����ΰ�ξ�ü�ι���
    HeightCoolDown = 0.0d0   !����ΰ�β�ü�ι���
    FlagDExnerDtRad = .true. !���������������ͤˤ���Ǯ��Ѥδ�Ϳ���θ���뤫�Υ����å�
  /

=== FlagRadMethod = "Baker1998"

Baker et al. (1998) �ζ����׻��Ǥ����Ͳ���. ���ꤹ����ܤ�̵��. 

=== FlagRadMethod = "Sounding"

  &radiation_sounding_nml
    SoundingFile    = ""       ! ������ǥ��󥰥ե�����
    AltCol          = 0        !�ֹ��١פ����ֹ� (������ǥ��󥰥ե�������)
    SWaveCol        = 0        !��û�����ͤˤ���ǮΨ�פ����ֹ� (������ǥ��󥰥ե�������)   
    LWaveCol        = 0        !��Ĺ�����ͤˤ���ǮΨ�פ����ֹ� (������ǥ��󥰥ե�������)   
    FlagDExnerDtRad = .true.   !���������������ͤˤ���Ǯ��Ѥδ�Ϳ���θ���뤫�Υ����å�
  /

== ��ɽ�̲���
  
=== FlagSurfaceMethod = "Diff"  

���������γȻ��������ᤦ������.   

  &surfaceflux_diff_nml
     Kappa = 800.0d0               ! ���������Ǥ���ή�Ȼ�����
     FlagDExnerDtSurf = .true.     ! Flag for diabatice heating term in pressure equation
  /


===  FlagSurfaceMethod = "Bulk"  

  &surfaceflux_bulk_nml
    FlagConstBulkCoef
                            ! Flag for using constant bulk coefficient
    FlagUseOfBulkCoefInNeutralCond
                            ! Flag for using bulk coefficient in neutral condition
    FlagDExnerDtSurf = .true.  
                            ! Flag for diabatice heating term in pressure equation
    ConstBulkCoef           
                            ! �Х륯����������. 
                            ! Steady value of bulk coefficient
    VelMinForRi = 1.0d-8    ! ����㡼�ɿ��׻���®�ٲ�����
                            ! Lower limit of velocity for Ri
    SfcRoughLength = 1.0d-2 ! ����Ĺ��
                            ! Roughness length
    Vel0 = 0.0d0            ! ���ؤǤο�ʿ®�ٿ�夲��
                            ! 
    VelBulkCoefMin = 0.0d0  ! $ u $ �Х륯�����Ǿ���. 
                            ! Minimum value of $ u $ bulk coefficient
    TempBulkCoefMin = 0.0d0 ! $ T $ �Х륯�����Ǿ���. 
                            ! Minimum value of $ T $ bulk coefficient
    QmixBulkCoefMin = 0.0d0 ! $ q $ �Х륯�����Ǿ���. 
                            ! Minimum value of $ q $ bulk coefficient
    VelBulkCoefMax = 1.0d2  ! $ u $ �Х륯����������. 
                            ! Maximum value of $ u $ bulk coefficient
    TempBulkCoefMax = 1.0d2 ! $ T $ �Х륯����������. 
                            ! Maximum value of $ T $ bulk coefficient
    QmixBulkCoefMax = 1.0d2 ! $ q $ �Х륯����������. 
                            ! Maximum value of $ q $ bulk coefficient
  /

===  FlagSurfaceMethod = "Const"

����, ������������ΰ����Ǯ����ư�̡�ʪ���ե�å�����Ϳ�������ˤ�
�ʲ��ι��ܤ����ꤹ��. 

  &surfaceflux_const_nml
     SfcXMomFluxBtm = 0.0d0    ! X �����α�ư�̥ե�å��� (��������)
     SfcXMomFluxTop = 0.0d0    ! X �����α�ư�̥ե�å��� (��������)
     SfcYMomFluxBtm = 0.0d0    ! Y �����α�ư�̥ե�å��� (��������)
     SfcYMomFluxTop = 0.0d0    ! Y �����α�ư�̥ե�å��� (��������)
     SfcHeatFluxBtm = 0.0d0    ! Ǯ�ե�å��� (��������)
     SfcHeatFluxTop = 0.0d0    ! Ǯ�ե�å��� (��������)
     SfcQmixFluxBtm = 0.0d0    ! ʪ���ե�å��� (��������)
     SfcQmixFluxTop = 0.0d0    ! ʪ���ե�å��� (��������)
     FlagDExnerDtSurf = .true. ! Flag for diabatice heating term in pressure equation
  /

===  FlagSurfaceMethod = "Baker1998"

Baker et al. (1998) �ζ����׻��Ǥ���ɽ��(?)����. ���ꤹ����ܤ�̵��. 


=end
