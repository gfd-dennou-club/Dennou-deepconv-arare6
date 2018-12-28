#= Totorial of deepconv/arare: rakuraku deepconv
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: namelist_init.rd,v 1.1 2014/03/04 04:44:05 sugiyama Exp $
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


== ����������

���ܾ� (FlagBasic) �Ⱦ���� (FlagDisturb...) �����ꤹ��. 
���ܾ�������ɬ��. ������ 1 �İʾ����ꤹ��ɬ�פ�����. 

  &initialdata_nml
    !
    ! [ɬ��] ���ܾ������. �ʲ��������򤹤뤳��. 
    !
    ! FlagBasic = "IsoThermal"   ! �׻��ΰ�β��ٰ���
    ! FlagBasic = "Dry"          ! ������Ǯ�ʴ��ܾ�
    ! FlagBasic = "Sounding"     ! ������ǥ��󥰥ե�������Ϳ����
    ! FlagBasic = "Yamasaki1983" ! Yamasaki et al (1983) ���Ϥ������ܾ�
    ! FlagBasic = "Takemi2007"   ! Takemi (2007) ���Ϥ������ܾ�
    ! FlagBasic = "Toon2002"     ! Toon et al. (2002) ���Ϥ������ܾ�
    ! FlagBasic = "Baker1998"    ! deepconv/arare4 �ν��Ϥ�����
    ! FlagBasic = "Arare4"       ! deepconv/arare4 �ν��Ϥ�����
    ! FlagBasic = "Arare4mmc"    ! deepconv/arare4 (�����׻�) �ν��Ϥ�����
    ! 
    ! ���̤ξ����. �ʲ��������ǽ. 
    !
    ! FlagDisturbPTemp = "GaussXY"   ! ���������� (XY ʿ��, Z ��������)
    ! FlagDisturbPTemp = "GaussXZ"   ! ���������� (XZ ʿ��, Y ��������)
    ! FlagDisturbPTemp = "GaussYZ"   ! ���������� (YZ ʿ��, X ��������)
    ! FlagDisturbPTemp = "GaussXYZ"  ! ���������� (XY ʿ��, 3 ����Ū��ʬ��)
    ! FlagDisturbPTemp = "Random"    ! ������
    ! FlagDisturbPTemp = "Rectangle" ! ���
    ! FlagDisturbPTemp = "CosXY"     ! cos (XY ʿ��, Z ��������)
    ! FlagDisturbPTemp = "CosXZ"     ! cos (XZ ʿ��, Y ��������)
    ! FlagDisturbPTemp = "CosYZ"     ! cos (YZ ʿ��, X ��������)
    ! FlagDisturbPTemp = "CosXYZ"    ! cos 
    ! FlagDisturbPTemp = "ConeXY"    ! �߿�
    ! FlagDisturbPTemp = "ConeXZ"    ! �߿�
    ! FlagDisturbPTemp = "ConeYZ"    ! �߿�
    ! FlagDisturbPTemp = "tanh"      ! tanh ���Υ���
    !
    ! �������ʡ��ؿ��ξ����. �ʲ��������ǽ. 
    ! 
    ! FlagDisturbExner = "GaussXY"   ! ���������� (XY ʿ��, Z ��������)
    ! FlagDisturbExner = "GaussXZ"   ! ���������� (XZ ʿ��, Y ��������)
    ! FlagDisturbExner = "GaussYZ"   ! ���������� (YZ ʿ��, X ��������)
    ! FlagDisturbExner = "GaussXYZ"  ! ���������� 
    ! FlagDisturbExner = "CosXY"     ! cos (XY ʿ��, Z ��������)
    ! FlagDisturbExner = "CosXZ"     ! cos (XZ ʿ��, Y ��������)
    ! FlagDisturbExner = "CosYZ"     ! cos (YZ ʿ��, X ��������)
    ! FlagDisturbExner = "CosXYZ"    ! cos
    ! FlagDisturbExner = "ConeXY"    ! �߿�
    ! FlagDisturbExner = "ConeXZ"    ! �߿�
    ! FlagDisturbExner = "ConeYZ"    ! �߿�
    !
    ! �ŷ���ʬ�κ�����. �ʲ��������ǽ. 
    !
    ! FlagDisturbQMix = "GaussXY"    ! ���������� (XY ʿ��, Z ��������)
    ! FlagDisturbQMix = "GaussXZ"    ! ���������� (XZ ʿ��, Y ��������)
    ! FlagDisturbQMix = "GaussYZ"    ! ���������� (YZ ʿ��, X ��������)
    ! FlagDisturbQMix = "GaussXYZ"   ! ����������
    ! FlagDisturbQMix = "Dryreg"     ! �����ΰ������
#    ! FlagDisturbQMix = "Moist"      ! ?
    ! FlagDisturbQMix = "CosXY"      ! cos (XY ʿ��, Z ��������)
    ! FlagDisturbQMix = "CosXZ"      ! cos (XZ ʿ��, Y ��������)
    ! FlagDisturbQMix = "CosYZ"      ! cos (YZ ʿ��, X ��������)
    ! FlagDisturbQMix = "CosXYZ"     ! cos
    ! FlagDisturbQMix = "ConeXY"     ! �߿�
    ! FlagDisturbQMix = "ConeXZ"     ! �߿�
    ! FlagDisturbQMix = "ConeYZ"     ! �߿�
    !
    ! ®�پ�κ�����. �ʲ��������ǽ. 
    !
    ! FlagDisturbWind = "Takemi2007" ! Takemi (2007) ���Ϥ���®�پ�
    ! FlagDisturbWind = "Sounding"   ! ������ǥ��󥰥ե����뤫��®�پ��Ϳ����
    ! FlagDisturbWind = "Zonal"      ! ��ʿ���ͤ�®�پ�
    ! FlagDisturbWind = "tanh"       ! tanh ���Υ���
    !
    ! FlagDisturb = "Arare4"         ! deepconv/arare4 �ν��Ϥ�����
    ! FlagDisturb = "Arare4mmc"      ! deepconv/arare4 (�����׻�) �ν��Ϥ�����
  /

== ���ܾ������κݤ�ɬ�פȤʤ����

=== FlagBasic = "IsoThermal" �ξ��

�äˤʤ�. 

=== FlagBasic = "Dry" �ξ��

  &initialdata_basic_nml
     Humidity
     TempTr
     DHeight
     HeightTr
  /

=== FlagBasic = "Sounding"  �ξ��

  &initialdata_sounding_nml
     SoundingFile = ""     ! ������ǥ��󥰥ե�����
     AltCol       = 0      ! �ֹ��١פ����ֹ� (������ǥ��󥰥ե�������)
     TempCol      = 0      ! �ֲ��١פ����ֹ� (������ǥ��󥰥ե�������)
     PressCol     = 0      ! �ְ��ϡפ����ֹ� (������ǥ��󥰥ե�������)
     VelXCol      = 0      ! ��X ������®�١פ����ֹ� (������ǥ��󥰥ե�������)
     VelYCol      = 0      ! ��Y ������®�١פ����ֹ� (������ǥ��󥰥ե�������)
     !
     ! �ʲ��ξ����, ���ꤷ�����٤�������ٰ�������ط������ꤹ�����ɬ��.
     !
     AltTr        = 0      ! ��ή�����̤ι���
     DelAlt       = 4.0d3  ! ���ٸ�Ψ�򥼥�˴��¤���ޤǤε�Υ�����
  /

=== FlagBasic = "Yamasaki1983" 

�ä����ꤹ��ɬ�פϤʤ�. 
����Ū���ͤ�, ((<����ͤη׻���|URL:http://www.gfd-dennou.org/library/deepconv/sample/2011-06-23_sugiyama/Init/exp.htm>)) ����.


=== FlagBasic = "Takemi2007" 

���ܾ����ꤹ��. 
����Ū���ͤ�, ((<����ͤη׻���|URL:http://www.gfd-dennou.org/library/deepconv/sample/2011-06-23_sugiyama/Init/exp.htm>)) ����.

  &initialdata_takemi2007_nml
    ! 
    ! ���ܾ������
    !
    ! FlagEnv = "MidLat_Q10"
    ! FlagEnv = "MidLat_Q12"
    ! FlagEnv = "MidLat_Q14"
    ! FlagEnv = "MidLat_Q16"
    ! FlagEnv = "MidLat_Q16DRY1"
    ! FlagEnv = "MidLat_Q16DRY2"
    ! FlagEnv = "MidLat_Q18"
    ! FlagEnv = "Tropic_Q18"
    ! FlagEnv = "Tropic_Q18DRY1"
    ! FlagEnv = "Tropic_Q18DRY2"
    ! FlagEnv = "Tropic_Q18DRY3"
    !
    !! ®�پ� (��ľ�ץ�ե�����) ������
    !!
    !! FlagWind = "LowLevel"
    !! FlagWind = "MiddleLevel"
    !! FlagWind = "HighLevel"
    !!
    !! ��ɽ���ն��®��
    !!
    !! VelXSfc = 0.0d0
  /

=== FlagBasic = "Toon2002"

�ä����ꤹ��ɬ�פϤʤ�. 
����Ū���ͤ�, ((<����ͤη׻���|URL:http://www.gfd-dennou.org/library/deepconv/sample/2011-06-23_sugiyama/Init/exp.htm>)) ����.

=== FlagBasic = "Baker1998"

�ä����ꤹ��ɬ�פϤʤ�. 

=== FlagBasic = "Arare4" or "Arare4mmc" �ξ��

deepconv/arare4 �ν��Ϥ򸵤˴��ܾ졦�����������ˤϰʲ��������Ԥ�. 

  &arare4fileio_nml
    Arare4Prefix = ""   ! deepconv/arare4 �Υե��������Ƭ������
  /


== �����������ɬ�פȤʤ���� 

=== ���������� [����, �������ʡ��ؿ�, ������]

FlagDisturbPTemp, FlagDisturbExner, FlagDisturbQMix �ʤɤ�, 
(({"Gauss..."})) �����򤷤�����, 
initialdata_disturb_gauss_nml �������Ԥ�. 

  &initialdata_disturb_gauss_nml
     PTempMax = 1.0d0      ! ���̾���κ�����
     ExnerMax = 0.01d0     ! �������ʡ��ؿ��κ�����
     QMixMax  = 1.0d-2     ! ������ (����) �κ�����
     Xc       = 1000.0d0   ! ��ɸ���濴 (X ����)
     Xr       = 100.0d0    ! Ⱦ�� (X ����)
     Yc       = 1000.0d0   ! ��ɸ���濴 (Y ����)
     Yr       = 100.0d0    ! Ⱦ�� (Y ����)
     Zc       = 1000.0d0   ! ��ɸ���濴 (Z ����)
     Zr       = 100.0d0    ! Ⱦ�� (Z ����)
  /

��: ɬ�פ��ͤΤ����ꤹ����ɤ�. �㤨��, (({FlagDisturbPTemp = "GaussXY"})) ��
���ꤷ�����ˤ�, �ʲ��Τ褦�����ꤹ��н�ʬ�Ǥ���. 

  &initialdata_disturb_gauss_nml
     PTempMax = 1.0d0      ! ���̾���κ�����
     Xc       = 1000.0d0   ! ��ɸ���濴 (X ����)
     Xr       = 100.0d0    ! Ⱦ�� (X ����)
     Yc       = 1000.0d0   ! ��ɸ���濴 (Y ����)
     Yr       = 100.0d0    ! Ⱦ�� (Y ����)
  /


=== ���Ѵؿ� [����, �������ʡ��ؿ�, ������]

FlagDisturbPTemp, FlagDisturbExner, FlagDisturbQMix �ʤɤ�, 
(({"Cos..."})) �����򤷤�����, 
initialdata_disturb_cos_nml �������Ԥ�. 

  &initialdata_disturb_cos_nml
     PTempMax = 1.0d0      ! ���̾���κ�����
     ExnerMax = 0.01d0     ! �������ʡ��ؿ��κ�����
     QMixMax  = 1.0d-2     ! ������ (����) �κ�����
     Xc       = 1000.0d0   ! ��ɸ���濴 (X ����)
     Xr       = 100.0d0    ! Ⱦ�� (X ����)
     Yc       = 1000.0d0   ! ��ɸ���濴 (Y ����)
     Yr       = 100.0d0    ! Ⱦ�� (Y ����)
     Zc       = 1000.0d0   ! ��ɸ���濴 (Z ����)
     Zr       = 100.0d0    ! Ⱦ�� (Z ����)
  /

��: ɬ�פ��ͤΤ����ꤹ����ɤ�. �㤨��, (({FlagDisturbPTemp = "CosXY"})) ��
���ꤷ�����ˤ�, �ʲ��Τ褦�����ꤹ��н�ʬ�Ǥ���. 

  &initialdata_disturb_cos_nml
     PTempMax = 1.0d0      ! ���̾���κ�����
     Xc       = 1000.0d0   ! ��ɸ���濴 (X ����)
     Xr       = 100.0d0    ! Ⱦ�� (X ����)
     Yc       = 1000.0d0   ! ��ɸ���濴 (Y ����)
     Yr       = 100.0d0    ! Ⱦ�� (Y ����)
  /


=== �߿� [����, �������ʡ��ؿ�, ������]

FlagDisturbPTemp, FlagDisturbExner, FlagDisturbQMix �ʤɤ�, 
(({"Cone..."})) �����򤷤�����, 
initialdata_disturb_cone_nml �������Ԥ�. 

  &initialdata_disturb_cone_nml
     PTempMax = 1.0d0      ! ���̾���κ�����
     ExnerMax = 0.01d0     ! �������ʡ��ؿ��κ�����
     QMixMax  = 1.0d-2     ! ������ (����) �κ�����
     Xc       = 1000.0d0   ! ��ɸ���濴 (X ����)
     Xr       = 100.0d0    ! Ⱦ�� (X ����)
     Yc       = 1000.0d0   ! ��ɸ���濴 (Y ����)
     Yr       = 100.0d0    ! Ⱦ�� (Y ����)
     Zc       = 1000.0d0   ! ��ɸ���濴 (Z ����)
     Zr       = 100.0d0    ! Ⱦ�� (Z ����)
  /

��: ɬ�פ��ͤΤ����ꤹ����ɤ�. �㤨��, (({FlagDisturbPTemp = "ConeXY"})) ��
���ꤷ�����ˤ�, �ʲ��Τ褦�����ꤹ��н�ʬ�Ǥ���. 

  &initialdata_disturb_cone_nml
     PTempMax = 1.0d0      ! ���̾���κ�����
     Xc       = 1000.0d0   ! ��ɸ���濴 (X ����)
     Xr       = 100.0d0    ! Ⱦ�� (X ����)
     Yc       = 1000.0d0   ! ��ɸ���濴 (Y ����)
     Yr       = 100.0d0    ! Ⱦ�� (Y ����)
  /

=== ������ [����]

�Ȥ�����٤˥�����ʾ�������ꤹ���� (FlagDisturbPTemp = "Random"), 
initialdata_disturb_random_nml �������Ԥ�. 

  &initialdata_disturb_random_nml
     PTempMax = 1.0d0      ! ����κ�����
     Zpos     = 1.0d3      ! ������֤�����
  /

=== ��� [����], DryReg [������]

����ΰ�Τ߲��٤��Ѳ��������� (FlagDisturbPTemp = "Rectangle"), 
initialdata_disturb_rectangle_nml �������Ԥ�. 

  &initialdata_disturb_rectangle_nml
     PTempMax = 1.0d0      ! ���̤�����ʬ
     XposMin  = 0.0d0      ! X �����λ���
     XposMax  = 1.0d2      ! X �����ν���
     YposMin  = 0.0d0      ! Y �����λ���
     YposMax  = 1.0d2      ! Y �����ν���
     ZposMin  = 0.0d0      ! Z �����λ���
     ZposMax  = 1.0d2      ! Z �����ν���
  /

����ΰ�ζŷ�ʪ�κ�����򥼥�ˤ����� (FlagDisturbQMix = "DryReg") ��, 
initialdata_disturb_rectangle_nml �������Ԥ�. 

  &initialdata_disturb_rectangle_nml
     XposMin  = 0.0d0      ! X �����λ���
     XposMax  = 1.0d2      ! X �����ν���
     YposMin  = 0.0d0      ! Y �����λ���
     YposMax  = 1.0d2      ! Y �����ν���
     ZposMin  = 0.0d0      ! Z �����λ���
     ZposMax  = 1.0d2      ! Z �����ν���
  /

=== tanh [����, ®�پ�] 

����ӥ󡦥إ��ۥ���԰���ν���ͤ��Ѥ���褦�� tanh ���Υ�����Ϳ������ 
(FlagDisturbPTemp = "tanh" & FlagDisturbWind = "tanh"), 
initialdata_disturb_tanh_nml �����ꤹ��. 

  &initialdata_disturb_tanh_nml
     PTempMean  = ʿ�Ѳ���
     PTempDel   = �����ؤξ岼�β��̺�     
     VelMean    = ʿ����®
     Zc         = �����ؤ��濴�ι���
     Zr         = �����ؤθ���
  /

=== Zonal [®�پ�]

���ͤ���®���Ϳ������ (FlagDisturbWind = "Zonal"), 
initialdata_disturb_wind_zonal_nml �����ꤹ��. 

  &initialdata_disturb_wind_zonal_nml/ &
    VelX0 = 10.0d0
    VelY0 = 10.0d0
    VelZ0 = 10.0d0
  /

=== Sounding [®�پ�]

  &initialdata_sounding_nml
     SoundingFile = ""     ! ������ǥ��󥰥ե�����
     VelXCol      = 0      ! ��X ������®�١פ����ֹ� (������ǥ��󥰥ե�������)
     VelYCol      = 0      ! ��Y ������®�١פ����ֹ� (������ǥ��󥰥ե�������)
  /

=== Takemi2007 [®�پ�]

®�پ�����ꤹ��. 

  &initialdata_takemi2007_nml
    !! 
    !! ���ܾ������
    !!
    !! FlagEnv = "MidLat_Q10"
    !! FlagEnv = "MidLat_Q12"
    !! FlagEnv = "MidLat_Q14"
    !! FlagEnv = "MidLat_Q16"
    !! FlagEnv = "MidLat_Q16DRY1"
    !! FlagEnv = "MidLat_Q16DRY2"
    !! FlagEnv = "MidLat_Q18"
    !! FlagEnv = "Tropic_Q18"
    !! FlagEnv = "Tropic_Q18DRY1"
    !! FlagEnv = "Tropic_Q18DRY2"
    !! FlagEnv = "Tropic_Q18DRY3"
    !!
    ! ®�پ� (��ľ�ץ�ե�����) ������
    !
    ! FlagWind = "LowLevel"
    ! FlagWind = "MiddleLevel"
    ! FlagWind = "HighLevel"
    !
    ! ��ɽ���ն��®��
    !
    ! VelXSfc = 0.0d0
  /

=== "Arare4" or "Arare4mmc" �ξ��

deepconv/arare4 �ν��Ϥ򸵤˴��ܾ졦�����������ˤϰʲ��������Ԥ�. 

  &arare4fileio_nml
    Arare4Prefix = ""   ! deepconv/arare4 �Υե��������Ƭ������
  /

=end
