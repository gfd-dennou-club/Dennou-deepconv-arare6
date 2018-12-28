#= Totorial of deepconv/arare: rakuraku deepconv
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: namelist.rd,v 1.4 2014/03/04 04:44:05 sugiyama Exp $
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

== ����¾������

�ʲ��Ǥ�, ����� (���ܾ졦�����) ��ʪ������������ʳ���������ܤ��������. 

=== �����ϥե�����˴ؤ�������

�ե�����˽񤭹������. 

  &fileset_nml
    !�ե������ǽ�Ū���ѹ�������, �ȿ�. 
    FileInstitution = 'GFD Dennou Club (http://www.gfd-dennou.org)'
    ! �¸�̾ 
    FileTitle       = 'cloud moist convection experiment' (�ǥե������)
    ! �ǡ����ե���������ץ����̾
    FileSource      =  'deepconv/arare5 (http://www.gfd-dennou.org/library/deepconv)' (�ǥե������)
   /

�ǥե�����ͤΤޤޤ��ɤ������ѿ����ɬ�פϤʤ�. �㤨��, �¸�̾��
�ץ����̾�ϥǥե���ȤΤޤޤ��ɤ����ˤ�, FileInstitution �Τ�
���ꤹ����ɤ�. 

  &fileset_nml
    !�ե������ǽ�Ū���ѹ�������, �ȿ�. 
    FileInstitution = 'XXXX@yyy.org'
   /

=== �ꥹ�����ȥե�����˴ؤ�������

�ꥹ�����ȤΤ�������ϡ����ϥե��������ꤹ��. 

  &restartfileio_nml
    InputFile  = "arare_restart.nc",   ! ���ϥե�����̾
    OutputFile = "arare_restart2.nc",  ! ���ϥե�����̾
  /

=== ��ʬ���֤˴ؤ�������

���ֹ����ʬ����, �ꥹ�����Ȼ��֤����ꤹ��. 
�ꥹ�����Ȥ�����ν�����, ((<������|URL:./settings5.htm>)) �򻲾Ȳ�����. 

  &timeset_nml
    DelTimeLong   = 4.0d0      !Ĺ�������ॹ�ƥå� [sec]
    DelTimeShort  = 4.0d-1     !û�������ॹ�ƥå�(���ȴ�Ϣ��) [sec]
    RestartTime   = 0.0d0      !���ϻ��� [sec]
    IntegPeriod   = 4800.0d0   !��ʬ���� [sec]
    DelTimeOutput = 240.0d0    !�ꥹ�����ȥե�����˽񤭽Ф����ֹ�� [sec]
  /

=== �ʻ������˴ؤ�������

�ʻ����˴ؤ���������ꤷ�ޤ�. 
2 �����׻���Ԥ����ν�����, ((<������|URL:./calc_2d.htm>)) �򻲾Ȳ�����. 
MPI ���Ѥ�������׻��Τ������, ((<������|URL:./calc_mpi.htm>)) �򻲾Ȳ�����. 

  &gridset_nml
    xsub  = 1                 ! X ����������� (�ǥե����)
    ysub  = 1                 ! Y ����������� (�ǥե����)
    xdim  = 10                ! X ����������� (�ǥե����)
    ydim  = 10                ! Y ����������� (�ǥե����)
    zdim  = 10                ! Z ����������� (�ǥե����)
    NCMAX = 1                 ! �ŷ���ʬ�ο� (�ǥե����)
    Xmg   = 2                 ! X �����Υޡ����� (�ǥե����)
    Ymg   = 2                 ! Y �����Υޡ����� (�ǥե����)
    Zmg   = 2                 ! Z �����Υޡ����� (�ǥե����)
  /

=== ��ɸ���˴ؤ�������

��ɸ���˴ؤ���������ꤷ�ޤ�. 

  &axesset_nml
    Xmin  = 0.0d0              ! X ��ɸ�λ��� (�ǥե����)
    Xmax  = 1.0d4              ! X ��ɸ�ν��� (�ǥե����)
    Ymin  = 0.0d0              ! X ��ɸ�λ��� (�ǥե����)
    Ymax  = 1.0d4              ! X ��ɸ�ν��� (�ǥե����)
    Zmin  = 0.0d0              ! Z ��ɸ�λ��� (�ǥե����)
    Zmax  = 1.0d4              ! Z ��ɸ�ν��� (�ǥե����)
  /

�ǥե�����ͤϻ��ꤷ�ʤ��Ƥ��ɤ��Ǥ�. 

  &axesset_nml
    Xmax  = 1.0d4              ! X ��ɸ�ν��� (�ǥե����)
    Ymax  = 1.0d4              ! X ��ɸ�ν��� (�ǥե����)
    Zmax  = 1.0d4              ! Z ��ɸ�ν��� (�ǥե����)
  /

=== �¸��ѥ�᥿(ʪ��������)������

�갵��Ǯ��������Ǯ��ʬ���̡������������ꤹ����ˡ�ϰʲ��� 3 �̤ꤢ��ޤ�. 
�ܤ���������ˡ�ˤĤ��Ƥ� ((<������|URL:./settings3.htm>)) �򻲾Ȳ�����. 

* �갵��Ǯ��ʬ���̤��ۤ�Ϳ������
* �갵��Ǯ�ȵ��������Ϳ������
* ʪ��̾�Ȥ��Υ�����Ϳ������

�ǥե�����ͤϰʲ����̤�Ǥ�. 

  &constants_nml
    Grav = 9.8d0          !���� [m/s^2]
    PressBasis = 965.0d0  !���̤δ�వ�� [Pa]
    TempSfc = 0.0d0       !��ɽ�̲��� [K]
    PressSfc = 0.0d0      !��ɽ�̰��� [Pa]
    TempTop = 0.0d0       !���������β��� [K]
    PressTop = 0.0d0      !���������Ǥΰ��� [Pa]
    CpDry  = 0.0d0        !������ʬ���갵��Ǯ [J/K kg]
    CpDryMol = 0.0d0      !������ʬ���갵��Ǯ [J/K kg]
    CvDry = 0.0d0         !������ʬ��������Ǯ [J/K kg]
    MolWtDry = 0.0d0      !������ʬ��ʬ����   [kg/mol]
    GasRDry  = 0.0d0      !������ʬ�ε������ [J/K kg]
    DayTime = 86400.0d0   ! 1 ����Ĺ�� [s]
  /


=== �����˴ؤ�������

ɬ�פ˱�����Ŭ������������. �ʲ����ϵ��絤�ξ��Ǥ���.

  &composition_nml
    SpcWetSymbol(1)  = 'H2O-g',  !������ʬ
    SpcWetSymbol(2)  = 'H2O-s-Cloud', !������ʬ
    SpcWetSymbol(3)  = 'H2O-s-Rain',  !������ʬ
    SpcWetMolFr(1)   = 1.0d-2,  !������ʬ��¸����(�����)
    SpcWetMolFr(2)   = 0.0d0,  !������ʬ��¸����(�����)
    SpcWetMolFr(3)   = 0.0d0,   !������ʬ��¸����(�����)
  /


=== ���그��������(���ݥ���) 

���ݥ��ؤθ��그������������ꤹ��. DepthVb �϶����¸��Τ褦��, 
���������˥��ݥ��ؤ����ꤹ��������Ѥ���. 

  &damping_nml
    EFTime = 3.0d2             !���ݥ��ؤθ��그���� e-folding time
    DepthH = 0.0d0             !���ݥ��ؤθ���(��ʿ����)
    DepthV = 0.0d0             !���ݥ��ؤθ���(��ľ����) [��������]
    DepthVb= 0.0d0             !���ݥ��ؤθ���(��ľ����) [��������] 
  /

=== ���그��������(���ȸ����)

�̾���ѹ�����ɬ�פϤʤ�.

  &dynamics_nml
   AlphaSound      = 2.0e-7    !���ȸ����η���
   AlphaNDiff      = 1.0d-4    !���ͳȻ���η���
   NDiffRatio      = 1.0d0     !®�٤��Ф���Ǵ����夲����Ͽ����� 1 �ʾ�ˤ���. 
   beta            = 1.0       !��ľ������������ˡ�ˤ�����ˤ� 1.0
                               !����󥯥˥��륽��ˡ�ˤ������ 0.5 �Ȥ���. 
   FactorBuoyTemp  = 1.0d0     !���� (���٤δ�Ϳ) ��̵ͭ
                               !��θ���ʤ������ͤ򥼥�ˤ���.
   FactorBuoyMolWt = 1.0d0     !���� (ʬ���̸���) ��̵ͭ
                               !��θ���ʤ������ͤ򥼥�ˤ���.
   FactorBuoyLoading = 1.0d0   !���� (�ٽŸ���) ��̵ͭ
                               !��θ���ʤ������ͤ򥼥�ˤ���.
  /



=== ���Ϥ���������

�ܺ٤�((<������|URL:settings4.htm>)) ����������.

  &gtool_historyauto_nml
    FilePrefix = "thermal-moist_" ! ���ϥե��������Ƭ��
    IntValue = 10.0,              ! ���ϴֳ֤ο���
    IntUnit = 'sec',              ! ���ϴֳ֤�ñ��
  /

=== �ǡ������Ϥθ�������

�����ѿ���θ��̤������Ԥ�. �ܺ٤�((<������|URL:settings4.htm>)) ����������.

  &gtool_historyauto_nml
    Name = ''                     ! �����ѿ�̾
    IntValue = 10.0,              ! ���ϴֳ֤ο���
    IntUnit = 'sec',              ! ���ϴֳ֤�ñ��
  !  TimeAverage = .true.,        ! ����ʿ�Ѥ�̵ͭ
  !  SpaceAverage = .true.,       ! ��ʿʿ�Ѥ�̵ͭ
  /


=end
