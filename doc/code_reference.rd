#= Makefile for deepconv/arare reference mannual.
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: code_reference.rd,v 1.2 2011/12/19 08:05:04 odakker Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
########################################################################
#

=begin JA

= deepconv/arare �����������ɥ�ե����

#* ����
#  * 2011/12/19(��������) ����
#  * 2008/06/19(����ã��) ����
#  * 2007/10/19(��������) ����
#  * 2006/11/10(��������) ����
#  * 2006/11/09(��������) ��������

=end JA
=begin EN

= Deepconv/arare source code reference

#* History
#  * 2008/06/19 (Tatsuya Yamashita) Update
#  * 2007/10/19 (Masatsugu Odaka) Update
#  * 2006/11/10 (Masatsugu Odaka) Update
#  * 2006/11/09 (Masatsugu Odaka) Initial release

=end EN


=begin JA

== �¹ԥץ����ΰ���

((<��ץ����(main)|URL:../src/main>)):
* ((<arare.f90|URL:code_reference/files/__/src/main/arare_f90.html>)): 
  ��ץ����
* ((<arare_init-data.f90|URL:code_reference/files/__/src/main/arare_init-data_f90.html>)): 
  ����ͷ׻��ѥץ����



=end JA
=begin EN

== Execute programs

((<Main program (main)|URL:../src/main>)):
* ((<arare.f90|URL:code_reference/files/main/arare_f90.html>)): 
  Main program
* ((<arare_init-data.f90|URL:code_reference/files/__/src/main/arare_init-data_f90.html>)): 
  Main program for initial data calculation.


=end EN


=begin JA
== ���֥롼����ȥ⥸�塼��ΰ���

((<���ز���(chemdata)|URL:../src/chemdata>))
* ((<ChemData|URL:code_reference/classes/ChemData.html>)):
  ���ؼ�ǡ����ݴɥ⥸�塼��


((<�ϳز���(dynamic)|URL:../src/dynamics>))
* ((<DynamicsHEVI|URL:code_reference/classes/DynamicsHEVI.html>)):
  �ϳز����׻��ѥ⥸�塼��


((<���������(env)|URL:../src/env/>))
* ((<initialdata_basic|URL:code_reference/classes/initialdata_basic.html>)):
  ���ܾ������ѥ⥸�塼��
* ((<initialdata_disturb|URL:code_reference/classes/initialdata_disturb.html>)):
  ������������ѥ⥸�塼��

* ((<initialdata_yamasaki1983|URL:code_reference/classes/initialdata_yamasaki1983.html>)):
  Yamasaki (1983) �ν���������ѥ⥸�塼��
* ((<initialdata_Toon2002|URL:code_reference/classes/initialdata_Toon2002.html>)):
  Toon (2002) �ν���������ѥ⥸�塼��
* ((<initialdata_takemi2007|URL:code_reference/classes/initialdata_takemi2007.html>)):
  Takemi (2007) �ν���������ѥ⥸�塼��


((<������(io)|URL:../src/io/>))
* ((<ReStartFileIO|URL:code_reference/classes/ReStartFileIO.html>)):
  �ꥹ�����ȥե�����ؤν��ϥ⥸�塼��
* ((<HistoryFileIO|URL:code_reference/classes/HistoryFileIO.html>)):
  �ҥ��ȥ�ե�����ؤν��ϥ⥸�塼��
* ((<Arare4fileio|URL:code_reference/classes/Arare4fileio.html>)):
  deepconv/arare4 �Ǻ��������ꥹ�����ȥե���������ϥ⥸�塼��


((<ʪ������(physics)|URL:../src/physics/>))
* ((<cloudphys_k1969|URL:code_reference/classes/cloudphys_k1969.html>)):
  Kessler (1969) �ˤ�����ʪ���ѥ�᥿�ꥼ�������׻��ѥ⥸�塼��
* ((<Cloudphys_MarsCond|URL:code_reference/classes/Cloudphys_MarsCond.html>)):
  �����絤����ʪ���ѥ�᥿�ꥼ�������׻��ѥ⥸�塼��
* ((<ECCM|URL:code_reference/classes/ECCM.html>)):
  �徺������Ǯ�����β��ٸ�Ψ��, ���ξ���ʿ���絤��¤��׻�����⥸�塼��
* ((<MoistAdjust|URL:code_reference/classes/MoistAdjust.html>)):
  ����˰��Ĵ��ˡ�⥸�塼��
* ((<Radiation_HeatBalance|URL:code_reference/classes/Radiation_HeatBalance.html>)):
  ���ͤ��Ϥ���Ǯ�����׻��ѥ⥸�塼��
* ((<Radiation_Simple|URL:code_reference/classes/Radiation_Simple.html>)):
  ���ͤ��Ϥ���Ǯ�����׻��ѥ⥸�塼��(ñ��С������)
* ((<Surfaceflux_bulk|URL:code_reference/classes/Surfaceflux_bulk.html>)):
  �Х륯ˡ�ˤ����ɽ�ե�å����׻��ѥ⥸�塼��
* ((<Surfaceflux_diff|URL:code_reference/classes/Surfaceflux_diff.html>)):
  �Ȼ��ˤ����ɽ�ե�å����׻��ѥ⥸�塼��
* ((<Turbulence_kw1978|URL:code_reference/classes/Turbulence_kw1978.html>)):
  Klemp & Wilhelmson (1978) ����ή�ѥ�᥿�ꥼ�������׻��ѥ⥸�塼��


((<�������(setup)|URL:../src/setup/>))
* ((<argset|URL:code_reference/classes/argset.html>)):
  ���ޥ�ɥ饤���������ѥ⥸�塼��
* ((<axesset|URL:code_reference/classes/axesset.html>)):
  3 �������ֳָ�߳ʻҳʻ�������⥸�塼��
* ((<basicset|URL:code_reference/classes/basicset.html>)):
  ���ܾ�����⥸�塼��
* ((<ChemCalc|URL:code_reference/classes/ChemCalc.html>)):
  ���ش�Ϣ�̷׻��⥸�塼��
* ((<clockset|URL:code_reference/classes/clockset.html>)):
  �׻����־�������ѥ⥸�塼��
* ((<composition|URL:code_reference/classes/composition.html>)):
  �����������⥸�塼��
* ((<constats|URL:code_reference/classes/constants.html>)):
  ��������ѥ⥸�塼��  
* ((<constats0|URL:code_reference/classes/constants0.html>)):
  ʪ���������������⥸�塼��
* ((<dataset|URL:code_reference/classes/dataset.html>)):
  ʪ��������Ū�ѥ�᡼������⥸�塼�� 
* ((<fileset|URL:code_reference/classes/fileset.html>)):
  �����ϥե�����̾����⥸�塼��
* ((<gridset|URL:code_reference/classes/gridset.html>)):
  �ʻ������󥵥�������⥸�塼��
* ((<mpi_wrapper|URL:code_reference/classes/mpi_wrapper.html>)):
  MPI ��åѡ��⥸�塼��
* ((<namelist_util|URL:code_reference/classes/namelist_util.html>)):
  NAMELIST �ե��������Ϥ˴ؤ���⥸�塼��
* ((<timset|URL:code_reference/classes/timeset.html>)):
  ������ʬ�ѥѥ�᡼������⥸�塼��


((<�������⥸�塼��(util)|URL:../src/util/>))
* ((<CFLCheck|URL:code_reference/classes/CFLCheck.html>)):
  CFL ����ǧ�⥸�塼��
* ((<Damping|URL:code_reference/classes/Damping.html>)):
  ���ȸ����ȥ��ݥ��ؤǤ��໤��η׻��⥸�塼��
* ((<FillNegative|URL:code_reference/classes/FillNegative.html>)):
  �����̤ʤɤ��������̤η����׻��⥸�塼��
* ((<setmargin|URL:code_reference/classes/setmargin.html>)):
  �����ΰ����󥵥�������⥸�塼��
* ((<TimeFilter|URL:code_reference/classes/TimeFilter.html>)):
  ���֥ե��륿���׻��⥸�塼��
* ((<xyz_bc_module|URL:code_reference/classes/xyz_bc_module.html>)):
  �����������⥸�塼��
* ((<xyz_deriv_c4_module|URL:code_reference/classes/xyz_deriv_c4_module.html>))
  4 �������濴��ʬ�׻��⥸�塼��
* ((<xys_deriv_module|URL:code_reference/classes/xyz_deriv_module.html>)):
  2 �������濴��ʬ�׻��⥸�塼��


=end JA
=begin EN

== Subroutines and modules

((<Chemical process (chemdata)|URL:../src/chemdata>))
* ((<ChemData|URL:code_reference/classes/ChemData.html>)):
  Chemical data module


((<Dynamics (dynamic)|URL:../src/dynamics>))
* ((<DynamicsHEVI|URL:code_reference/classes/DynamicsHEVI.html>)):
  Module for dynamical process

((<Initial environment setup (env)|URL:../src/env/>))
* ((<initialdata_basic|URL:code_reference/classes/initialdata_basic.html>)):
  Basic state set up module
* ((<initialdata_disturb|URL:code_reference/classes/initialdata_disturb.html>)):
  Initial disturbance set up module

* ((<initialdata_yamasaki1983|URL:code_reference/classes/initialdata_yamasaki1983.html>)):
  Initial value used by Yamasaki (1983) setup module
* ((<initialdata_Toon2002|URL:code_reference/classes/initialdata_Toon2002.html>)):
  Initial value used by Toon et al. (2002) setup module
* ((<initialdata_takemi2007|URL:code_reference/classes/initialdata_takemi2007.html>)):
  Initial value used by Takemi (2007) setup module


((<"Input/Output (io)"|URL:../src/io/>))
* ((<HistoryFileIO|URL:code_reference/classes/HistoryFileIO.html>)):
  I/O module of history files
* ((<ReStartFileIO|URL:code_reference/classes/ReStartFileIO.html>)):
  I/O module of restart file
* ((<Arare4fileio|URL:code_reference/classes/Arare4fileio.html>)):
  I/O module of restart file generated by deepconv/arare4


((<Physics (physics)|URL:../src/physics/>))
* ((<cloudphys_k1969|URL:code_reference/classes/cloudphys_k1969.html>)):
  Kessler (1969) cloud parameterization module
* ((<Cloudphys_MarsCond|URL:code_reference/classes/Cloudphys_MarsCond.html>)):
  Cloud parameterization module for the Martian atmosphere
* ((<ECCM|URL:code_reference/classes/ECCM.html>)):
  ECCM (Ensemble Cloud Condensation Model) module
* ((<MoistAdjust|URL:code_reference/classes/MoistAdjust.html>)):
  Moist adjustment module
* ((<Radiation_HeatBalance|URL:code_reference/classes/Radiation_HeatBalance.html>)):
  Thermal forcing module associated with atmospheric radiation
* ((<Radiation_Simple|URL:code_reference/classes/Radiation_Simple.html>)):
  Thermal forcing module associated with atmospheric radiation (simple verion)
* ((<Surfaceflux_bulk|URL:code_reference/classes/Surfaceflux_bulk.html>)):
  Surface flux calculation module by using bulk method
* ((<Surfaceflux_diff|URL:code_reference/classes/Surfaceflux_diff.html>)):
  Surface flux calculation module by using diffusion
* ((<Turbulence_kw1978|URL:code_reference/classes/Turbulence_kw1978.html>)):
  Klemp & Wilhelmson (1978) turbulent parameterization module

((<Set up (setup)|URL:../src/setup/>))
* ((<argset|URL:code_reference/classes/argset.html>)):
  Command line argument module
* ((<axesset|URL:code_reference/classes/axesset.html>)):
  Grid arrangement set up module
* ((<basicset|URL:code_reference/classes/basicset.html>)):
  Basic state set up module
* ((<ChemCalc|URL:code_reference/classes/ChemCalc.html>)):
  Chemical process module
* ((<clockset|URL:code_reference/classes/clockset.html>)):
  Clock information set up module
* ((<composition|URL:code_reference/classes/composition.html>)):
  Chemical constants set up module
* ((<constats|URL:code_reference/classes/constants.html>)):
  Constant values set up module
* ((<constats0|URL:code_reference/classes/constants0.html>)):
  Physical and mathmatical constants set up module
* ((<dataset|URL:code_reference/classes/dataset.html>)):
  Physical and chemical parameter set up module
* ((<fileset|URL:code_reference/classes/fileset.html>)):
  I/O file names set up module
* ((<gridset|URL:code_reference/classes/gridset.html>)):
  Grid array size set up module
* ((<mpi_wrapper|URL:code_reference/classes/mpi_wrapper.html>)):
  MPI wrappaer module
* ((<namelist_util|URL:code_reference/classes/namelist_util.html>)):
  NAMELIST file name parameter set up module
* ((<timset|URL:code_reference/classes/timeset.html>)):
  Time integration parameters set up module


((<Utility (util)|URL:../src/util/>))
* ((<CFLCheck|URL:code_reference/classes/CFLCheck.html>)):
  CFL condition check 
* ((<Damping|URL:code_reference/classes/Damping.html>)):
  Sound wave damping term and Rayleigh damping near the upper boundary
* ((<FillNegative|URL:code_reference/classes/FillNegative.html>)):
  Fulfill negative value of positive definite variables
* ((<setmargin|URL:code_reference/classes/setmargin.html>)):
  Array size of margine area set up module
* ((<TimeFilter|URL:code_reference/classes/TimeFilter.html>)):
  Time filter for time integration
* ((<xyz_bc_module|URL:code_reference/classes/xyz_bc_module.html>)):
  Adapting boundary condition 
* ((<xyz_deriv_c4_module|URL:code_reference/classes/xyz_deriv_c4_module.html>))
  4th order centered differentiate scheme
* ((<xys_deriv_module|URL:code_reference/classes/xyz_deriv_module.html>)):
  2'nd order centered differentiate scheme 

=end EN
