!= Module DataSet
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: dataset.f90,v 1.2 2011/06/17 19:04:00 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]
!
!== Overview 
!
! ��ǥ�����Ѥ���ʪ��������Ū���������뤿����ѿ����ȷ��⥸�塼��
!
!== Error Handling
!
!== Known Bugs
!
!== Note
!
!== Future Plans
!

module dataset
  !
  ! ��ǥ�����Ѥ���ʪ��������Ū���������뤿����ѿ����ȷ��⥸�塼��
  !
  
  !�⥸�塼���ɤ߹���
  
  !���ۤη�����ػ�
  implicit none
  
  !save °��
  save
  
  !�����ѿ�
  
contains
  
  subroutine dataset_init(namelist_filename)
    !
    !=����
    !
    !NameList �ե����뤫�������������.
    !���Υ��֥롼�������, ���ؾ���ν������ԤäƤ���
    !
    !=�Ž���ʬ�μ�갷���ˤĤ���
    !
    !�׻������Ѥ���Ž���ʬ��, NAMELIST �ե�������� SpcWetSymbol ��
    !�񤫤�Ƥ���. �����ϰʲ����̤�. 
    !  
    !  H2O-g, H2O-l-Cloud, H2O-l-Rain, NH3-g, H2S-g, NH4SH-s-Cloud
    !
    !�ݥ���Ȥ�, ����ȸ���ˤ�, �� or ����ɽ��ʸ������ղä��뤳�ȤǤ���.
    !
    !���줾��β��ؼ�ˤ�, ChemData.f90 ��������� ID �ֹ�(ChemData_SpcID)��
    !�����Ƥ���. ���� ID ��, SpcWetID ���� 1 �������ݴɤ���. 
    !����ȶŽ�����б���, ChemData.f90 ��������� ChemData_SameSpc ��
    !���Ѥ���. �嵭����Ǥ�, SpcWetID �ϰʲ��Τ褦�ˤʤ�. 
    !
    ! Symbol:   H2O-g, H2O-l-Cloud, H2O-l-Rain, NH3-g, H2S-g, NH4SH-s-Cloud
    ! ID(:,1):      5,           6,          6,     8,     10,           11
    ! ID(:,2):      0,           5,          5,     0,      0,            8
    ! ID(:,3):      0,           0,          0,     0,      0,           10
    !
    !���Ѥ��ʤ���ʬ�ˤϥ�����������Ƥ���. 
    !
    
    
    !���ۤη�����ػ�
    implicit none

    

  end subroutine dataset_init
  
end module dataset
