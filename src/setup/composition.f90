!= �ŷ���ʬ�˴ؤ����������뤿����ѿ����ȷ��⥸�塼��
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: composition.f90,v 1.7 2014/07/08 01:05:32 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module composition
  !
  != �ŷ���ʬ�˴ؤ����������뤿����ѿ����ȷ��⥸�塼��
  !
  
  !�⥸�塼���ɤ߹���
  !
  use dc_types,      only: DP

  !���ۤη�����ػ�
  !
  implicit none
  
  ! �ǥե���Ȥ�°��
  !
  private
  
  ! �ѿ����
  !
  integer, save, public       :: GasNum       = 0   ! ����ο�
  integer, save, public       :: CloudNum     = 0   ! ���ο�  
  integer, save, public       :: RainNum      = 0   ! ���ο�
  integer, save, public       :: IdxG(10)     = 0   ! ���Τ�����ź����
  integer, save, public       :: IdxC(10)     = 0   ! ��������ź����
  integer, save, public       :: IdxR(10)     = 0   ! ��������ź����  
  integer, save, public       :: CondNum      = 0   ! �ŷ�����ο�
  integer, save, public       :: IdxCG(10)    = 0   ! �ŷ����(����)������ź����
  integer, save, public       :: IdxCC(10)    = 0   ! �ŷ����(��)������ź����
  integer, save, public       :: IdxCR(10)    = 0   ! �ŷ����(��)������ź����
  integer, save, public       :: RactNum      = 0   ! ����ȿ���ο�
  integer, save, public       :: IdxNH3       = 0   ! NH3 (����)������ź����
  integer, save, public       :: IdxH2S       = 0   ! H2S (����)������ź����
  integer, save, public       :: IdxNH4SHc    = 0   ! NH4SH (��)������ź����
  integer, save, public       :: IdxNH4SHr    = 0   ! NH4SH (��)������ź����

  real(DP), save              :: SpcWetMolFr(20)    !������ʬ�β��ؼ��¸����
  character(20), save         :: SpcWetSymbol(20)   !������ʬ�β��ؼ�̾  

  integer, allocatable, save  :: SpcWetID(:)        !������ʬ�β��ؼ��ID
  real(DP), allocatable, save :: MolWtWet(:)        !������ʬ��ʬ����  

  integer, allocatable, save  :: IDGas(:)           !�����β��ؼ��ID
  integer, allocatable, save  :: IDCloud(:)         !���β��ؼ��ID
  integer, allocatable, save  :: IDRain(:)          !���β��ؼ��ID

  ! ���֥롼����θ���
  !
  public SpcWetID, MolWtWet, SpcWetMolFr, SpcWetSymbol
  public IDGas, IDCloud, IDRain
  public composition_init
  
contains
  
  subroutine composition_init
    !
    !=����
    !
    !NameList �ե����뤫�������������.
    !���Υ��֥롼�������, ���ؾ���ν������ԤäƤ���
    !
    !=�Ž���ʬ�μ�갷���ˤĤ���
    !
    !�׻������Ѥ���Ž���ʬ�ξ���� basicset.f90 ����������
    !SpcWetSymbol �� SpcWetID ���ݴɤ���Ƥ���
    ! 
    ! Symbol:  H2O-g, NH3-g, H2S-g, H2O-l-Cloud, H2O-l-Rain, NH4SH-s-Cloud, NH4SH-s-Rain
    ! ID:      5,     8,     10,    7,           7,          11,            11
    !
    !ID �ֹ�(ChemData_SpcID)�� ChemData.f90 ��������Ƥ���
    !
    !�嵭�ξ���򸵤�, ���Υ롼����Ǥϰʲ��ξ������.  
    !
    !  * �ƥ��ƥ��꡼�˴ޤޤ��ʪ���ο�
    !
    !    GasNum = 3,  CloudNum = 2, RainNum = 2
    !  
    !  * �ƥ��ƥ��꡼������ź����. ���Τ�����������������������Ѥ���.
    !
    !    IdxG = 1, 2, 3, 0, 0, 0, ...
    !    IdxC = 4, 6, 0, 0, 0, 0, ...
    !    IdxR = 5, 7, 0, 0, 0, 0, ...
    !
    !  * �ŷ�(Condensation)��������ʪ���ο���, ����������ź����. 
    !    �嵭����Ǥ� H2O �ζŷ�Τߤ�������
    !
    !    CondNum = 1
    !    IdxCG = 1, 0, 0, 0, 0, 0, ...
    !    IdxCC = 4, 0, 0, 0, 0, 0, ...
    !    IdxCR = 5, 0, 0, 0, 0, 0, ...
    !
    !  * NH4SH ������ȿ���˴�Ϳ����ʪ��������ź����
    !
    !    IdxNH3    = 2
    !    IdxH2S    = 3
    !    IdxNH4SHc = 6
    !    IdxNH4SHr = 7
    !
    !���Ѥ��ʤ���ʬ�ˤϥ�����������Ƥ���. 
    !
    
    !�⥸�塼���ɤ߹���
    !
    use dc_types,      only: STRING
    use dc_iounit,     only: FileOpen
    use dc_message,    only: MessageNotify
    use ChemData,      only: ChemData_OneSpcID, &!���ؼ�� ID
      &                      ChemData_MolWt      !ʬ����
    use gridset,       only: ncmax  !���ؼ�ο�
    use namelist_util, only: namelist_filename
    
    !���ۤη�����ػ�
    implicit none
      
    !�ѿ����
    character(20), allocatable   :: Symbol(:)   !�������
    integer                      :: SpcWetNum   !������ʬ�β��ؼ�ο�
    integer                      :: s, s1, s2
    integer                      :: n1, n2, n3
    integer                      :: unit
    integer                      :: num

    !-----------------------------------------------------------------
    ! NAMELIST �����������
    !
    NAMELIST /composition_nml/ SpcWetSymbol, SpcWetMolFr

    SpcWetSymbol = '' 
    SpcWetMolFr  = 0.0d0
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=composition_nml)
    close(unit)
      
    !----------------------------------------------------------
    ! ������ʬ�� ID ������
    !
    !������ʬ�θĿ��������
    SpcWetNum = count(SpcWetSymbol /= "")
    if (SpcWetNum /= ncmax) then 
      call MessageNotify( "E", "basicset: ", "SpcWetNum is not equal to ncmax." )
    end if
    
    !����γ������
    allocate(SpcWetID(SpcWetNum), Symbol(SpcWetNum), MolWtWet(SpcWetNum))

    !SpcWetSymbol ��ʸ���󤫤�, -Rain, -Cloud ���������Τ� Symbol �Ȥ����ݴ�
    do s = 1, SpcWetNum
      n1 = index(SpcWetSymbol(s), '-Cloud' )
      n2 = index(SpcWetSymbol(s), '-Rain' )
      n3 = max(n1, n2)
      if (n3 == 0) then
        Symbol(s) = SpcWetSymbol(s)
      else
        Symbol(s) = SpcWetSymbol(s)(1:n3-1)
      end if
    end do
    
    !���ؼ�� ID �����
    do s =1, SpcWetNum
      SpcWetID(s) = ChemData_OneSpcID( Symbol(s) )
    end do
    
    !ʬ���̤��ݴ�
    do s = 1, SpcWetNum
      MolWtWet(s) = ChemData_MolWt(SpcWetID(s))
    end do


    !-----------------------------------------------------------
    ! ��γ�ȵ��Τ� ID ���Ȥ���
    !
    !����, ��, ���Ȥ�ʬΥ����. 
    SelectCloud: do s = 1, ncmax
      
      !'-g' �Ȥ���ʸ���󤬴ޤޤ���ΤθĿ��������
      n1 = index(SpcWetSymbol(s), '-g' )
      if (n1 /= 0) then
        GasNum        = GasNum + 1
        IdxG(GasNum)   = s
      end if
      
      !'Cloud' �Ȥ���ʸ���󤬴ޤޤ���ΤθĿ��������
      n2 = index(SpcWetSymbol(s), '-Cloud' )
      if (n2 /= 0) then
        CloudNum         = CloudNum + 1
        IdxC(CloudNum)  = s
      end if

      !'Rain' �Ȥ���ʸ���󤬴ޤޤ���ΤθĿ��������
      n3 = index(SpcWetSymbol(s), '-Rain' )
      if (n3 /= 0) then
        RainNum         = RainNum + 1
        IdxR(RainNum)   = s
      end if

    end do SelectCloud


    !�ŷ�������Ф���, �����ȱ��Ȥ��Ф��������. 
    SelectCond: do s = 1, ncmax
      
      ! NH4SH ��¸�ߤ�����
      if ( trim(SpcWetSymbol(s)) == 'NH4SH-s-Cloud' ) then 
        RactNum           = 1
        cycle SelectCond
      end if
      
      !'Cloud' �Ȥ���ʸ���󤬴ޤޤ���ΤθĿ��������
      n2 = index(SpcWetSymbol(s), '-Cloud' )
      if (n2 /= 0) then
        CondNum          = CondNum  + 1
        IdxCC(CondNum)   = s

        do s1 = 1, ncmax
          if ( trim(SpcWetSymbol(s1)) == trim(SpcWetSymbol(s)(1:n2-3)//'-g') ) then 
            IdxCG(CondNum)   = s1
          end if
        end do
        
        do s2 = 1, ncmax
          if ( trim(SpcWetSymbol(s2)) == trim(SpcWetSymbol(s)(1:n2-1)//'-Rain') ) then 
            IdxCR(CondNum)   = s2
          end if
        end do
      end if
      
    end do SelectCond
    
    !-----------------------------------------------------------
    ! β�������˥���, ����ӥ����˥���β�����Ǥ� ID �����
    !
    do s = 1, ncmax
      if ( trim(SpcWetSymbol(s)) == 'NH4SH-s-Cloud' ) then 

        IdxNH4SHc = s

        do s1 = 1, ncmax
          if ( trim(SpcWetSymbol(s1)) == 'NH3-g' ) then 
            IdxNH3 = s1
          end if
        end do

        do s2 = 1, ncmax
          if ( trim(SpcWetSymbol(s2)) == 'H2S-g' ) then 
            IdxH2S = s2
          end if
        end do

      end if
      
      if ( trim(SpcWetSymbol(s)) == 'NH4SH-s-Rain' ) then 
        IdxNH4SHr = s
      end if

    end do
        
    !-----------------------------------------------------------
    ! ID ���Ȥ���
    !
    allocate(IDGas(CondNum), IDCloud(CondNum), IDRain(CondNum))
    do s = 1, CondNum
      IDGas(s)   = SpcWetID(IdxCG(s))
      IDCloud(s) = SpcWetID(IdxCC(s))
      IDRain(s)  = SpcWetID(IdxCR(s))
    end do

    !-----------------------------------------------------------
    ! ��ǧ
    !
    call MessageNotify( "M", &
         &  "composition_init","GasNum   = %d", i=(/GasNum/)   )
    call MessageNotify( "M", &
         & "composition_init", "CloudNum = %d", i=(/CloudNum/) )    
    call MessageNotify( "M", &
        & "composition_init", "RainNum  = %d", i=(/RainNum/)  ) 
    call MessageNotify( "M", &
         & "composition_init", "CondNum  = %d", i=(/CondNum/)  )    
    call MessageNotify( "M", &
         & "composition_init", "RactNum  = %d", i=(/RactNum/)  ) 
    call MessageNotify( "M", &
         & "composition_init", "IdxNH3 = %d",   i=(/IdxNH3/)   )
    call MessageNotify( "M", &
         & "composition_init", "IdxH2S = %d",   i=(/IdxH2S/)   )
    call MessageNotify( "M", &
         & "composition_init", "IdxNH4SHc = %d", i=(/IdxNH4SHc/) )
    call MessageNotify( "M", &
        & "composition_init", "IdxNH4SHr = %d", i=(/IdxNH4SHr/) )

    Num = count(IdxG /= 0)
    call MessageNotify( "M", &
      & "composition_init", "IdxG = %d %d %d %d %d", i=(/IdxG(1:Num)/) )
    call MessageNotify( "M", &
      & "composition_init", "IdxC = %d %d %d %d %d", i=(/IdxC(1:Num)/) )
    call MessageNotify( "M", &
      & "composition_init", "IdxR = %d %d %d %d %d", i=(/IdxR(1:Num)/) )
    Num = count(IdxCG /= 0)
    call MessageNotify( "M", &
      & "composition_init", "IdxCG = %d %d %d %d %d", i=(/IdxCG(1:Num)/) )
    call MessageNotify( "M", &
      & "composition_init", "IdxCC = %d %d %d %d %d", i=(/IdxCC(1:Num)/) )
    call MessageNotify( "M", &
      & "composition_init", "IdxCR = %d %d %d %d %d", i=(/IdxCR(1:Num)/) )

    do s = 1, SpcWetNum
       call MessageNotify( "M", &
            & "composition_init", "SpcWetID = %d",     i=(/SpcWetID(s)/) )
       call MessageNotify( "M", &
            & "composition_init", "SpcWetSymbol = %c", c1=trim(SpcWetSymbol(s)) )
       call MessageNotify( "M", &
            & "composition_init", "SpcWetMolFr = %f",  d=(/SpcWetMolFr(s)/) )
       call MessageNotify( "M", &
            & "composition_init", "MolWtWet = %f",     d=(/MolWtWet(s)/) )
    end do

  end subroutine composition_init

end module composition
