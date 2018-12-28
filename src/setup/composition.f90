!= 凝結成分に関する定数を決めるための変数参照型モジュール
!
! Authors::   SUGIYAMA Ko-ichiro
! Version::   $Id: composition.f90,v 1.7 2014/07/08 01:05:32 sugiyama Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
! License::   See COPYRIGHT[link:../../COPYRIGHT]

module composition
  !
  != 凝結成分に関する定数を決めるための変数参照型モジュール
  !
  
  !モジュール読み込み
  !
  use dc_types,      only: DP

  !暗黙の型宣言禁止
  !
  implicit none
  
  ! デフォルトの属性
  !
  private
  
  ! 変数定義
  !
  integer, save, public       :: GasNum       = 0   ! 気相の数
  integer, save, public       :: CloudNum     = 0   ! 雲の数  
  integer, save, public       :: RainNum      = 0   ! 雨の数
  integer, save, public       :: IdxG(10)     = 0   ! 気体の配列添え字
  integer, save, public       :: IdxC(10)     = 0   ! 雲の配列添え字
  integer, save, public       :: IdxR(10)     = 0   ! 雨の配列添え字  
  integer, save, public       :: CondNum      = 0   ! 凝結過程の数
  integer, save, public       :: IdxCG(10)    = 0   ! 凝結過程(気体)の配列添え字
  integer, save, public       :: IdxCC(10)    = 0   ! 凝結過程(雲)の配列添え字
  integer, save, public       :: IdxCR(10)    = 0   ! 凝結過程(雨)の配列添え字
  integer, save, public       :: RactNum      = 0   ! 化学反応の数
  integer, save, public       :: IdxNH3       = 0   ! NH3 (気体)の配列添え字
  integer, save, public       :: IdxH2S       = 0   ! H2S (気体)の配列添え字
  integer, save, public       :: IdxNH4SHc    = 0   ! NH4SH (雲)の配列添え字
  integer, save, public       :: IdxNH4SHr    = 0   ! NH4SH (雨)の配列添え字

  real(DP), save              :: SpcWetMolFr(20)    !湿潤成分の化学種の存在度
  character(20), save         :: SpcWetSymbol(20)   !湿潤成分の化学種名  

  integer, allocatable, save  :: SpcWetID(:)        !湿潤成分の化学種のID
  real(DP), allocatable, save :: MolWtWet(:)        !湿潤成分の分子量  

  integer, allocatable, save  :: IDGas(:)           !蒸気の化学種のID
  integer, allocatable, save  :: IDCloud(:)         !雲の化学種のID
  integer, allocatable, save  :: IDRain(:)          !雨の化学種のID

  ! サブルーチンの公開
  !
  public SpcWetID, MolWtWet, SpcWetMolFr, SpcWetSymbol
  public IDGas, IDCloud, IDRain
  public composition_init
  
contains
  
  subroutine composition_init
    !
    !=概要
    !
    !NameList ファイルから情報を取得する.
    !このサブルーチン内で, 化学情報の初期化を行っている
    !
    !=凝縮成分の取り扱いについて
    !
    !計算に利用する凝縮成分の情報は basicset.f90 で定義される
    !SpcWetSymbol と SpcWetID に保管されている
    ! 
    ! Symbol:  H2O-g, NH3-g, H2S-g, H2O-l-Cloud, H2O-l-Rain, NH4SH-s-Cloud, NH4SH-s-Rain
    ! ID:      5,     8,     10,    7,           7,          11,            11
    !
    !ID 番号(ChemData_SpcID)は ChemData.f90 で定義している
    !
    !上記の情報を元に, このルーチンでは以下の情報を作る.  
    !
    !  * 各カテゴリーに含まれる物質の数
    !
    !    GasNum = 3,  CloudNum = 2, RainNum = 2
    !  
    !  * 各カテゴリーの配列添え字. 気体だけに操作したい場合等々で利用する.
    !
    !    IdxG = 1, 2, 3, 0, 0, 0, ...
    !    IdxC = 4, 6, 0, 0, 0, 0, ...
    !    IdxR = 5, 7, 0, 0, 0, 0, ...
    !
    !  * 凝結(Condensation)を生じる物質の数と, それらの配列添え字. 
    !    上記の例では H2O の凝結のみが生じる
    !
    !    CondNum = 1
    !    IdxCG = 1, 0, 0, 0, 0, 0, ...
    !    IdxCC = 4, 0, 0, 0, 0, 0, ...
    !    IdxCR = 5, 0, 0, 0, 0, 0, ...
    !
    !  * NH4SH の生成反応に関与する物質の配列添え字
    !
    !    IdxNH3    = 2
    !    IdxH2S    = 3
    !    IdxNH4SHc = 6
    !    IdxNH4SHr = 7
    !
    !利用しない部分にはゼロを代入しておく. 
    !
    
    !モジュール読み込み
    !
    use dc_types,      only: STRING
    use dc_iounit,     only: FileOpen
    use dc_message,    only: MessageNotify
    use ChemData,      only: ChemData_OneSpcID, &!化学種の ID
      &                      ChemData_MolWt      !分子量
    use gridset,       only: ncmax  !化学種の数
    use namelist_util, only: namelist_filename
    
    !暗黙の型宣言禁止
    implicit none
      
    !変数定義
    character(20), allocatable   :: Symbol(:)   !作業配列
    integer                      :: SpcWetNum   !湿潤成分の化学種の数
    integer                      :: s, s1, s2
    integer                      :: n1, n2, n3
    integer                      :: unit
    integer                      :: num

    !-----------------------------------------------------------------
    ! NAMELIST から情報を取得
    !
    NAMELIST /composition_nml/ SpcWetSymbol, SpcWetMolFr

    SpcWetSymbol = '' 
    SpcWetMolFr  = 0.0d0
    
    call FileOpen(unit, file=namelist_filename, mode='r')
    read(unit, NML=composition_nml)
    close(unit)
      
    !----------------------------------------------------------
    ! 湿潤成分の ID を得る
    !
    !湿潤成分の個数を数える
    SpcWetNum = count(SpcWetSymbol /= "")
    if (SpcWetNum /= ncmax) then 
      call MessageNotify( "E", "basicset: ", "SpcWetNum is not equal to ncmax." )
    end if
    
    !配列の割り当て
    allocate(SpcWetID(SpcWetNum), Symbol(SpcWetNum), MolWtWet(SpcWetNum))

    !SpcWetSymbol の文字列から, -Rain, -Cloud を除いたものを Symbol として保管
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
    
    !化学種の ID を取得
    do s =1, SpcWetNum
      SpcWetID(s) = ChemData_OneSpcID( Symbol(s) )
    end do
    
    !分子量を保管
    do s = 1, SpcWetNum
      MolWtWet(s) = ChemData_MolWt(SpcWetID(s))
    end do


    !-----------------------------------------------------------
    ! 雲粒と気体の ID の組を作る
    !
    !蒸気, 雲, 雨とに分離する. 
    SelectCloud: do s = 1, ncmax
      
      !'-g' という文字列が含まれるものの個数を数える
      n1 = index(SpcWetSymbol(s), '-g' )
      if (n1 /= 0) then
        GasNum        = GasNum + 1
        IdxG(GasNum)   = s
      end if
      
      !'Cloud' という文字列が含まれるものの個数を数える
      n2 = index(SpcWetSymbol(s), '-Cloud' )
      if (n2 /= 0) then
        CloudNum         = CloudNum + 1
        IdxC(CloudNum)  = s
      end if

      !'Rain' という文字列が含まれるものの個数を数える
      n3 = index(SpcWetSymbol(s), '-Rain' )
      if (n3 /= 0) then
        RainNum         = RainNum + 1
        IdxR(RainNum)   = s
      end if

    end do SelectCloud


    !凝結過程に対して, 蒸気と雲との対を作成する. 
    SelectCond: do s = 1, ncmax
      
      ! NH4SH が存在する場合
      if ( trim(SpcWetSymbol(s)) == 'NH4SH-s-Cloud' ) then 
        RactNum           = 1
        cycle SelectCond
      end if
      
      !'Cloud' という文字列が含まれるものの個数を数える
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
    ! 硫化アンモニウム, およびアンモニアと硫化水素の ID を取得
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
    ! ID の組を作る
    !
    allocate(IDGas(CondNum), IDCloud(CondNum), IDRain(CondNum))
    do s = 1, CondNum
      IDGas(s)   = SpcWetID(IdxCG(s))
      IDCloud(s) = SpcWetID(IdxCC(s))
      IDRain(s)  = SpcWetID(IdxCR(s))
    end do

    !-----------------------------------------------------------
    ! 確認
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
