#!/usr/bin/ruby -Ke
# -*- coding: utf-8 -*-
#
# テストラン実行用のスクリプト
# 

require "numru/ggraph"
include NumRu

# デフォルトファイル
#
prefix = "sound-wave"
ofile  = "#{prefix}.conf"

# 実行ファイル置き場
bindir = "../../bin"

# パラメタ情報
#
dim   = [ [50, 1, 50], [50, 5, 50], [50, 50, 5], [5, 50, 50] ]
reg   = [ ["10.0d3", "2.0d2", "10.0d3"],  ["10.0d3", "1.0d3", "10.0d3"], ["10.0d3", "10.0d3", "1.0d2"], ["10.0d2", "10.0d3", "10.0d3"] ]
type  = ["XZ", "XZ", "XY", "YZ"]
flag1 = ["std", "Center2"]

# 描画時刻
#
disptime = 100

# 保管用の配列
#
prefixs = Array.new
comm = Array.new

###
### 設定ファイルの作成
###
dim.size.times{ |i|
  flag1.size.times{ |j|
    
    prefix2 = "#{prefix}_#{type[i]}_#{dim[i].join('x')}_#{flag1[j]}"
    cfile   = "#{prefix2}.conf" 
    prefixs.push( prefix2 )
    
    orig = open( ofile )      
    conf = open( cfile, "w" )
    while line = orig.gets
      line.chomp!
      
      if (/xdim/ =~ line)
        conf.write( "  xdim = #{dim[i][0]} \n" )
      elsif (/ydim/ =~ line)
        conf.write( "  ydim = #{dim[i][1]} \n" )
      elsif (/zdim/ =~ line)
        conf.write( "  zdim = #{dim[i][2]} \n" )
      elsif (/Xmax/ =~ line)
        conf.write( "  Xmax = #{reg[i][0]} \n" )
      elsif (/Ymax/ =~ line)
        conf.write( "  Ymax = #{reg[i][1]} \n" )
      elsif (/Zmax/ =~ line)
        conf.write( "  Zmax = #{reg[i][2]} \n" )
      elsif (/FlagAcousticmode/ =~ line)
        conf.write( "  FlagAcousticmode = \"#{flag1[j]}\"\n" )
      elsif (/InitialFile/ =~ line)
        conf.write( "  InitialFile  = \"#{prefix2}_ini.nc\" \n" )
      elsif (/OutputFile/ =~ line)
        conf.write( "  OutputFile   = \"#{prefix2}_res.nc\" \n" )
      elsif (/FilePrefix/ =~ line)
        conf.write( "  FilePrefix   = \"#{prefix2}_\" \n" )
      elsif (/FlagDisturbExner/ =~ line) 
        conf.write( "  FlagDisturbExner = \"Gauss#{type[i]}\" \n")
      else
        conf.write( "#{line} \n" )
      end
      
    end
    conf.close
    orig.close
  }
}

###
### 計算の実行
###
prefixs.each{ |prefix2|
  p "#{bindir}/arare_init-data -N=#{prefix2}.conf" 
  p "#{bindir}/arare -N=#{prefix2}.conf" 
  
  system( "#{bindir}/arare_init-data -N=#{prefix2}.conf" )
  system( "#{bindir}/arare -N=#{prefix2}.conf" )
}

###
### 描画設定
###
DCL.swpset('IWIDTH',  1200)
DCL.swpset('IHEIGHT', 600)
DCL.sgscmn(4)
#DCL.sgscmn(10)

DCL.swlset("lwnd",false)
DCL.gropn(4)
DCL.sgpset('lcntl', false)   # 制御文字を解釈しない
DCL.sgpset('lfull',true)     # 全画面表示
DCL.sgpset('lcorner',false)  # コーナーマークを書かない
DCL.uzfact(0.5)              # 座標軸の文字列サイズを定数倍
DCL.sgpset('lfprop',true)    # プロポーショナルフォントを使う
#DCL.udpset('lmsg',false)     # コンター間隔非表示
#DCL.sgpset('lclip',true)     # 枠からはみ出した分を描画しない.

#rmiss = DCL.glpget('rmiss')

###
### 描画
###
i = 0
prefixs.each{ |prefix2|
    
  url = "#{prefix2}_Exner.nc@Exner"
  
  if    ( /_XY_/ =~ prefix2 ) 
    exner1 = GPhys::IO.open_gturl(url).cut('z'=>0).mean('x').cut('t'=>0..disptime) * 1.0
    exner2 = GPhys::IO.open_gturl(url).cut('z'=>0).mean('y').cut('t'=>0..disptime) * 1.0
    t1 = "Exner func. (x-mean)"
    t2 = "Exner func. (y-mean)"
  elsif ( /_XZ_/ =~ prefix2 ) 
    exner1 = GPhys::IO.open_gturl(url).cut('y'=>0).mean('x').cut('t'=>0..disptime) * 1.0
    exner2 = GPhys::IO.open_gturl(url).cut('y'=>0).mean('z').cut('t'=>0..disptime) * 1.0
    t1 = "Exner func. (x-mean)"
    t2 = "Exner func. (z-mean)"
  elsif ( /_YZ_/ =~ prefix2 ) 
    exner1 = GPhys::IO.open_gturl(url).cut('x'=>0).mean('y').cut('t'=>0..disptime) * 1.0
    exner2 = GPhys::IO.open_gturl(url).cut('x'=>0).mean('z').cut('t'=>0..disptime) * 1.0
    t1 = "Exner func. (y-mean)"
    t2 = "Exner func. (z-mean)"
  end
  
  GGraph.set_fig('viewport'=>[0.1, 0.45, 0.15, 0.4], 'new_frame'=>true)   
  GGraph.tone( exner1,
               true, 
               'annot'=>false,
               'title'=> t1,
               'max' => 5.0e-3,
               'min' => 0.0
               )    
  
  GGraph.set_fig('viewport'=>[0.55, 0.9, 0.15, 0.4], 'new_frame'=>false)  
  GGraph.tone( exner2,
               true, 
               'annot'=>false,
               'title'=> t2,
               'max' => 5.0e-3,
               'min' => 0.0
               )    
  GGraph.color_bar('vwidth'=>0.01, 'tickintv'=>0, 'voff'=>0.02)

  DCL.uxsttl('b', ".", 1) 
  DCL.uxsttl('b', ".", 1) 
  DCL.uxsttl('b', ".", 1) 
  DCL.uxsttl('b', "#{prefix2}", 1) 
  
  comm.push( "mv dcl_#{sprintf('%03d',i+1)}.png #{prefix2}_Exner.png" )
  i = i + 1

}
DCL.grcls

comm.size.times{|i|
  p "#{comm[i]}" 
  system( "#{comm[i]}" ) 
}
