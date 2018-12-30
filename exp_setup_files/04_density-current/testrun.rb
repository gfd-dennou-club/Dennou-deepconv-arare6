#!/usr/bin/env ruby 
# -*- coding: utf-8 -*-
#
# テストラン実行用のスクリプト
# 

require "numru/ggraph"
include NumRu

# デフォルトファイル
#
prefix = "density-current"
ofile  = "#{prefix}.sample"

# 実行ファイル置き場
bindir = "../../bin"

# パラメタ情報
#
#dim   = [ [512, 1, 64], [512, 5, 64], [5, 512, 64], [1024, 1, 128], [256, 1, 32], [128, 1, 16]]
dim   = [ [512, 1, 64] ]
#reg   = [ ["51.2d3", "100.0d0", "6.4d3"], ["51.2d3", "500.0d0", "6.4d3"], ["500.0d0", "51.2d3", "6.4d3"], ["51.2d3", "100.0d0", "6.4d3"], ["51.2d3", "100.0d0", "6.4d3"], ["51.2d3", "100.0d0", "6.4d3"]]
reg   = [ ["51.2d3", "100.0d0", "6.4d3"] ]
ini   = [ "CosXZ", "CosXZ", "CosYZ", "CosXZ", "CosXZ", "CosXZ" ]
#flag1 = ["std", "Center4"]
#flag2 = ["std", "Center2"]
flag1 = ["Center4"]
flag2 = ["Center2"]

# 描画時刻
#
disptime = [0, 300, 600, 900]
#disptime = [900]

# 保管用の配列
#
prefixs = Array.new
comm = Array.new

###
### 設定ファイルの作成
###
dim.size.times{ |i|
  flag1.size.times{ |j|
    flag2.size.times{ |k|

      prefix2 = "#{prefix}_#{dim[i].join('x')}_#{flag1[j]}_#{flag2[k]}"
      cfile   = "#{prefix2}.conf" 
      prefixs.push( prefix2 )

      orig = open( ofile )      
      conf = open( cfile, "w" )
      while line = orig.gets
        line.chomp!
#        p line

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
        elsif (/FlagAdvection/ =~ line)
          conf.write( "  FlagAdvection = \"#{flag1[j]}\"\n" )
        elsif (/FlagAcousticmode/ =~ line)
          conf.write( "  FlagAcousticmode = \"#{flag2[k]}\"\n" )
        elsif (/InitialFile/ =~ line)
          conf.write( "  InitialFile  = \"#{prefix2}_ini.nc\" \n" )
        elsif (/OutputFile/ =~ line)
          conf.write( "  OutputFile   = \"#{prefix2}_res.nc\" \n" )
        elsif (/FilePrefix/ =~ line)
          conf.write( "  FilePrefix   = \"#{prefix2}_\" \n" )
        elsif (/FlagDisturbPTemp/ =~ line)
          conf.write( "  FlagDisturbPTemp = \"#{ini[i]}\" \n" )
        else
          conf.write( "#{line} \n" )
        end
      end
      conf.close
    }
  }
}

###
### 計算の実行
###
prefixs.each{ |prefix2|
  p "#{bindir}/arare_init-data -N=#{prefix2}.conf" 
  p "#{bindir}/arare -N=#{prefix2}.conf" 

  unless (FileTest.exist?("#{prefix2}_PTemp.nc"))
    system( "#{bindir}/arare_init-data -N=#{prefix2}.conf" )
    system( "#{bindir}/arare -N=#{prefix2}.conf" )
  end
}

###
### 描画設定
###
DCL.swpset('IWIDTH',  800)
DCL.swpset('IHEIGHT', 400)
#DCL.sgscmn(4)
#DCL.sgscmn(10)

DCL.swlset("lwnd",false)
DCL.gropn(4)
DCL.sgpset('lcntl', false)   # 制御文字を解釈しない
DCL.sgpset('lfull',true)     # 全画面表示
DCL.sgpset('lcorner',false)  # コーナーマークを書かない
DCL.uzfact(0.6)              # 座標軸の文字列サイズを定数倍
DCL.sgpset('lfprop',true)    # プロポーショナルフォントを使う
#DCL.udpset('lmsg',false)     # コンター間隔非表示
#DCL.sgpset('lclip',true)     # 枠からはみ出した分を描画しない.

#rmiss = DCL.glpget('rmiss')

###
### 描画
###
i = 0
prefixs.each{ |prefix2|
  disptime.each{ |t0|
    
    url = "#{prefix2}_PTemp.nc@PTemp"
    p prefix2
    
    if ( /_\d\d+x/ =~ prefix2 ) 
      ptemp = GPhys::IO.open_gturl(url).cut('y'=>0).cut('x'=>25.5e3..45e3) * 1.0
    else
      ptemp = GPhys::IO.open_gturl(url).cut('x'=>0).cut('y'=>25.5e3..45e3) * 1.0
    end
    
    GGraph.set_fig('viewport'=>[0.15, 0.8, 0.15, 0.4], 'new_frame'=>true)  
    GGraph.tone( ptemp.cut( 't'=>t0 ),
                 true, 
                 'annot'=>false,
                 'title'=>"Pot. Temp. Anomaly (t=#{t0})", 
                 'max'=> 0.0,
                 'min'=>-20.0,
                 )    
#    GGraph.color_bar('vwidth'=>0.01, 'tickintv'=>0, "landscape"=>true)
    GGraph.color_bar('vwidth'=>0.01, 'tickintv'=>0, 'voff'=>0.02)
    GGraph.contour( ptemp.cut( 't'=>t0 ),
                    false, 
                    'annot'=>false,
                    'int' => 0.5,
                    'nozero'=>true
                    )    
    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', "#{prefix2}", 1) 
    
    comm.push( "mv dcl_#{sprintf('%03d',i+1)}.png #{prefix2}_PTemp_#{sprintf('%03d',t0)}.png" )
    i = i + 1
  }
}
DCL.grcls

comm.size.times{|i|
  p "#{comm[i]}" 
  system( "#{comm[i]}" ) 
}
