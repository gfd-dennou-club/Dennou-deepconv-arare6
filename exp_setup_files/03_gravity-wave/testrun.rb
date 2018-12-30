#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
#
# テストラン実行用のスクリプト
# 

require "numru/ggraph"
include NumRu

# デフォルトファイル
#
prefix = "gravity-wave"
ofile  = "#{prefix}.sample"

# 実行ファイル置き場
bindir = "../../bin"

# パラメタ情報
#
dim   = [ [300, 1, 10], [300, 5, 10], [5, 300, 10]]
reg   = [ ["3.0d5", "1.0d3", "1.0d4"],  ["3.0d5", "5.0d3", "1.0d4"],  ["5.0d3", "3.0d5", "1.0d4"] ]
ini   = [ "SK1994XZ", "SK1994XZ", "SK1994YZ"]
#flag1 = ["std", "Center4"]
flag1 = ["Center4"]
#flag1 = ["std"]
#flag2 = ["std", "Center2"]
#flag2 = ["std"]
flag2 = ["Center2"]
vel0  = [["20.0d0", "0.0d0"], ["20.0d0", "0.0d0"], ["0.0d0", "20.0d0"]]

# 描画時刻
#
disptime = [3000]

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
        elsif (/VelX0/ =~ line)
          conf.write( "  VelX0 = #{vel0[i][0]} \n" )
        elsif (/VelY0/ =~ line)
          conf.write( "  VelY0 = #{vel0[i][1]} \n" )
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
DCL.sgscmn(4)
#DCL.sgscmn(10)

DCL.swlset("lwnd",false)
DCL.gropn(1)
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
      ptemp = GPhys::IO.open_gturl(url).cut('y'=>0) * 1.0
    else
      ptemp = GPhys::IO.open_gturl(url).cut('x'=>0) * 1.0
    end
    
    GGraph.set_fig('viewport'=>[0.15, 0.8, 0.15, 0.4], 'new_frame'=>true)  
    GGraph.tone( ptemp.cut( 't'=>t0 ),
                 true, 
                 'annot'=>false,
                 'title'=>"Pot. Temp. Anomaly (t=#{t0})", 
                 'max'=> 0.005,
                 'min'=>-0.005,
                 )    
#    GGraph.color_bar('vwidth'=>0.01, 'tickintv'=>0, "landscape"=>true)
    GGraph.color_bar('vwidth'=>0.01, 'tickintv'=>0, 'voff'=>0.02)
    GGraph.contour( ptemp.cut( 't'=>t0 ),
                    false, 
                    'annot'=>false,
                    'int' => 5.0e-4,
                    'nozero'=>true
                    )    
    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', "#{prefix2}", 1) 
    
    comm.push( "mv dcl_#{sprintf('%04d',i+1)}.png #{prefix2}_PTemp_#{sprintf('%03d',t0)}.png" )
    i = i + 1
  }
}
DCL.grcls

comm.size.times{|i|
  p "#{comm[i]}" 
  system( "#{comm[i]}" ) 
}
