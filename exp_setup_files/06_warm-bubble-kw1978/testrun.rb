#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
#
# テストラン実行用のスクリプト
# 

require "numru/ggraph"
include NumRu

# デフォルトファイル
#
prefix = "warm-bubble-kw1978"
ofile  = "#{prefix}.sample"

# 実行ファイル置き場
bindir = "../../bin"

# パラメタ情報
#
dim   = [
         [160, 1,   80], 
         [160, 5,   80], 
         [5,   160, 80],
         [40,  1,   20], 
        ]
reg   = [ 
         ["20.0d3",  "125.0d0", "10.0d3"], 
         ["20.0d3",  "625.0d0", "10.0d3"], 
         ["625.0d0", "20.0d3",  "10.0d3"],
         ["20.0d3",  "125.0d0", "10.0d3"]
        ]
ini   = [ "CosXZ", "CosXZ", "CosYZ", "CosXZ"]
flag1 = ["std", "Center4"]
flag2 = ["std", "Center2"]
flag3 = ["std", "Center2"]
vel0  = [0.0, 20.0]

# 描画時刻
#
disptime = [1020]

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
      vel0.size.times{ |l|
        flag2.size.times{ |m|
          
          prefix2 = "#{prefix}_#{dim[i].join('x')}_#{flag1[j]}_#{flag2[k]}_Tub-#{flag3[m]}_#{vel0[l]}"
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
            elsif (/FlagTurbulence/ =~ line)
              conf.write( "  FlagTurbulence = \"#{flag3[m]}\" \n" )
            elsif (/VelX0/ =~ line)
              if dim[i][0] > 10 
                conf.write( "  VelX0 = #{vel0[l]} \n" )
              else
                conf.write( "  VelX0 = 0.0d0 \n" )
              end
            elsif (/VelY0/ =~ line)
              if dim[i][1] > 10 
                conf.write( "  VelY0 = #{vel0[l]} \n" )
              else
                conf.write( "  VelY0 = 0.0d0 \n" )
              end
            else
              conf.write( "#{line} \n" )
            end
          end
          conf.close
        }
      }
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
DCL.swpset('IWIDTH',  1600)
DCL.swpset('IHEIGHT', 800)
DCL.sgscmn(4)
#DCL.sgscmn(10)

DCL.swlset("lwnd",false)
DCL.gropn(1)
DCL.sgpset('lcntl', false)   # 制御文字を解釈しない
DCL.sgpset('lfull',true)     # 全画面表示
DCL.sgpset('lcorner',false)  # コーナーマークを書かない
DCL.uzfact(0.4)              # 座標軸の文字列サイズを定数倍
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
    
    url1 = "#{prefix2}_PTemp.nc@PTemp"
    url2 = "#{prefix2}_VelZ.nc@VelZ"
    p prefix2

    next unless FileTest.exist?( "#{prefix2}_PTemp.nc" )
    
    if ( /_\d\d+x/ =~ prefix2 ) 
      ptemp = GPhys::IO.open_gturl(url1).cut('y'=>0) * 1.0
      velz  = GPhys::IO.open_gturl(url2).cut('y'=>0) * 1.0
    else
      ptemp = GPhys::IO.open_gturl(url1).cut('x'=>0) * 1.0
      velz  = GPhys::IO.open_gturl(url2).cut('x'=>0) * 1.0
    end
    
    GGraph.set_fig('viewport'=>[0.1, 0.45, 0.15, 0.4], 'new_frame'=>true)  
    GGraph.tone( ptemp.cut( 't'=>t0 ),
                 true, 
                 'annot'=>false,
                 'title'=>"theta (t=#{t0})", 
                 'max'=>  1.5,
                 'min'=> -0.1,
                 )    
#    GGraph.color_bar('vwidth'=>0.01, 'tickintv'=>0, "landscape"=>true)
    GGraph.color_bar('vwidth'=>0.01, 'tickintv'=>0, 'voff'=>0.02)
    GGraph.contour( ptemp.cut( 't'=>t0 ),
                    false, 
                    'annot'=>false,
                    'int' => 0.25,
                    'nozero'=>true
                    )    
#    DCL.uxsttl('t', "#{prefix2}", 1) 
    
    GGraph.set_fig('viewport'=>[0.55, 0.9, 0.15, 0.4], 'new_frame'=>false)  
    GGraph.tone( velz.cut( 't'=>t0 ),
                 true, 
                 'annot'=>false,
                 'title'=>"w (t=#{t0})", 
                 'max'=>  15.0,
                 'min'=> -15.0,
                 )    
#    GGraph.color_bar('vwidth'=>0.01, 'tickintv'=>0, "landscape"=>true)
    GGraph.color_bar('vwidth'=>0.01, 'tickintv'=>0, 'voff'=>0.02)
    GGraph.contour( velz.cut( 't'=>t0 ),
                    false, 
                    'annot'=>false,
                    'int' => 1.5,
                    'nozero'=>true
                    )    
    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', "#{prefix2}", 1) 
    
    comm.push( "mv dcl_#{sprintf('%04d',i+1)}.png #{prefix2}_#{sprintf('%03d',t0)}.png" )
    i = i + 1
  }
}
DCL.grcls

comm.size.times{|i|
  p "#{comm[i]}" 
  system( "#{comm[i]}" ) 
}
