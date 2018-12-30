#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
#
# テストラン実行用のスクリプト
# 

require "numru/ggraph"
include NumRu

# デフォルトファイル
#
prefix = "advection"
ofile  = "#{prefix}.sample"

# 実行ファイル置き場
bindir = "../../bin"

# パラメタ情報
#
dim   = [ [200, 1, 5], [200, 5, 5], [5, 200, 5]]
reg   = [ ["400.0d3", "2.0d3", "10.0d3"],  ["400.0d3", "10.0d3", "10.0d3"],  ["10.0d3", "400.0d3", "10.0d3"] ]
alpha = [ "0.0d0", "1.0d-3" ]
flag1 = ["std", "Center4"]
flag2 = ["std", "Center2"]
vel0  = [["20.0d0", "0.0d0"], ["20.0d0", "0.0d0"], ["0.0d0", "20.0d0"]]
ratio = [32, 20, 16, 12, 8, 4, 2]

# 描画時刻
#
disptime = [0,40000]

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
      alpha.size.times{ |l| 
        ratio.size.times{ |m| 
        
          prefix2 = "#{prefix}_#{dim[i].join('x')}_#{flag1[j]}_#{flag2[k]}_#{alpha[l]}_#{sprintf('%02d',ratio[m])}Del"
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
            elsif (/AlphaNDiff/ =~ line)
            conf.write( "  AlphaNDiff = #{alpha[l]} \n" )
            elsif (/VelX0/ =~ line)
              conf.write( "  VelX0 = #{vel0[i][0]} \n" )
            elsif (/VelY0/ =~ line)
              conf.write( "  VelY0 = #{vel0[i][1]} \n" )
            elsif (/XposMin/ =~ line)
              if dim[i][0] > dim[i][1] 
                conf.write( "  XposMin = #{101.0 - 1.0 * ratio[m]}d3 \n" )
              else
                conf.write( "  XposMin = 0.0d0 \n" )
              end
            elsif (/XposMax/ =~ line)
              if dim[i][0] > dim[i][1] 
                conf.write( "  XposMax = #{101.0 + 1.0 * ratio[m]}d3 \n" )
              else
                conf.write( "  XposMax = 10.0d3 \n" )
              end
            elsif (/YposMin/ =~ line)
              if dim[i][0] > dim[i][1] 
                conf.write( "  YposMin = 0.0d0 \n" )
              else
                conf.write( "  YposMin = #{101.0 - 1.0 * ratio[m]}d3 \n" )
              end
            elsif (/YposMax/ =~ line)
              if dim[i][0] > dim[i][1] 
                conf.write( "  YposMax = 10.0d3 \n" )
              else
                conf.write( "  YposMax = #{101.0 + 1.0 * ratio[m]}d3 \n" )
              end
            else
              conf.write( "#{line} \n" )
            end
          end
          conf.close
          orig.close
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
   
  unless (FileTest.exist?("#{prefix2}_H2O-g.nc"))
    system( "#{bindir}/arare_init-data -N=#{prefix2}.conf" )
    system( "#{bindir}/arare -N=#{prefix2}.conf" )
  end
}

###
### 描画設定
###
DCL.swpset('IWIDTH',  800)
DCL.swpset('IHEIGHT', 800)
DCL.sgscmn(4)
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
    
  url = "#{prefix2}_H2O-g.nc@H2O-g"
  p prefix2
    
  if ( /_\d\d+x/ =~ prefix2 ) 
    ptemp = GPhys::IO.open_gturl(url).cut('y'=>0) * 1.0
  else
    ptemp = GPhys::IO.open_gturl(url).cut('x'=>0) * 1.0
  end
  
  GGraph.set_fig('viewport'=>[0.15, 0.8, 0.15, 0.8], 'new_frame'=>true)  
  GGraph.line( ptemp.cut( 't'=>0 ),
               true, 
               'annot'=>false,
               'title'=>"Mixing ratio (t = 0, 4e4 sec)", 
               'max'=> 1.5e-3,
               'min'=>-5.0e-4,
               )    
  
  GGraph.line( ptemp.cut( 't'=>4e4 ),
               false, 
               'annot'=>false
               )    
  DCL.uxsttl('b', ".", 1) 
  DCL.uxsttl('b', "#{prefix2}", 1) 
  
  comm.push( "mv dcl_#{sprintf('%03d',i+1)}.png #{prefix2}_H2O-g.png" )
  i = i + 1

}
DCL.grcls

comm.size.times{|i|
  p "#{comm[i]}" 
  system( "#{comm[i]}" ) 
}
