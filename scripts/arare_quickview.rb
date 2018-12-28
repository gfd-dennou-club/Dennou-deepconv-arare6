#!/usr/bin/ruby -Ke
# -*- coding: utf-8 -*-
#
# おまかせ作図用スクリプト
# 
require 'yaml'
require 'pp'
require "numru/ggraph"
include NumRu
require 'optparse'

# オプション解析
params = ARGV.getopts('y:','d:')
dir   = params["d"]
yfile = params["y"]

###
### データ置き場が指定されていない場合は help を出す.
###
unless dir
  print "Please set options \n"
  print "ruby arare_quickview.rb -d <data dir> -y <YAML file> \n"
  exit
end

###
### YAML ファイルが指定されていない場合は雛形を作る.
###
unless yfile

  # 標準出力
  print "Input no YAML file....\n"
  print "Now making sample YAML file....\n"

  # ファイル出力
  sfile = open("sample.yml", 'w')
  sfile.print( "Vars: \n" )

  # netCDF ファイル名を取り出す
  vars = Array.new
  nfiles = Dir.glob("#{dir}*.nc").sort
  
  # netcdf を開いて変数を取り出す
  GDir.top='/'
  
  # 初期化
  t = Array.new

  nfiles.size.times{|i|
    nfile = nfiles[i]
    next if /restart/ =~ nfile
    next if /Analysis/ =~ nfile
    next if /BasicZ/ =~ nfile
    next if /CAPE/ =~ nfile
    next if /DownFlow/ =~ nfile
    next if /Energy/ =~ nfile

    # ファイル出力
    sfile.print( " #{nfile}: \n" )

    print nfile,":\n"
    gdir = GDir.new(File.expand_path(nfile))
    
    # 変数名を保管
    vars[i] = Array.new
    gdir.list_data.each{|s| 
      vars[i].push( s )
    }
    
    # 変数を開いて最大値・最小値を取り出す. 
    #
    axis = ""
    vars[i].each{|var|
      
      # 軸の場合は名称を保管しておく.
      next if var == "x" || var == "y" || var == "z" || var == "s" || var == "t"
      
      # ファイルのオープン
      gphys = GPhys::IO.open( nfile, var )

      # 軸情報の取り出し
      # (netCDF 内で定義されている軸が必ずしも使われているわけでない)
      rank = gphys.rank
      axis = ""
      rank.times{|j|
        a0 = gphys.coord(j).name
        axis = "#{axis}#{a0}"
        if a0 == "t"
          t.push( gphys.coord(j).to_a )
        end
      }
      p "#{var} (#{axis})"
      
      # ファイル出力
      sfile.print( "  #{var}: \n" )
      sfile.print( "   axis: #{axis} \n" )

      if /xyzt/ =~ axis || /xzt/ =~ axis || /xyt/ =~ axis || /xt/ =~ axis 
        # スナップショット用
        min1 = gphys.min.to_f.abs
        max1 = gphys.max.to_f.abs

        # 空間平均用
        if /xyzt/ =~ axis || /xyt/ =~ axis
          min2 = gphys.mean('x').mean('y').min.to_f.abs
          max2 = gphys.mean('x').mean('y').max.to_f.abs
        elsif /xzt/ =~ axis || /xt/ =~ axis 
          min2 = gphys.mean('x').min.to_f.abs
          max2 = gphys.mean('x').max.to_f.abs
        end
        
        if min1 > max1 
          min_snp = min1 * -1.0
          max_snp = min1
        else
          min_snp = max1 * -1.0
          max_snp = max1
        end

        if min2 > max2
          min_mean = min2 * -1.0
          max_mean = min2
        else
          min_mean = max2 * -1.0
          max_mean = max2
        end       

        if /Cloud$/ =~ var
          sfile.print( "   min_snp: 1.0e-8 \n" )
          sfile.print( "   max_snp: #{sprintf('%.e', max_snp)}\n" )
          sfile.print( "   min_mean: 1.0e-8 \n" )
          sfile.print( "   max_mean: #{sprintf('%.e', max_mean)}\n" )
          sfile.print( "   itr: 3 \n" ) 
        elsif /Rain$/ =~ var
          sfile.print( "   min_snp: 1.0e-8 \n" )
          sfile.print( "   max_snp: #{sprintf('%.e', max_snp)}\n" )
          sfile.print( "   min_mean: 1.0e-8 \n" )
          sfile.print( "   max_mean: #{sprintf('%.e', max_mean)}\n" )
          sfile.print( "   itr: 3 \n" )
        else
          sfile.print( "   min_snp: #{sprintf('%.e', min_snp)}\n" )
          sfile.print( "   max_snp: #{sprintf('%.e', max_snp)}\n" )
          sfile.print( "   min_mean: #{sprintf('%.e', min_mean)}\n" )
          sfile.print( "   max_mean: #{sprintf('%.e', max_mean)}\n" )
          sfile.print( "   itr: 1 \n" )
        end

      elsif /^zt/ =~ axis || /^t/ =~ axis 
        # x, y 方向が含まれない場合

        min1 = gphys.min.to_f.abs
        max1 = gphys.max.to_f.abs

        if min1 > max1 
          min_mean = min1 * -1.0
          max_mean = min1
        else
          min_mean = max1 * -1.0
          max_mean = max1
        end
        sfile.print( "   min_mean: #{sprintf('%.e', min_mean)}\n" )
        sfile.print( "   max_mean: #{sprintf('%.e', max_mean)}\n" )
        sfile.print( "   itr: 1 \n" )
      end
    }
  }
  # ファイル出力
  sfile.print( "Time: \n" )
  sfile.print( " min: #{t[0].min} \n" )
  sfile.print( " max: #{t[0].max} \n" )
  del0 = t[0][2] - t[0][1]
  fct = (t[0].max - t[0].min) / del0 / 20
  sfile.print( " del: #{del0 * fct.to_i} \n" )

  sfile.print( "Axis: \n" )
  sfile.print( " x: \n" )
  sfile.print( " y: [0.0] \n" )  # 決めうち
  sfile.print( " z:  \n" )
  exit
end


###
### yml ファイルの指定に従って作図する
###

# yml ファイルの読み込み
data  = YAML.load_file( yfile )
#pp data

###
### 描画設定
###
DCL.swpset('IWIDTH',  800)
DCL.swpset('IHEIGHT', 800)
#DCL.sgscmn(4)
DCL.sgscmn(10)

DCL.swlset("lwnd",false)
DCL.gropn(4)
DCL.sgpset('lcntl', false)   # 制御文字を解釈しない
DCL.sgpset('lfull',true)     # 全画面表示
DCL.sgpset('lcorner',false)  # コーナーマークを書かない
DCL.uzfact(0.7)              # 座標軸の文字列サイズを定数倍
DCL.sgpset('lfprop',true)    # プロポーショナルフォントを使う
#DCL.udpset('lmsg',false)     # コンター間隔非表示
#DCL.sgpset('lclip',true)     # 枠からはみ出した分を描画しない.

#rmiss = DCL.glpget('rmiss')

###
### 描画
###
comm = Array.new

x = data["Axis"]["x"]
y = data["Axis"]["y"]
z = data["Axis"]["z"]

t0 = data["Time"]["min"]
t1 = data["Time"]["max"]
dt = data["Time"]["del"]
num= ((t1 - t0) / dt).to_i
i = 0

data["Vars"].keys.each{|nfile|
  p nfile
  data["Vars"]["#{nfile}"].keys.each{|var|
    p var 

    # ファイルのオープン
    gphys = GPhys::IO.open( nfile, var )
    
    itr  = data["Vars"]["#{nfile}"]["#{var}"]["itr"].to_i
    axis = data["Vars"]["#{nfile}"]["#{var}"]["axis"]
    mmax = data["Vars"]["#{nfile}"]["#{var}"]["max_mean"].to_f + 1.0e-8
    mmin = data["Vars"]["#{nfile}"]["#{var}"]["min_mean"].to_f - 1.0e-8
    smax = data["Vars"]["#{nfile}"]["#{var}"]["max_snp"].to_f + 1.0e-8
    smin = data["Vars"]["#{nfile}"]["#{var}"]["min_snp"].to_f - 1.0e-8  
    
    if itr == 3
      gphys = ( gphys + 1.0e-8 ).log10
      smin = -8
      mmin = -8
      smax = -2
      mmax = -2
    end

    ## 
    ## スナップショット
    ## 
    num.times{|j|
      tn = t0 + dt * j
      
      y.each{|y0|
        GGraph.set_fig('viewport'=>[0.15, 0.8, 0.3, 0.8], 'new_frame'=>true)  
        
        if axis == "xyzt"
          gphys0 = gphys.cut( 'y'=> y0 ).cut( 't'=> tn )
        elsif axis == "xzt" || axis == "xt"
          gphys0 = gphys.cut( 't'=> tn )
        elsif axis == "xyt"
          gphys0 = gphys.cut( 'y'=> y0 ).cut( 't'=> tn )
        else
          break
        end
        
        if gphys0.rank == 1
          GGraph.line( gphys0, true, 'annot'=>false, 'max'=> smax, 'min'=> smin)
        else
          GGraph.tone( gphys0, true, 'annot'=>false, 'max'=> smax, 'min'=> smin)
          GGraph.color_bar('tickintv'=>0, 'voff'=>0.02)
        end
        DCL.uxsttl('t', "t=#{tn} y=#{y0}", 1) 
        DCL.uxsttl('b', ".", 1) 
        DCL.uxsttl('b', ".", 1) 
        DCL.uxsttl('b', ".", 1) 
        DCL.uxsttl('b', "#{nfile}", 1)   
        comm.push("mv dcl_#{sprintf('%03d',i+1)}.png #{var}-snap_y#{sprintf('%08d',y0)}_t#{sprintf('%08d',tn)}.png")
        i = i + 1
      }
      
      
      if axis == "xyzt" && gphys.coord(1).to_a.size > 1
        z.each{|z0|
          GGraph.set_fig('viewport'=>[0.2, 0.8, 0.2, 0.8], 'new_frame'=>true)  
          GGraph.tone( gphys.cut( 'z'=> z0 ).cut( 't'=> tn ),
                       true, 
                       'annot'=>false,
                       'max'=> smax,
                       'min'=> smin
                       )          
          GGraph.color_bar('tickintv'=>0, 'voff'=>0.02)
          DCL.uxsttl('t', "t=#{tn} z=#{z0}", 1) 
          DCL.uxsttl('b', ".", 1) 
          DCL.uxsttl('b', ".", 1) 
          DCL.uxsttl('b', "#{nfile}", 1)   
          comm.push("mv dcl_#{sprintf('%03d',i+1)}.png #{var}-snap_z#{sprintf('%08d',z0)}_t#{sprintf('%08d',tn)}.png")
          i = i + 1
        }
      end
    }
    
    ##
    ## 時系列
    ##
    GGraph.set_fig('viewport'=>[0.15, 0.8, 0.3, 0.8], 'new_frame'=>true)  

    if axis == "xyzt" || axis == "xyt" 
      gphys0 = gphys.mean('x').mean('y')
    elsif axis == "xzt" || axis == "xt"
      gphys0 = gphys.mean('x')
    else
      gphys0 = gphys
    end
    
    var0 = var
    if axis == "t" || axis == "zt"
      var0 = "ZZ_#{var}"
    end

    if gphys0.rank == 1
      GGraph.line( gphys0, true, 'annot'=>false, 'max'=> mmax, 'min'=> mmin)      
    else
      GGraph.tone( gphys0, true, 'annot'=>false, 'max'=> mmax, 'min'=> mmin, 'exchange'=>true )
      GGraph.color_bar('tickintv'=>0, 'voff'=>0.02)
    end

    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', ".", 1) 
    DCL.uxsttl('b', "#{nfile}", 1)   
    comm.push("mv dcl_#{sprintf('%03d',i+1)}.png #{var0}-mean.png")    
    i = i + 1
  }  
}
DCL.grcls

comm.size.times{|i|
  p "#{comm[i]}" 
  system( "#{comm[i]}" ) 
}

