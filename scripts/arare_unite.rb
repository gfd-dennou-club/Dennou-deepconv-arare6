#!/usr/bin/ruby -Ke
# -*- coding: euc-jp -*-
#
# MPI でコア毎に出力されたファイルをひとまとめにするためのスクリプト. 
# 設定ファイルより, MPI のコア数とファイル名の接頭詞を読み込む. 
# なお, 出力先は, 時刻名のディレクトリである. 
#
# USAGE:
#
#   $ ruby arare_unite.rb <設定ファイル> <分割数>
#
# 想定するディレクトリ構造
#   ./arare.conf         設定ファイル
#   ./data/              生データ置き場
#          VelX_rank000000.nc
#          VelX_rank000001.nc
#          VelX_rank000002.nc
#          ....
#   ./TIME_XXXX-XXXXX/   結合したファイル置き場. 
#          VelX.nc


require "numru/ggraph"
include NumRu

dir     = ""
maxsize = 512 * 512 * 250 * 8

###
### 設定ファイルを読み込む. ファイル名の接頭詞と並列数を取り出す.
###
xprocs = 1
yprocs = 1

conf = open( ARGV[0] )
while line = conf.gets
  line.chomp!
  if (/FilePrefix/ =~ line)
    prefix = line.split("=")[1]
    prefix = prefix.gsub("\"","").gsub(" ","")
    prefix0 = File.basename(prefix)
  elsif (/xsub/ =~ line)
    xprocs = line.split("=")[1].to_i
  elsif (/ysub/ =~ line) 
    yprocs = line.split("=")[1].to_i
  end
end
conf.close


###
### 変数毎にファイルを取り出してまとめる. 1 変数 1 ファイルを仮定. 
###
Dir.glob("#{prefix}*rank000000.nc").sort.each{|nc|
  
  # ファイル名から変数名を判断
  tmp = nc.gsub("#{prefix}", "")
  /^([A-Z].+)_rank000000.nc$/ =~ tmp 

  # 変数に対する処理.
  if $1
    var = $1.to_s
    p var

    ####
    #### rank 0 のファイルをオープンして座標の情報を取り出した後に, 
    #### ファイル名の配列を作り, オープンし, 出力する. 
    ####
    gphys0 = GPhys::IO.open( nc, var )  
    tmp2 = gphys0.coord(0).name

    if tmp2 == "t"

      # 平均操作を行う
      files = Dir.glob( "#{prefix}#{var}_rank*.nc" )
      gphys = GPhys::IO.open( files[0], var )  / files.size.to_f
      (files.size - 1).times{|kk|
        gphys = gphys + GPhys::IO.open( files[kk+1], var ) / files.size.to_f
      }

      # 時刻情報の取り出し
      t = gphys.coord(0).val      

      # ディレクトリ作成
      dir = "TIME_#{sprintf('%08d',t.min)}-#{sprintf('%08d',t.max)}"
      Dir::mkdir( dir ) unless FileTest.exist?( dir )
      
      # ファイルの書き出し. 
      newfile = "#{dir}/#{prefix0}#{var}.nc"
      p newfile
      unless FileTest.exist?( newfile )
        outfile = NetCDF.create(newfile)
        GPhys::NetCDF_IO.write( outfile, gphys ) 
        outfile.close
        p "WRITE: #{newfile}"
      end

    else

      # ファイル名の配列を作成. 
      files = Array.new
      yprocs.times{|j|
        files[j] = Array.new      
        xprocs.times{|i|
          k = xprocs * yprocs - (j + 1) * xprocs + i
          files[j].push( "#{prefix}#{var}_rank#{sprintf("%06d", k)}.nc")
        }
      }
      
      # ファイルのオープン
      if yprocs == 1 
        gphys = GPhys::IO.open( files[0], var )
      else
        gphys = GPhys::IO.open( files, var )
      end

      # 初期化
      rank = gphys.rank
      xyzsize = 1
      t = gphys.coord(0).val

      # 軸情報
      rank.times{|j|
        a0 = gphys.coord(j).name
        if a0 == "t"
          t = gphys.coord(j).val
          break 
        else
          xyzsize = xyzsize * gphys.coord(j).to_a.size
        end
      }

      # 1 ファイルに入れる時間方向のグリッド数を決める
      tsize = maxsize / xyzsize 
      taxis = t.size
      div   = taxis / tsize  + 1
      p "div: #{div}; n: #{taxis / div}"

      # 2 GB の壁問題があるので, 分割出力する
      #
      if div > 1
        
        n = taxis / div  # 1 ファイルあたりのデータ数
        
        # 割り切れないときのための処理. 
        if ( n % div > 0 )           
          n = n + 1 
        end

        div.times{|pp|
          n1 = n * pp
          n2 = [n * ( pp + 1 ) - 1, t.size - 1].min
#          p "#{pp}, #{n1} -> #{n2}"

          # ディレクトリ作成
          dir = "TIME_#{sprintf('%08d',t[n1])}-#{sprintf('%08d',t[n2])}"
          Dir::mkdir( dir ) unless FileTest.exist?( dir )

          # ファイル書き出し
          newfile = "#{dir}/#{prefix0}#{var}.nc"
          unless FileTest.exist?( newfile )
            outfile = NetCDF.create(newfile)
            GPhys::NetCDF_IO.write( outfile, gphys[true,true,true,n1..n2] ) 
            outfile.close
            p "time : #{t[n1]}--#{t[n2]}"
            p "WRITE: #{newfile}"
          end
        }

      else

        # ディレクトリ作成
        dir = "TIME_#{sprintf('%08d',t.min)}-#{sprintf('%08d',t.max)}"
        Dir::mkdir( dir ) unless FileTest.exist?( dir )

        # ファイルの書き出し. 
        newfile = "#{dir}/#{prefix0}#{var}.nc"
        unless FileTest.exist?( newfile )
          outfile = NetCDF.create(newfile)
          GPhys::NetCDF_IO.write( outfile, gphys ) 
          outfile.close
          p "WRITE: #{newfile}"
        end
      end

    end
  else
    p "SKIP: #{tmp}"
  end  
}

