#!/usr/bin/ruby -Ke
# -*- coding: euc-jp -*-
#
# MPI �ǥ�����˽��Ϥ��줿�ե������ҤȤޤȤ�ˤ��뤿��Υ�����ץ�. 
# ����ե�������, MPI �Υ������ȥե�����̾����Ƭ����ɤ߹���. 
# �ʤ�, �������, ����̾�Υǥ��쥯�ȥ�Ǥ���. 
#
# USAGE:
#
#   $ ruby arare_unite.rb <����ե�����> <ʬ���>
#
# ���ꤹ��ǥ��쥯�ȥ깽¤
#   ./arare.conf         ����ե�����
#   ./data/              ���ǡ����֤���
#          VelX_rank000000.nc
#          VelX_rank000001.nc
#          VelX_rank000002.nc
#          ....
#   ./TIME_XXXX-XXXXX/   ��礷���ե������֤���. 
#          VelX.nc


require "numru/ggraph"
include NumRu

dir     = ""
maxsize = 512 * 512 * 250 * 8

###
### ����ե�������ɤ߹���. �ե�����̾����Ƭ������������Ф�.
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
### �ѿ���˥ե��������Ф��ƤޤȤ��. 1 �ѿ� 1 �ե��������. 
###
Dir.glob("#{prefix}*rank000000.nc").sort.each{|nc|
  
  # �ե�����̾�����ѿ�̾��Ƚ��
  tmp = nc.gsub("#{prefix}", "")
  /^([A-Z].+)_rank000000.nc$/ =~ tmp 

  # �ѿ����Ф������.
  if $1
    var = $1.to_s
    p var

    ####
    #### rank 0 �Υե�����򥪡��ץ󤷤ƺ�ɸ�ξ������Ф������, 
    #### �ե�����̾���������, �����ץ�, ���Ϥ���. 
    ####
    gphys0 = GPhys::IO.open( nc, var )  
    tmp2 = gphys0.coord(0).name

    if tmp2 == "t"

      # ʿ������Ԥ�
      files = Dir.glob( "#{prefix}#{var}_rank*.nc" )
      gphys = GPhys::IO.open( files[0], var )  / files.size.to_f
      (files.size - 1).times{|kk|
        gphys = gphys + GPhys::IO.open( files[kk+1], var ) / files.size.to_f
      }

      # �������μ��Ф�
      t = gphys.coord(0).val      

      # �ǥ��쥯�ȥ����
      dir = "TIME_#{sprintf('%08d',t.min)}-#{sprintf('%08d',t.max)}"
      Dir::mkdir( dir ) unless FileTest.exist?( dir )
      
      # �ե�����ν񤭽Ф�. 
      newfile = "#{dir}/#{prefix0}#{var}.nc"
      p newfile
      unless FileTest.exist?( newfile )
        outfile = NetCDF.create(newfile)
        GPhys::NetCDF_IO.write( outfile, gphys ) 
        outfile.close
        p "WRITE: #{newfile}"
      end

    else

      # �ե�����̾����������. 
      files = Array.new
      yprocs.times{|j|
        files[j] = Array.new      
        xprocs.times{|i|
          k = xprocs * yprocs - (j + 1) * xprocs + i
          files[j].push( "#{prefix}#{var}_rank#{sprintf("%06d", k)}.nc")
        }
      }
      
      # �ե�����Υ����ץ�
      if yprocs == 1 
        gphys = GPhys::IO.open( files[0], var )
      else
        gphys = GPhys::IO.open( files, var )
      end

      # �����
      rank = gphys.rank
      xyzsize = 1
      t = gphys.coord(0).val

      # ������
      rank.times{|j|
        a0 = gphys.coord(j).name
        if a0 == "t"
          t = gphys.coord(j).val
          break 
        else
          xyzsize = xyzsize * gphys.coord(j).to_a.size
        end
      }

      # 1 �ե�������������������Υ���åɿ������
      tsize = maxsize / xyzsize 
      taxis = t.size
      div   = taxis / tsize  + 1
      p "div: #{div}; n: #{taxis / div}"

      # 2 GB �������꤬����Τ�, ʬ����Ϥ���
      #
      if div > 1
        
        n = taxis / div  # 1 �ե����뤢����Υǡ�����
        
        # ����ڤ�ʤ��Ȥ��Τ���ν���. 
        if ( n % div > 0 )           
          n = n + 1 
        end

        div.times{|pp|
          n1 = n * pp
          n2 = [n * ( pp + 1 ) - 1, t.size - 1].min
#          p "#{pp}, #{n1} -> #{n2}"

          # �ǥ��쥯�ȥ����
          dir = "TIME_#{sprintf('%08d',t[n1])}-#{sprintf('%08d',t[n2])}"
          Dir::mkdir( dir ) unless FileTest.exist?( dir )

          # �ե�����񤭽Ф�
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

        # �ǥ��쥯�ȥ����
        dir = "TIME_#{sprintf('%08d',t.min)}-#{sprintf('%08d',t.max)}"
        Dir::mkdir( dir ) unless FileTest.exist?( dir )

        # �ե�����ν񤭽Ф�. 
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

